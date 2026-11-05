import numpy as np
import polars as pl
import math
import os
from polars import col
import logging as log
import itertools
import copy
from itertools import count
from itertools import takewhile
from scipy.stats import entropy
import re
import multiprocessing
from multiprocessing import get_context
import progressbar
import nanomotif as nm
from nanomotif.constants import MOD_TYPE_TO_CANONICAL, MOD_CODE_TO_PRETTY, BASES
from nanomotif.model import BetaBernoulliModel
from nanomotif.utils import subseq_indices, calculate_match_length
from nanomotif.seq import EqualLengthDNASet, DNAsequence, DNAarray
from nanomotif.motif import Motif, MotifTree, MotifSearchResult
from nanomotif.logger import configure_logger
from nanomotif.parallel import update_progress_bar
from nanomotif.seed import set_seed
import heapq as hq
import networkx as nx
import warnings
import time
from typing import Generator, Optional, List, Tuple
from dataclasses import dataclass
import pysam
import pyfastx
try:
    import nanomotif_rs
except ImportError:
    nanomotif_rs = None
import shutil
import sys
import subprocess
import threading
from pyinstrument import Profiler


def start_profiler():
    """
    Start profiler if enabled
    """
    if os.getenv("PROFILE_WORKERS") == "1":
        profiler = Profiler(interval = 0.0001)
        profiler.start()
        return profiler
    return None

def stop_profiler(profiler: Profiler, profile_name: str, outdir: str):
    """
    Stop profiler and save output
    """
    if profiler is not None:
        os.makedirs(outdir, exist_ok=True)
        profiler.stop()
        profiler.write_html(f"{outdir}/{profile_name}.html")

def load_bin_assembly(assembly_path: str, contigs: list):
    """The bin's contigs from the indexed FASTA (pyfastx), as an Assembly."""
    fasta = pyfastx.Fasta(assembly_path, build_index=True)
    return nm.seq.Assembly({contig: str(fasta[contig]) for contig in contigs if contig in fasta})


def set_polars_env():
    os.environ["POLARS_MAX_THREADS"] = "1"

@dataclass
class ProcessorConfig:
    assembly: object
    pileup_path: str
    bin_contig: pl.DataFrame
    threads: int
    search_frame_size: int
    methylation_threshold_low: float
    methylation_threshold_high: float
    minimum_kl_divergence: float
    score_threshold: float
    log_dir: str
    seed: int
    output_dir: str
    verbose: bool = False
    min_valid_coverage: int = 5
    # FASTA the assembly was read from; workers then fetch their bin's contigs from its index
    assembly_path: Optional[str] = None
    # file that keeps the filtered pileup partitions for later runs (--pileup_cache); None: temporary
    pileup_cache: Optional[str] = None

    def __post_init__(self):
        if self.assembly is None:
            raise ValueError("assembly cannot be None")
        
        if self.pileup_path is None:
            raise ValueError("pileup_path cannot be None")

        if self.bin_contig is None or self.bin_contig.is_empty():
            raise ValueError("bin_contig cannot be None or empty")
        if "contig" not in self.bin_contig.columns or "bin" not in self.bin_contig.columns:
            raise ValueError("bin_contig must contain 'contig' and 'bin' columns")
        assembly_contigs = (
            set(self.assembly.assembly.keys())
            if hasattr(self.assembly, "assembly")
            else set(self.assembly.keys())
        )
        contigs_in_bin_file = set(self.bin_contig.get_column("contig").to_list())
        missing_contigs = sorted(contigs_in_bin_file - assembly_contigs)
        if missing_contigs:
            preview = ", ".join(missing_contigs[:5])
            suffix = "..." if len(missing_contigs) > 5 else ""
            log.warning(
                "Removing %d contig-bin assignments because the contigs are "
                "absent from the assembly: %s%s",
                len(missing_contigs),
                preview,
                suffix,
            )
            self.bin_contig = self.bin_contig.filter(pl.col("contig").is_in(list(assembly_contigs)))
        if self.bin_contig.is_empty():
            raise ValueError("No contigs remain in bin_contig after filtering against the assembly")
        
        if self.threads <= 0:
            raise ValueError("threads must be greater than 0")

        if self.min_valid_coverage < 0:
            raise ValueError("min_valid_coverage must be >= 0")

        if self.search_frame_size <= 1:
            raise ValueError("search_frame_size must be greater than 1")

        if not (0 <= self.methylation_threshold_high <= 1):
            raise ValueError("methylation_threshold_high must be in [0,1]")

        if not (0 <= self.methylation_threshold_low <= 1):
            raise ValueError("methylation_threshold_low must be in [0,1]")

        if self.methylation_threshold_high <= self.methylation_threshold_low:
            raise ValueError("methylation_threshold_high must be greater than methylation_threshold_low")

        if self.minimum_kl_divergence <= 0:
            raise ValueError("minimum_kl_divergence must be > 0")
        
        if not (0 <= self.score_threshold):
            raise ValueError("score_threshold must be > 0")

    @property
    def padding(self):
        return self.search_frame_size // 2

class TaskStrategy:
    # The pool pickles worker_function, and with it the strategy, for every task. Workers get their
    # bin's data in the task itself, so the partitions and the full assembly are left out.
    _not_sent_to_workers = ("partitioned_assembly", "partitioned_pileup")

    def __getstate__(self):
        state = self.__dict__.copy()
        for name in self._not_sent_to_workers:
            state[name] = None
        state["config"] = copy.copy(self.config)
        state["config"].assembly = None
        # No polars objects either: the pool pickles tasks in one thread while another unpickles results
        # (polars DataFrames). Polars (de)serialisation in two threads at once can deadlock (GIL held while
        # waiting for polars' thread pool): ZymoFecal hung for 30 min at 522 of 1470 tasks. Workers use
        # the bin_contig dict, not this DataFrame.
        state["config"].bin_contig = None
        return state

    def n_tasks(self):
        raise NotImplementedError
    def cleanup(self):
        pass
    def task_generator(self):
        raise NotImplementedError
    def worker_function(self, args):
        raise NotImplementedError

class PileupStrategyPlain(TaskStrategy):
    def __init__(
            self, 
            partitioned_assembly,
            partitioned_pileup: dict,
            bin_contig: dict[str, list[str]],
            config: ProcessorConfig
        ):
        self.partitioned_assembly = partitioned_assembly
        self.partitioned_pileup = partitioned_pileup
        self.bin_contig = bin_contig
        self.config = config

    def n_tasks(self):
        return len(self.partitioned_assembly) * len(nm.constants.MOD_CODE_TO_PRETTY)

    def cleanup(self):
        # a pileup cache (tuple values) is kept; temporary partition files (lists) are removed
        for paths in self.partitioned_pileup.values():
            for path in paths if isinstance(paths, list) else []:
                if os.path.exists(path):
                    os.remove(path)

    def task_generator(self, counter, lock):
        for bin_name in self.partitioned_assembly.keys():
            for modtype in nm.constants.MOD_CODE_TO_PRETTY.keys():
                if (bin_name, modtype) not in self.partitioned_pileup:
                    log.warning(f"Bin {bin_name} has no pileup data, skipping")
                    with lock:
                        counter.value += 1
                    continue
                if bin_name not in self.bin_contig.keys():
                    log.debug(f"Bin {bin_name} not in bin_contig, skipping")
                    with lock:
                        counter.value += 1
                    continue
                bin_pileup = self.partitioned_pileup[(bin_name, modtype)]
                bin_assembly = None if self.config.assembly_path else self.partitioned_assembly[bin_name]
                yield (bin_name, modtype, bin_assembly, bin_pileup)

    def worker_function(self, args):
        bin_name, modtype, assembly, bin_pileup = args
        if isinstance(bin_pileup, list):
            # Arrow IPC files written by the builder, read instead of sent through the pool
            bin_pileup = pl.concat([pl.read_ipc(path) for path in bin_pileup])
        elif isinstance(bin_pileup, tuple):
            # (cache file, [batch numbers]) in the pileup cache
            bin_pileup = nm.pileup_cache.read_partition(*bin_pileup)
        if assembly is None:
            assembly = load_bin_assembly(self.config.assembly_path, self.bin_contig[bin_name])
        profiler = start_profiler()
        set_seed(seed=self.config.seed)
        process_id = os.getpid()
        if self.config.log_dir:
            log_file = f"{self.config.log_dir}/find-motifs.{process_id}.log"
            configure_logger(log_file, verbose=self.config.verbose)
        log.info(f"[Worker {process_id}] started")
        result = process_subpileup(
            {bin_name: self.bin_contig[bin_name]},
            modtype,
            bin_pileup,
            assembly,
            self.config.minimum_kl_divergence,
            self.config.padding,
            self.config.methylation_threshold_low,
            self.config.methylation_threshold_high,
            self.config.score_threshold,
            output_dir = self.config.output_dir
        )
        stop_profiler(profiler, f"{bin_name}_worker_{process_id}_{modtype}", f"{self.config.output_dir}/profiles")
        return result
    
class PileupStrategyBgzip(TaskStrategy):
    def __init__(
            self, 
            partitioned_assembly,
            pileup_path: str, 
            bin_contig: dict[str, list[str]], 
            config: ProcessorConfig
        ):
        self.partitioned_assembly = partitioned_assembly
        self.pileup_path = pileup_path
        self.bin_contig = bin_contig
        self.config = config

    def n_tasks(self):
        return len(self.partitioned_assembly)

    def task_generator(self, counter, lock):
        for bin_name in self.partitioned_assembly.keys():
            if bin_name not in self.bin_contig.keys():
                log.info(f"Bin {bin_name} not in bin_contig, skipping")
                with lock:
                    counter.value += 1
                continue

            bin_assembly = None if self.config.assembly_path else self.partitioned_assembly[bin_name]
            yield (
                bin_name,
                bin_assembly
            )
    def worker_function(self, args):
        bin_name, assembly = args
        if assembly is None:
            assembly = load_bin_assembly(self.config.assembly_path, self.bin_contig[bin_name])
        profiler = start_profiler()
        set_seed(seed=self.config.seed)
        process_id = os.getpid()
        if self.config.log_dir is not None:
            log_file = f"{self.config.log_dir}/find-motifs.{process_id}.log"
            configure_logger(log_file, verbose=self.config.verbose)
        log.info(f"Worker {process_id} started")
        
        log.debug(f"Worker {process_id} processing bin {bin_name}")
        results = []

        log.debug(f"Loading pileup for {bin_name}")
        bin_pileup = nm.dataload.load_contigs_pileup_bgzip(self.pileup_path, self.bin_contig[bin_name])
        log.debug(f"Filtering pileup for {bin_name}")
        bin_pileup = nm.dataload.filter_pileup(bin_pileup, min_coverage=self.config.min_valid_coverage)
        if bin_pileup is None or len(bin_pileup) == 0:
            log.info(f"Bin {bin_name} has no pileup data after filtering, skipping")
            return None
        bin_pileup = nm.dataload.filter_pileup_minimummod_frequency(bin_pileup, methylation_threshold=self.config.methylation_threshold_high)
        if bin_pileup is None or len(bin_pileup) == 0:
            log.info(f"Bin {bin_name} has no pileup data after filtering by minimum mod frequency, skipping")
            return None
        bin_pileup = nm.dataload.filter_pileup_adjacency_filter(bin_pileup, methylation_threshold=self.config.methylation_threshold_high)
        if bin_pileup is None or len(bin_pileup) == 0:
            log.info(f"Bin {bin_name} has no pileup data after adjacency filtering, skipping")
            return None
        log.debug(f"Loaded pileup for {bin_name}")

        for modtype in MOD_CODE_TO_PRETTY.keys():
            log.debug(f"Procesing modtype {modtype}")
            bin_pileup_mod_type = bin_pileup.filter(pl.col("mod_type") == modtype)
            if bin_pileup_mod_type is None or len(bin_pileup_mod_type) == 0:
                log.info(f"Bin {bin_name} modtype {modtype} has no pileup data, skipping")
                continue

            
            
            contigs_in_pileup = bin_pileup_mod_type.get_column("contig").unique().to_list()
            contigs_in_assembly = set(assembly.assembly.keys())
            contigs_to_use = [
                contig for contig in self.bin_contig[bin_name]
                if contig in contigs_in_pileup and contig in contigs_in_assembly
            ]
            dropped_from_assembly = [
                contig for contig in self.bin_contig[bin_name]
                if contig in contigs_in_pileup and contig not in contigs_in_assembly
            ]
            if dropped_from_assembly:
                log.warning(
                    "Skipping %d contigs from bin %s for modtype %s because they are "
                    "missing in the assembly partition: %s%s",
                    len(dropped_from_assembly),
                    bin_name,
                    modtype,
                    ", ".join(dropped_from_assembly[:5]),
                    "..." if len(dropped_from_assembly) > 5 else "",
                )
            if len(contigs_to_use) == 0:
                log.warning(
                    "Bin %s modtype %s has no contigs present in both pileup and assembly, skipping",
                    bin_name,
                    modtype,
                )
                continue
            filtered_bin_contig = {bin_name: contigs_to_use}
            log.debug(f"Contigs in pileup for {bin_name} {modtype}: {contigs_to_use}")
            assembly_modtype = nm.seq.Assembly(
                {contig: assembly[contig].sequence for contig in contigs_to_use}
            )

            results.append(process_subpileup(
                filtered_bin_contig,
                modtype,
                bin_pileup_mod_type,
                assembly_modtype,
                self.config.minimum_kl_divergence,
                self.config.padding,
                self.config.methylation_threshold_low,
                self.config.methylation_threshold_high,
                self.config.score_threshold,
                output_dir = self.config.output_dir
            ))
        results = [result for result in results if result is not None]
        if len(results) == 0:
            log.info(f"No motifs found for bin {bin_name}")
            return None
        
        motifs = pl.concat(results, rechunk=True, parallel=False)
        motifs = nm.motif.MotifSearchResult(motifs)
        
        stop_profiler(profiler, f"{bin_name}_worker_{process_id}", f"{self.config.output_dir}/profiles")
        return motifs


class BinMotifProcessor:
    def __init__(
        self,
        task_strategy: TaskStrategy,
        config: ProcessorConfig,
    ):
        self.task_strategy = task_strategy
        self.config = config
        self.ctx = get_context("spawn")
        self.manager = self.ctx.Manager()
        self.counter = self.manager.Value("i", 0)
        self.lock = self.manager.Lock()



    def run(self):
        # Create a pool of workers
        pool = self.ctx.Pool(processes=self.config.threads, initializer=set_polars_env)

        # Create a process for the progress bar
        log.debug("Counting number of tasks")
        total_tasks = self.task_strategy.n_tasks()

        log.debug("Starting progress bar")
        progress_bar_process = self.ctx.Process(target=update_progress_bar, args=(self.counter, total_tasks, True))
        progress_bar_process.start()

        # Put them workers to work
        log.debug("Starting workers")
        results = []
        for r in pool.imap_unordered(
            self.task_strategy.worker_function, 
            self.task_strategy.task_generator(self.counter, self.lock), 
            chunksize=1
            ):
            with self.lock:
                self.counter.value += 1
            results.append(r)
        log.debug("Joining results")
        results = [nm.motif.MotifSearchResult(result) for result in results if result is not None]
        log.debug("Closing pool")
        pool.close()
        pool.join()
        self.task_strategy.cleanup()

        log.debug("Joining progress bar")
        progress_bar_process.join()

        log.debug("Creating final dataframe")
        if len(results) == 0:
            log.debug("No motifs found")
            return None

        motifs = pl.concat(results, rechunk=True, parallel=False)
        motifs = nm.motif.MotifSearchResult(motifs)
        log.debug("Returning final dataframe")
        return motifs

class BinMotifProcessorBuilder:
    def __init__(self, config: ProcessorConfig):
        self.config = config
        self.bin_contig_dict = self._contig_bin_dict()
    
    def build(self) -> BinMotifProcessor:
        pileup_path = self.config.pileup_path

        if pileup_path.endswith(".gz") and nanomotif_rs is None:
            if self.config.pileup_cache:
                log.warning("--pileup_cache needs nanomotif-rs for bgzipped pileups; not using it")
            # without the compiled reader, each worker reads its bin through the tabix index
            if not os.path.exists(pileup_path + ".tbi"):
                raise FileNotFoundError(f"Tabix index for {pileup_path} not found.")
            partitioned_assembly = self._partition_assembly()
            # Remove contigs not in tabix index from contig_bin_dict
            with pysam.TabixFile(pileup_path) as tbx:
                indexed_contigs = set(tbx.contigs)
            for bin_name, contigs in self.bin_contig_dict.items():
                self.bin_contig_dict[bin_name] = [c for c in contigs if c in indexed_contigs]
            strategy = PileupStrategyBgzip(
                partitioned_assembly=partitioned_assembly,
                pileup_path=pileup_path,
                bin_contig=self.bin_contig_dict,
                config=self.config,
            )
        else:
            # One lz4 Arrow IPC file (or several) per bin and mod type, read by the workers; with
            # --pileup_cache, all of them in one file that later runs reuse
            contig_bin = dict(self.config.bin_contig.select("contig", "bin").iter_rows())
            cache = self.config.pileup_cache
            partitioned_pileup = None
            if cache:
                manifest = nm.pileup_cache.make_manifest(
                    pileup_path, contig_bin, self.config.min_valid_coverage, self.config.methylation_threshold_high
                )
                partitioned_pileup = nm.pileup_cache.load(cache, manifest)
            if partitioned_pileup is None:
                partition_dir = os.path.join(self.config.output_dir, "temp", "pileup_partitions")
                os.makedirs(partition_dir, exist_ok=True)
                if nanomotif_rs is not None and cache:
                    log.info(f"Partitioning the pileup into {cache} (streaming, nanomotif_rs)")
                    partitioned_pileup = nm.pileup_cache.write_with_rust(
                        nanomotif_rs, cache, manifest, pileup_path, contig_bin,
                        self.config.min_valid_coverage, self.config.methylation_threshold_high, self.config.threads,
                    )
                elif nanomotif_rs is not None:
                    log.info("Partitioning the pileup (streaming, nanomotif_rs)")
                    partitioned_pileup = nanomotif_rs.partition_pileup(
                        pileup_path,
                        contig_bin,
                        partition_dir,
                        min_coverage=self.config.min_valid_coverage,
                        threshold=self.config.methylation_threshold_high,
                        threads=self.config.threads,
                    )
                else:
                    partitioned_pileup = self._partition_plain_python(pileup_path, partition_dir, cache, manifest if cache else None)
            if not partitioned_pileup:
                log.info(f"No pileup data after filtering, skipping")
                return None
            partitioned_assembly = self._partition_assembly()
            strategy = PileupStrategyPlain(
                partitioned_assembly=partitioned_assembly,
                partitioned_pileup=partitioned_pileup,
                bin_contig=self.bin_contig_dict,
                config=self.config,
            )

        return BinMotifProcessor(task_strategy=strategy, config=self.config)

    def _partition_plain_python(self, pileup_path, partition_dir, cache=None, manifest=None):
        """Pure-Python version of nanomotif_rs.partition_pileup: load the whole pileup, filter, split."""
        log.warning("Pileup is not bgzipped and nanomotif_rs is not installed, loading full pileup into memory.")
        pileup = nm.dataload.load_pileup(
            pileup_path,
            min_coverage=self.config.min_valid_coverage,
            contigs=self.config.bin_contig.get_column("contig").to_list()
        )
        pileup = nm.dataload.filter_pileup_minimummod_frequency(pileup, methylation_threshold=self.config.methylation_threshold_high)
        pileup = nm.dataload.filter_pileup_adjacency_filter(pileup, methylation_threshold=self.config.methylation_threshold_high)
        pileup = pileup.join(self.config.bin_contig, left_on="contig", right_on="contig", how="inner")
        if cache:
            parts = pileup.partition_by(["bin", "mod_type"], as_dict=True)
            return nm.pileup_cache.write_from_frames(cache, manifest, parts)
        partitioned_pileup = {}
        for i, ((bin_name, mod_type), part) in enumerate(pileup.partition_by(["bin", "mod_type"], as_dict=True).items()):
            path = os.path.join(partition_dir, f"{i}.arrow")
            part.write_ipc(path, compression="lz4")
            partitioned_pileup[(bin_name, mod_type)] = [path]
        return partitioned_pileup

    def _partition_assembly(self):
        if self.config.assembly_path:
            # workers read their contigs from the FASTA index
            return {bin_name: None for bin_name in self.bin_contig_dict}
        return {
            bin_name: nm.seq.Assembly(
                {contig: self.config.assembly[contig].sequence for contig in contigs if contig in self.config.assembly}
            )
            for bin_name, contigs in self.bin_contig_dict.items()
        }
    def _contig_bin_dict(self):
        bin_contig_dict = {}
        for row in self.config.bin_contig.iter_rows(named=True):
            bin_name = row["bin"]
            contig_name = row["contig"]
            if bin_name not in bin_contig_dict:
                bin_contig_dict[bin_name] = []
            bin_contig_dict[bin_name].append(contig_name)
        return bin_contig_dict























def process_subpileup(
        bin_contig: dict, modtype, bin_pileup, assembly, 
        min_kl_divergence, padding, 
        low_meth_threshold, high_meth_threshold,
        score_threshold,
        output_dir = None
    ):
    """
    Process a single subpileup for one bin and one modtype
    """
    bin_name = list(bin_contig.keys())[0]
    log.info(f"Processing {bin_name} {modtype}")
    log.debug(f"Contigs in bin: {list(bin_contig[bin_name])}")
    assert bin_pileup is not None, "Subpileup is None"
    assert len(bin_pileup) > 0, "Subpileup is empty"
    assert assembly is not None, "Assembly is None"
    assert min_kl_divergence >= 0, "min_kl_divergence must be greater than 0"
    assert bin_pileup.get_column("mod_type").unique().to_list() == [modtype], "subpileup modtype does not match modtype"

    log.debug(f"Subpileup for {bin_name} {modtype}: {bin_pileup}")
    assembly_store = assembly.assembly if hasattr(assembly, "assembly") else assembly
    available_contigs = set(assembly_store.keys())
    requested_contigs = bin_contig.get(bin_name, [])
    missing_contigs = [contig for contig in requested_contigs if contig not in available_contigs]
    if missing_contigs:
        log.warning(
            "Skipping %d contigs in bin %s because they are missing from the assembly: %s%s",
            len(missing_contigs),
            bin_name,
            ", ".join(missing_contigs[:5]),
            "..." if len(missing_contigs) > 5 else "",
        )
    contigs_to_use = [contig for contig in requested_contigs if contig in available_contigs]
    if len(contigs_to_use) == 0:
        log.warning(f"Bin {bin_name} {modtype} has no contigs present in the assembly, skipping")
        return None
    bin_sequences = {}
    for contig in contigs_to_use:
        try:
            bin_sequences[contig] = assembly[contig]
        except KeyError:
            log.warning(
                "Encountered missing contig %s in bin %s despite earlier filtering, skipping",
                contig,
                bin_name,
            )
    if len(bin_sequences) == 0:
        log.warning(f"Bin {bin_name} {modtype} has no contigs with sequences after filtering, skipping")
        return None

    # Find motifs
    site_index = MotifSiteIndex(bin_pileup, bin_sequences, low_meth_threshold, high_meth_threshold, half_width=padding)
    candidates = find_best_candidates(
        bin_pileup, 
        bin_sequences, 
        modtype,
        bin_name,
        output_dir=output_dir,
        site_index=site_index,
        low_meth_threshold=low_meth_threshold,
        high_meth_threshold=high_meth_threshold,
        padding=padding,
        min_kl=min_kl_divergence,
        max_dead_ends=25,
        max_rounds_since_new_best=30,
        score_threshold=score_threshold
    )
    if candidates is None:
        log.info("No motifs found")
        return None
    motif_graph, best_candidates = candidates
    identified_motifs = nxgraph_to_dataframe(motif_graph) \
        .filter(col("sequence").is_in(best_candidates))
    
    if len(identified_motifs) == 0:
        log.info("No motifs found")
        return None
    
    identified_motifs = identified_motifs.with_columns(
        pl.lit(padding).alias("mod_position"),
        pl.lit(bin_name).alias("reference"),
        pl.lit(modtype).alias("mod_type"),
    ).rename({"sequence": "motif"})
    motifs = nm.motif.MotifSearchResult(identified_motifs)

    # Post-process motifs
    log.info("Refining motifs")
    if output_dir is None:
        output_dir = "./"
    motifs_outdir = output_dir + "/precleanup-motifs/" + bin_name + "-" + modtype
    os.makedirs(motifs_outdir, exist_ok=True)
    motifs_file_name = motifs_outdir + "/motifs"
    motifs.write_motifs(motifs_file_name + ".tsv")

    log.info(" - Removing noisy motifs")
    motifs = nm.postprocess.remove_noisy_motifs(motifs)
    if motifs.is_empty():
        log.info("No motifs found")
        return None
    motifs_file_name = motifs_file_name +   "-noise"
    motifs.write_motifs(motifs_file_name + ".tsv")

    log.info(" - Merging motifs")
    motifs = nm.find_motifs_bin.merge_motifs_in_df(
        motifs, bin_pileup, assembly, bin_contig,
        low_meth_threshold=low_meth_threshold, high_meth_threshold=high_meth_threshold,
        site_index=site_index
    )
    motifs = motifs.unique()
    if motifs.is_empty():
        log.info("No motifs found")
        return None
    motifs_file_name = motifs_file_name +   "-merge"
    motifs.write_motifs(motifs_file_name + ".tsv")

    log.info(" - Removing sub motifs")
    motifs = nm.postprocess.remove_sub_motifs(motifs)
    motifs = motifs.unique()
    if motifs.is_empty():
        log.info("No motifs found")
        return None
    motifs_file_name = motifs_file_name +   "-sub"
    motifs.write_motifs(motifs_file_name + ".tsv")

    log.info(" - Joining motif complements")
    motifs = nm.postprocess.join_motif_complements(motifs)
    motifs = motifs.unique()
    if motifs.is_empty():
        log.info("No motifs found")
        return None
    motifs_file_name = motifs_file_name +   "-complement"
    motifs.write_motifs(motifs_file_name + ".tsv")


    return motifs






#########################################################################
# Motif candidate state space 

def find_best_candidates(        
        bin_pileup: pl.DataFrame, 
        bin_sequences: dict[str, DNAsequence],
        mod_type: str,
        bin_name: str,
        output_dir: str,
        low_meth_threshold: float,
        high_meth_threshold: float,
        padding: int,
        min_kl: float = 0.2, 
        max_dead_ends: int = 25, 
        max_rounds_since_new_best: int = 30,
        score_threshold: float = 0.2,
        remaining_sequences_threshold: float = 0.001,
        background_sampling_frequency: float = 0.01,
        site_index: Optional["MotifSiteIndex"] = None
    ) -> tuple[MotifTree, list[Motif]]:
    """
    Find the best motif candidates in a sequence.
    """
    if site_index is None:
        site_index = MotifSiteIndex(bin_pileup, bin_sequences, low_meth_threshold, high_meth_threshold, half_width=padding)
    subpileup_confident = bin_pileup.filter(pl.col("fraction_mod") >= high_meth_threshold) \
        .partition_by("contig", as_dict=True)
    no_sites = bin_pileup.clear()
    # collected as lists and turned into sets once: adding EqualLengthDNASets re-validates everything so far
    methylation_sequences = []
    background_sequences = []
    os.makedirs(f"{output_dir}/temp/{bin_name}", exist_ok=True)
    # sorted, so contigs draw their random background windows in the same order every run
    for contig in bin_pileup.get_column("contig").unique().sort():
        subpileup_confident_contig = subpileup_confident.get((contig,), no_sites)
        contig_sequence = bin_sequences[contig]
        contig_length = len(contig_sequence)
        n_samples = int(max(math.ceil(contig_length * background_sampling_frequency), 50))
        # Extract the sequences for confidently methylated positions
        index_plus = subpileup_confident_contig.filter(pl.col("strand") == "+").get_column("position").to_list()
        index_minus = subpileup_confident_contig.filter(pl.col("strand") == "-").get_column("position").to_list()
        
        # Sample sequences in bin to get background for KL-divergence
        log.debug(f"Sampling {n_samples} sequences from {contig}")
        contig_sequences_sample = contig_sequence.sample_n_subsequences_unique(
            padding * 2 + 1, n_samples, base=MOD_TYPE_TO_CANONICAL[mod_type]
        )

        background_sequences += contig_sequences_sample.sequences

        log.debug(f"Extracting {len(index_plus)} + strand and {len(index_minus)} - strand methylation sequences from {contig}")
        methylation_sequences_string = []
        if len(index_plus) >= 1:
            methylation_sequences_temp = contig_sequence.sample_at_indices(index_plus, padding)
            if methylation_sequences_temp is not None:
                methylation_sequences_string += methylation_sequences_temp.sequences
        if len(index_minus) >= 1:
            methylation_sequences_temp = contig_sequence.sample_at_indices(index_minus, padding)
            if methylation_sequences_temp is not None:
                methylation_sequences_string += methylation_sequences_temp.reverse_compliment().sequences
        if len(methylation_sequences_string) == 0:
            log.debug(f"No methylation sequences found in {contig}")
            continue
        methylation_sequences += methylation_sequences_string
    # Convert methylation sequences to DNAarray
    if len(methylation_sequences) == 0:
        log.info("No methylation sequences found")
        return None
    methylation_sequences = EqualLengthDNASet(methylation_sequences)
    background_sequences = EqualLengthDNASet(background_sequences) if background_sequences else None
    with open(f"{output_dir}/temp/{bin_name}/methylation_sequences_{mod_type}.fasta", "w") as f:
        for i, seq in enumerate(methylation_sequences.sequences):
            f.write(f">methylation_seq_{i}\n{seq.sequence}\n")
    if background_sequences is not None:
        with open(f"{output_dir}/temp/{bin_name}/background_sequences_{mod_type}.fasta", "w") as f:
            for i, seq in enumerate(background_sequences.sequences):
                f.write(f">background_seq_{i}\n{seq.sequence}\n")
    methylation_sequences = nm.seq.PackedDNAarray.from_sequences(methylation_sequences)
    total_seqs = methylation_sequences.shape[0]

    # Convert background sequences to DNAarray
    if background_sequences is None or len(background_sequences.sequences) == 0:
        log.info("No background sequences found")
        return None

    bin_pssm = background_sequences.pssm()
    np.savetxt(f"{output_dir}/temp/{bin_name}/background_pssm_{mod_type}.txt", bin_pssm, fmt="%.4f")

    root_motif = Motif("." * padding + MOD_TYPE_TO_CANONICAL[mod_type] + "." * padding, padding)
    methylation_sequences_clone = methylation_sequences
    best_candidates = []
    continue_search = True
    dead_ends = 0
    motif_graph = None

    while continue_search:
        log.debug(f"Searching for best candidates, {len(best_candidates)} found so far")
        if dead_ends >= max_dead_ends:
            log.debug("Stopping search, too many low scoring candidates")
            break
        # Find the initial guess within the tree
        searcher = MotifSearcher(
            root_motif, 
            bin_sequences, 
            bin_pssm,
            bin_pileup, 
            methylation_sequences_clone, 
            padding,
            high_meth_threshold,
            low_meth_threshold,
            motif_graph = motif_graph,
            min_kl = min_kl,
            max_rounds_since_new_best = max_rounds_since_new_best,
            site_index = site_index
        )
        motif_graph, naive_guess = searcher.run()

        # If there is no naive guess, stop the search
        if naive_guess == root_motif:
            log.debug("No naive guess found, stopping search")
            break

        # Prune the motif for low scoring positions
        pruning = True
        pruning_round = 1
        temp_motif = naive_guess
        indices_to_prune = set()
        pruned_to_single_base = False
        while pruning:
            log.debug(f"Pruning round {pruning_round} for motif {temp_motif}")
            parent_motif_scores = get_parent_scores(
                temp_motif,
                bin_pileup,
                bin_sequences,
                low_meth_threshold,
                high_meth_threshold,
                site_index=site_index
            )
            parent_scores = [d["score"] for d in parent_motif_scores.values()]
            mean_parent_score = np.mean(parent_scores)
            log.debug(f"Mean parent score: {mean_parent_score}")
            parent_pretty_print = ", ".join([f"{parent} (pos {data['motif_position']}, score {data['score']:.2f}, model {data['parent_model']})" for parent, data in parent_motif_scores.items()])
            log.debug(f"Parents: {parent_pretty_print}")
            # Find positions to prune
            for parent, data in parent_motif_scores.items():
                if data["score"] < 0.4:
                    log.debug(f"    Pruning position {data['motif_position']} from {temp_motif}, score: {data['score']}")
                    log.debug(f"    Poor parent motif: {parent}, model: {data['parent_model']}")
                    indices_to_prune.add(data["motif_position"])
            log.debug(f"Pruning round {pruning_round} complete, pruned {len(indices_to_prune)} positions.")

            if len(indices_to_prune) == 0:
                log.debug("No more positions to prune")
                pruning = False
                break

            motif_split = temp_motif.split()
            for i in indices_to_prune:
                motif_split[i] = "."
            pruned_motif_str = "".join(motif_split)
            pruned_motif = Motif(pruned_motif_str, temp_motif.mod_position)
            log.debug(f"Pruned motif to {pruned_motif}")
            if len(pruned_motif_str.replace(".", "")) == 1:
                pruned_to_single_base = True
                pruning = False
                break
            if pruned_motif == temp_motif:
                pruning = False
                break
            pruning_round += 1
            temp_motif = pruned_motif
        
        if pruned_to_single_base:
            log.debug("Motif reduced to single base, removing original motif from methylation sequences")
            log.debug("Updating score based on parents")
            mean_parent_score = np.mean([d["score"] for d in parent_motif_scores.values()])
            motif_graph.nodes[naive_guess]["score"] = mean_parent_score
        elif mean_parent_score < score_threshold:
            log.debug(f"Pruned motif {temp_motif} has low score {mean_parent_score}, discarding")
            mean_parent_score = np.mean([d["score"] for d in parent_motif_scores.values()])
            motif_graph.nodes[naive_guess]["score"] = mean_parent_score

        elif temp_motif != naive_guess:
            log.debug(f"Pruned motif from {naive_guess} to {temp_motif}")
            
            child_model = [d["child_model"] for _, d in parent_motif_scores.items()][0]
            # Add the pruned motif to the graph
            motif_graph.add_node(
                temp_motif,
                model=child_model,
                motif=temp_motif,
                visited=True,
                score=mean_parent_score,
                priority=0,
                depth=0
            )
            naive_guess = temp_motif
        else:
            log.debug(f"No pruning performed for {temp_motif}")
            log.debug("Updating score based on parents")
            mean_parent_score = np.mean([d["score"] for d in parent_motif_scores.values()])
            motif_graph.nodes[naive_guess]["score"] = mean_parent_score

        # Remove new candidate from methylation sequences
        seq_before = methylation_sequences_clone.shape[0]
        methylation_sequences_clone = methylation_sequences_clone.filter_sequence_matches(naive_guess.one_hot(), keep_matches = False)
        if methylation_sequences_clone is None:
            log.debug("No more sequences left")
            break

        # Check if we should continue the search
        seq_remaining = methylation_sequences_clone.shape[0]
        seq_remaining_percent = seq_remaining/total_seqs

        if motif_graph.nodes[naive_guess]["score"] < score_threshold:
            dead_ends += 1
            log.debug(f"Candidate has low score ({motif_graph.nodes[naive_guess]['score']}), {naive_guess}. {dead_ends} of {max_dead_ends} before termination. Represented in {seq_before-seq_remaining} seqs. model: {motif_graph.nodes[naive_guess]['model']}. ({100*seq_remaining_percent:.1f} % of sequences remaining)")
            continue

        log.info(f"Keeping {naive_guess}, represented in {seq_before-seq_remaining} seqs. model: {motif_graph.nodes[naive_guess]['model']}. ({100*seq_remaining_percent:.1f} % of sequences remaining)")

        best_candidates.append(naive_guess)
        
        if (seq_remaining/total_seqs) < remaining_sequences_threshold:
            log.debug("Stopping search, too few sequences remaining")
            break
        log.debug("Continuing search")
    if motif_graph is None or len(motif_graph.nodes) == 0:
        log.info("No motifs found")
        return None
    # Check if any high scoring candidates were found, best were not high scoring
    missed_candidates = motif_graph.get_missed_candidates(best_candidates, score_threshold)
    missed_candidates = [c for c in missed_candidates if not c.sub_motif_of_any(best_candidates) and not any(candidate.sub_motif_of(c) for candidate in best_candidates)]
    missed_candidates = sorted(missed_candidates, key=lambda c: (-motif_graph.nodes[c]["score"], c.string))
    log.debug(f"Found {len(missed_candidates)} missed candidates")
    for candidate in missed_candidates:
        best_candidates.append(candidate)
        log.info(f"Adding missed candidate {candidate} with score {motif_graph.nodes[candidate]['score']}")
    # Save graph to file
    motif_graph_file = f"{output_dir}/temp/{bin_name}/motif_graph_{mod_type}.gml"
    motif_graph.export_graph_gml(motif_graph_file)

    return motif_graph, best_candidates



class MotifSearcher:
    """
    Class to perform motif search for identifying enriched motifs in methylated sequences.
    """

    def __init__(
        self,
        root_motif: Motif,
        bin_sequence: DNAsequence,
        bin_pssm: DNAarray,
        bin_pileup: pl.DataFrame,
        methylation_sequences: DNAarray,
        padding: int,
        high_meth_threshold: float,
        low_meth_threshold: float,
        motif_graph: Optional[MotifTree] = None,
        min_kl: float = 0.1,
        freq_threshold: float = 0.15,
        max_rounds_since_new_best: int = 30,
        max_motif_length: int = 25,
        site_index: Optional["MotifSiteIndex"] = None
    ):
        """
        Initialize the MotifSearcher class.

        Parameters:
            root_motif (Motif): The initial motif to start the search from.
            bin_sequence: The bin sequence object with necessary methods.
            bin_pssm: The position-specific scoring matrix for the bin background.
            bin_pileup: The pileup data associated with the bin.
            methylation_sequences: The methylation sequences to be analyzed.
            padding (int): The padding size for sequence sampling.
            read_level_methylation (bool): Flag to use read-level methylation.
            motif_graph (Optional[MotifTree]): An optional pre-initialized motif graph.
            min_kl (float): Minimum Kullback-Leibler divergence threshold.
            freq_threshold (float): Frequency threshold for methylation.
            max_rounds_since_new_best (int): Max iterations without improvement.
            max_motif_length (int): Maximum allowed motif length.
        """
        # Input validation
        if bin_pileup.is_empty():
            raise ValueError("bin_pileup cannot be empty.")
        if padding < 0:
            raise ValueError("padding must be non-negative.")

        self.root_motif = root_motif
        self.bin_sequence = bin_sequence
        self.bin_pssm = bin_pssm
        self.bin_pileup = bin_pileup
        self.methylation_sequences = methylation_sequences
        self.padding = padding
        self.motif_graph = motif_graph or MotifTree()
        self.min_kl = min_kl
        self.freq_threshold = freq_threshold
        self.max_rounds_since_new_best = max_rounds_since_new_best
        self.max_motif_length = max_motif_length
        self.low_meth_threshold = low_meth_threshold
        self.high_meth_threshold = high_meth_threshold
        self.site_index = site_index or MotifSiteIndex(
            bin_pileup, bin_sequence, low_meth_threshold, high_meth_threshold, half_width=padding
        )

    def _priority_function(self, next_model, root_model) -> float:
        """
        Calculate the priority for a motif based on the change in model parameters.

        Parameters:
            next_model: The model of the next motif.
            root_model: The model of the root motif.

        Returns:
            float: The calculated priority.
        """
        # score = self._scoring_function(next_model, root_model)
        # return -score
        # odds_ratio = (next_model._alpha * root_model._beta) / (next_model._beta * root_model._alpha)
        try:
            d_alpha = 1 - (next_model._alpha / root_model._alpha)
        except ZeroDivisionError:
            d_alpha = 1
        try:
            d_beta = next_model._beta / root_model._beta
        except ZeroDivisionError:
            d_beta = 1
        priority = d_alpha * d_beta
        return priority

    def _scoring_function(self, next_model, current_model) -> float:
        """
        Calculate the score for a motif based on its model.

        Parameters:
            next_model: The model of the next motif.
            current_model: The model of the current motif.

        Returns:
            float: The calculated score.
        """
        mean_diff = next_model.mean() - current_model.mean()
        try:
            score = (next_model.mean() * -np.log10(next_model.standard_deviation())) * mean_diff
        except (ValueError, ZeroDivisionError):
            score = 0
        return score

    def _scoring_function_predictive_evaluation(self, next_model, current_model) -> float:
        """
        Calculate the score for a motif based on its predictive model.

        Parameters:
            next_model: The model of the next motif.
            current_model: The model of the current motif.

        Returns:
            float: The calculated score.
        """
        return predictive_evaluation_score(next_model, current_model)

    def _motif_child_nodes_kl_dist_max(
        self,
        motif: Motif,
        meth_pssm: np.ndarray,
        bin_pssm: np.ndarray
    ) -> Generator[Motif, None, None]:
        """
        Generate child motifs based on maximum KL divergence.

        Parameters:
            motif (Motif): The current motif.
            meth_pssm (np.ndarray): Position-specific scoring matrix for methylation.
            bin_pssm (np.ndarray): Position-specific scoring matrix for the bin.

        Yields:
            Motif: The next motif in the search.
        """
        kl_divergence = kl_divergence_columns(meth_pssm, bin_pssm)
        split_motif = motif.split()

        # Positions to evaluate (positions with '.')
        evaluated_positions = np.array([i for i, base in enumerate(split_motif) if base == "."])
        if evaluated_positions.size == 0:
            return

        kl_divergence_masked = kl_divergence.copy()
        kl_divergence_masked[~np.isin(np.arange(len(split_motif)), evaluated_positions)] = 0

        max_kl = np.max(kl_divergence_masked)
        if max_kl < self.min_kl:
            return

        pos = int(np.argmax(kl_divergence_masked))

        # Methylation frequency must be above bin frequency
        index_meth_freq_higher = meth_pssm[:, pos] > bin_pssm[:, pos] * 0.5

        # Methylation frequency must be above a threshold
        index_meth_freq_above_thresh = meth_pssm[:, pos] > self.freq_threshold

        # Combine the two filters
        index_position_filter = np.logical_and(index_meth_freq_higher, index_meth_freq_above_thresh)
        bases_indices = np.argwhere(index_position_filter).reshape(-1)
        bases_filtered = [BASES[int(i)] for i in bases_indices]

        if not bases_filtered:
            return

        # All combinations of the bases
        # max_combination_length = min(len(bases_filtered), 4)
        # for i in range(1, max_combination_length + 1):
        #    for base_tuple in itertools.combinations(bases_filtered, i):
        #        if len(base_tuple) > 2:
        #             continue
        #        if len(base_tuple) > 1:
        #            base_str = "[" + "".join(base_tuple) + "]"
        #        else:
        #            base_str = base_tuple[0]
        #        new_motif_sequence = split_motif.copy()
        #        new_motif_sequence[pos] = base_str
        #        new_motif_str = "".join(new_motif_sequence)
        #        yield Motif(new_motif_str, motif.mod_position)
        for base in bases_filtered:
            new_motif_sequence = split_motif.copy()
            new_motif_sequence[pos] = base
            new_motif_str = "".join(new_motif_sequence)
            yield Motif(new_motif_str, motif.mod_position)

    
    def run(self) -> tuple[MotifTree, Motif]:
        """
        Execute the motif search algorithm.

        Returns:
            Tuple[MotifTree, Motif]: The graph of motifs explored and the best motif found.
        """
        # Initialize variables
        best_guess = self.root_motif
        root_model = self.site_index.model(self.root_motif)

        best_score = self._scoring_function_predictive_evaluation(root_model, root_model)
        rounds_since_new_best = 0
        visited_nodes: set[Motif] = set()

        root_depth = 0
        if not self.motif_graph.has_node(self.root_motif):
            self.motif_graph.add_node(
                self.root_motif,
                model=root_model,
                motif=self.root_motif,
                visited=False,
                score=best_score,
                priority=0,
                depth=root_depth
            )

        # Initialize priority queue 
        priority_queue: list[tuple] = []
        hq.heappush(
            priority_queue,
            (0, root_depth,  self.root_motif)
        )

        # Search loop
        while priority_queue:
            # Get the current best candidate
            _, _, current_motif = hq.heappop(priority_queue)
            if current_motif in visited_nodes:
                continue

            current_attrs = self.motif_graph.nodes[current_motif]
            current_model = current_attrs["model"] 
            current_depth = current_attrs.get("depth", 0)

            current_n_mod, current_n_nomod = current_model.get_raw_counts()
            if current_n_mod + current_n_nomod < 10:
                log.debug(
                    f"{current_motif.string} expansion skipped due to low support ({current_n_mod} mod, {current_n_nomod} non-mod) | Model: {current_model} | "
                    f"Score: {current_attrs['score']:.2f} | "
                    f"Queue size: {len(priority_queue)} | "
                    f"Came from {', '.join(str(node) for node in list(self.motif_graph.predecessors(current_motif)))}"
                )
                continue
            # Enforce motif length safely (current_model is defined)
            if len(current_motif.strip()) > self.max_motif_length:
                log.debug(
                    f"{current_motif.string} expansion skipped due to length | Model: {current_model} | "
                    f"Score: {current_attrs['score']:.2f} | "
                    f"Queue size: {len(priority_queue)} | "
                    f"Came from {', '.join(str(node) for node in list(self.motif_graph.predecessors(current_motif)))}"
                )
                continue

            visited_nodes.add(current_motif)
            log.debug(
                f"{current_motif.string} expanded | Model: {current_model} | "
                f"Score: {current_attrs['score']:.2f} | "
                f"Queue size: {len(priority_queue)} | "
                f"Came from {', '.join(str(node) for node in list(self.motif_graph.predecessors(current_motif)))}"
            )
            self.motif_graph.nodes[current_motif]["visited"] = True
            rounds_since_new_best += 1

            # Filter sequences for this motif
            active_methylation_sequences, active_pssm = self.methylation_sequences.filter_and_pssm(
                current_motif.one_hot(), keep_matches=True
            )
            if active_methylation_sequences is None:
                log.debug("No more sequences left after filtering")
                continue

            # Generate and evaluate neighbor motifs
            neighbors = list(self._motif_child_nodes_kl_dist_max(
                current_motif,
                active_pssm,
                self.bin_pssm
            ))

            # Add neighbors to graph
            for next_motif in neighbors:
                next_depth = current_depth + 1

                # Build model for the next motif
                if next_motif not in self.motif_graph.nodes:
                    next_model = self.site_index.model(next_motif)
                else:
                    next_model = self.motif_graph.nodes[next_motif]["model"]

                # Score relative to current
                score = self._scoring_function_predictive_evaluation(next_model, current_model)

                n_isolated = next_motif.count_isolated_bases(isolation_size=1)
                priority = self._priority_function(next_model, root_model)
                if n_isolated > 0:
                    priority *= pow(10, n_isolated)  # Deprioritize motifs with isolated bases

                if next_motif in self.motif_graph.nodes:
                    # Update existing node if new score is better
                    if self.motif_graph.nodes[next_motif]["score"] < score:
                        self.motif_graph.nodes[next_motif]["score"] = score
                else:
                    self.motif_graph.add_node(
                        next_motif,
                        model=next_model,
                        motif=next_motif,
                        visited=False,
                        score=score,
                        priority=priority,
                        depth=next_depth
                    )
                if not self.motif_graph.has_edge(current_motif, next_motif):
                    self.motif_graph.add_edge(current_motif, next_motif)

                # Add neighbor to priority queue if not visited
                if next_motif not in visited_nodes:
                    attrs = self.motif_graph.nodes[next_motif]
                    hq.heappush(
                        priority_queue,
                        (attrs["priority"], attrs["depth"], next_motif)
                    )

                # Track best
                if score > best_score:
                    best_score = score
                    best_guess = next_motif
                    rounds_since_new_best = 0

            # Stopping criteria
            if rounds_since_new_best >= self.max_rounds_since_new_best:
                break

        return self.motif_graph, best_guess





def kl_divergence_columns(pk: np.ndarray, qk: np.ndarray) -> np.ndarray:
    """scipy.stats.entropy(pk, qk) along axis 0 (columns normalised to 1), without its per-call overhead."""
    pk = pk / pk.sum(axis=0, keepdims=True)
    qk = qk / qk.sum(axis=0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(pk > 0, np.where(qk > 0, pk * np.log(pk / qk), np.inf), np.where(qk >= 0, 0.0, np.inf))
    return terms.sum(axis=0)


def count_periods_at_start(string: str) -> int:
    """
    Count the number of periods at the start of a string

    Parameters:
    - string (str): The string to be processed.

    Returns:
    - int: The number of periods at the start of the string.

    Example:
    >>> count_periods_at_start("...ACGT")
    3
    """
    count = 0
    for char in string:
        if char == '.':
            count += 1
        else:
            break
    return count
def count_periods_at_end(string):
    string = string[::-1]
    count = 0
    for char in string:
        if char == '.':
            count += 1
        else:
            break
    return count
def nxgraph_to_dataframe(graph):
    return pl.DataFrame({
        "sequence":[i for i in graph.nodes],
        "model":[graph.nodes[i]["model"] for i in graph.nodes],
        "score":[graph.nodes[i]["score"] for i in graph.nodes]
    }).sort("score", descending=True)






##########################################
# Motif candidate scoring
##########################################

def methylated_motif_occourances(
        motif, 
        sequence, 
        methylated_positions,
        non_methylated_positions
    ) -> tuple:
    """
    Get occourances of a motif in a bin

    Parameters:
    - motif (str): The motif to search for.
    - seq (str): The bin to search in.
    - bin_meth_positions (list): The positions of the methylation sites in the bin.

    Returns:
    - tuple: A tuple of two numpy arrays, the first containing the positions of the methylated motifs, the second containing the positions of the non-methylated motifs.
    """
    assert len(motif) > 0, "Motif is empty"
    assert len(sequence) > 0, "Sequence is empty"
    assert type(motif) == Motif, "Motif is not a Motif type"
    # Get the index of the methylation in the motif in the bin
    motif_index = subseq_indices(motif.string, sequence) + motif.mod_position

    # Methylated motif positions
    meth_occurences = methylated_positions[np.isin(methylated_positions, motif_index, assume_unique=True)]

    # Non-methylated motif positions
    nonmeth_occurences =  non_methylated_positions[np.isin(non_methylated_positions, motif_index, assume_unique=True)]

    return meth_occurences, nonmeth_occurences

def _pack_bits(mask: np.ndarray) -> np.ndarray:
    """Boolean array -> bitset as uint64 words (padding bits are 0)."""
    packed = np.packbits(mask)
    return np.concatenate([packed, np.zeros(-len(packed) % 8, np.uint8)]).view(np.uint64)


_POPCOUNT_8 = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)


def _popcount(words: np.ndarray) -> int:
    if hasattr(np, "bitwise_count"):  # numpy >= 2.0
        return int(np.bitwise_count(words).sum())
    return int(_POPCOUNT_8[words.view(np.uint8)].sum())


class MotifSiteIndex:
    """
    Sequence windows around the confidently methylated and unmethylated positions of a bin,
    for counting a motif's methylated and unmethylated occurrences without scanning the contigs.

    motif_model_contig finds all occurrences of a motif with a regex over the whole contig and then
    intersects them with the methylated and unmethylated positions. Here the motif is instead tested at
    those positions only, against the +-half_width bases around each (one array gather per fixed motif
    position). The counts are the same; a motif that does not fit in the window is counted with
    motif_model_contig.
    """
    def __init__(self, pileup: pl.DataFrame, contigs: dict, low_meth_threshold: float, high_meth_threshold: float,
                 half_width: int = 20):
        self.pileup = pileup
        self.contigs = contigs
        self.low_meth_threshold = low_meth_threshold
        self.high_meth_threshold = high_meth_threshold
        self.half_width = half_width
        windows = {(strand, meth): [] for strand in "+-" for meth in (True, False)}
        sites = pileup.filter(
            (pl.col("fraction_mod") >= high_meth_threshold) | (pl.col("fraction_mod") <= low_meth_threshold)
        ).select("contig", "position", "strand", (pl.col("fraction_mod") >= high_meth_threshold).alias("methylated"))
        for (contig,), contig_sites in sites.partition_by("contig", as_dict=True).items():
            if contig not in contigs:
                continue
            sequence = np.frombuffer(contigs[contig].sequence.encode("ascii"), dtype=np.uint8)
            # 0 marks positions outside the contig; it matches no motif base, not even "."
            padded = np.concatenate([np.zeros(half_width, np.uint8), sequence, np.zeros(half_width, np.uint8)])
            # a strided view: indexing it copies 1 byte per window base (an index matrix would take 8)
            all_windows = np.lib.stride_tricks.sliding_window_view(padded, 2 * half_width + 1)
            positions = contig_sites.get_column("position").to_numpy()
            on_contig = (positions >= 0) & (positions < len(sequence))
            for strand in "+-":
                for meth in (True, False):
                    keep = on_contig \
                        & (contig_sites.get_column("strand").to_numpy() == strand) \
                        & (contig_sites.get_column("methylated").to_numpy() == meth)
                    windows[(strand, meth)].append(all_windows[positions[keep]])
        width = 2 * half_width + 1
        # One bitset over the sites per (window column, base): a motif's count is then the popcount of
        # the AND over its fixed positions of the OR over each position's allowed bases.
        self.n_sites = {}
        self.bits = {}
        for key in list(windows):
            w = windows.pop(key)  # freed as soon as its bitsets are built
            w = np.concatenate(w) if w else np.empty((0, width), np.uint8)
            self.n_sites[key] = len(w)
            self.bits[key] = {}
            if nanomotif_rs is not None and len(w):
                # one pass over the windows in Rust, same bit layout as _pack_bits
                codes, bits = nanomotif_rs.column_bitsets(np.ascontiguousarray(w))
                for column in range(width):
                    for i, code in enumerate(codes):
                        if bits[column, i].any():
                            self.bits[key][(column, int(code))] = bits[column, i]
                continue
            w = np.ascontiguousarray(w.T)  # one transpose instead of a strided copy per column
            for column in range(width):
                values = w[column]
                for code in np.flatnonzero(np.bincount(values, minlength=256)[1:]) + 1:
                    self.bits[key][(column, int(code))] = _pack_bits(values == code)

    def _count(self, key: tuple, motif: Motif) -> int:
        """Number of indexed sites at which the (stripped) motif occurs with its modified base at the site."""
        match = None
        for i, token in enumerate(motif.split()):
            if token == ".":
                continue
            column = self.half_width - motif.mod_position + i
            allowed = None
            for base in token.strip("[]"):
                bits = self.bits[key].get((column, ord(base)))
                if bits is not None:
                    allowed = bits if allowed is None else allowed | bits
            if allowed is None:
                return 0
            match = allowed if match is None else match & allowed
        return self.n_sites[key] if match is None else _popcount(match)

    def counts(self, motif: Motif) -> tuple:
        """(methylated, unmethylated) occurrences of the motif at the indexed positions."""
        stripped = motif.new_stripped_motif()
        reverse = stripped.reverse_compliment()
        length = stripped.length()
        if max(stripped.mod_position, length - 1 - stripped.mod_position) > self.half_width:
            model = BetaBernoulliModel(0, 0)
            for contig, contig_sequence in self.contigs.items():
                model = motif_model_contig(
                    self.pileup.filter(pl.col("contig") == contig), contig_sequence.sequence, model, motif,
                    low_meth_threshold=self.low_meth_threshold, high_meth_threshold=self.high_meth_threshold
                )
            return model._alpha, model._beta
        n_meth = self._count(("+", True), stripped) + self._count(("-", True), reverse)
        n_nonmeth = self._count(("+", False), stripped) + self._count(("-", False), reverse)
        return n_meth, n_nonmeth

    def model(self, motif: Motif, model=None):
        """Same result as motif_model_bin(pileup, contigs, motif, model, low, high)."""
        model = BetaBernoulliModel() if model is None else model
        model.update(*self.counts(motif))
        return model

def motif_model_bin(
        pileup,
        contigs: dict[str, DNAsequence],
        motif: str,
        model,
        low_meth_threshold,
        high_meth_threshold,
):
    for contig, contig_sequence in contigs.items():
        pileup_contig = pileup.filter(pl.col("contig") == contig)
        model = motif_model_contig(
            pileup_contig,
            contig_sequence.sequence,
            model,
            motif,
            low_meth_threshold=low_meth_threshold,
            high_meth_threshold=high_meth_threshold
        )
    return model

def motif_model_contig(
        pileup, 
        contig: str, 
        prior: BetaBernoulliModel,
        motif, 
        low_meth_threshold=0.3,
        high_meth_threshold=0.7,
        save_motif_positions=False,
    ):
    """
    Get the posterior for a single motif. Uses number of methylated motifs as methylation count.

    Parameters:
    - pileup (Pileup): The pileup to be processed.
    - contig (str): The contig to be processed.
    - motif (str): The motif to be processed.
    - motif_meth_index (int): The index of the methylation site in the motif.

    Returns:
    - BetaBernoulliModel: The model for the candidate.
    """
    # Strip motif of periods    # Strip motif of periods
    motif_stripped = motif.new_stripped_motif()
    pileup_high_meth = pileup.filter(pl.col("fraction_mod") >= high_meth_threshold)
    pileup_low_meth = pileup.filter(pl.col("fraction_mod") <= low_meth_threshold)

    meth_positions_fwd = pileup_high_meth.filter(pl.col("strand") == "+")["position"].to_numpy()
    non_meth_positions_fwd = pileup_low_meth.filter(pl.col("strand") == "+")["position"].to_numpy()
    meth_positions_rev = pileup_high_meth.filter(pl.col("strand") == "-")["position"].to_numpy()
    non_meth_positions_rev = pileup_low_meth.filter(pl.col("strand") == "-")["position"].to_numpy()

    index_meth_fwd, index_nonmeth_fwd = methylated_motif_occourances(motif_stripped, contig, meth_positions_fwd, non_meth_positions_fwd)
    index_meth_rev, index_nonmeth_rev = methylated_motif_occourances(motif_stripped.reverse_compliment(), contig, meth_positions_rev, non_meth_positions_rev)

    model = prior
    model.update(len(index_meth_fwd) + len(index_meth_rev), len(index_nonmeth_fwd) + len(index_nonmeth_rev))
    
    if save_motif_positions:
        motif_data = {
            'index_meth_fwd': index_meth_fwd,
            'index_nonmeth_fwd': index_nonmeth_fwd,
            'index_meth_rev': index_meth_rev,
            'index_nonmeth_rev': index_nonmeth_rev
        }
        return model, motif_data
    else:
        return model


def methylated_reads_counts(
        pileup: pl.DataFrame, 
        sequence: str, 
        motif: Motif
    ) -> tuple:
    """
    Get the number of methylated and non-methylated reads for a motif.

    Parameters:
    - pileup (Pileup): The pileup to be processed. 
    - sequence (str): The sequence to be processed.
    - motif (Motif): The motif to be processed.

    Returns:
    - tuple: A tuple of two integers, the first containing the number of methylated reads, the second containing the number of non-methylated reads.
    """
    motif_index = subseq_indices(motif.string, sequence) + motif.mod_position

    pileup = pileup.filter(pl.col("position").is_in(motif_index))
    counts = pileup \
        .select(["n_mod", "n_nonmod"]) \
        .sum()
    n_mod = counts.get_column("n_mod")[0]
    n_nonmod = counts.get_column("n_nonmod")[0]
    return n_mod, n_nonmod

def predictive_evaluation_score(next_model, current_model) -> float:
    """
    Calculate the score for a motif based on its predictive model.

    Parameters:
        next_model: The model of the next motif.
        current_model: The model of the current motif.

    Returns:
        float: The calculated score.
    """
    extra_positive, extra_negative = current_model._alpha - next_model._alpha, current_model._beta - next_model._beta
    model_extra = BetaBernoulliModel()
    model_extra.update(extra_positive, extra_negative)

    ppcp_next = next_model.posterior_predictive_per_obs(next_model._alpha, next_model._beta)
    ppcp_extra = next_model.posterior_predictive_per_obs(extra_positive, extra_negative)
    mean_ratio = next_model.mean()  / current_model.mean()

    return mean_ratio * (ppcp_next - ppcp_extra)


def get_parent_scores(
        motif: Motif,
        pileup: pl.DataFrame,
        contigs: dict[str, DNAsequence],
        low_meth_threshold: float,
        high_meth_threshold: float,
        site_index: Optional["MotifSiteIndex"] = None,
    ):
    """
    Get the scores for all parent motifs of a motif.
    Parameters:
    - motif (Motif): The motif to be processed.
    - pileup (Pileup): The pileup to be processed.
    - contigs (dict): Dictionary of contig sequences.
    - low_meth_threshold (float): The threshold for low methylation.
    - high_meth_threshold (float): The threshold for high methylation.
    Returns:
    - dict: A dictionary of parent motifs and their scores.
    """
    if site_index is None:
        site_index = MotifSiteIndex(pileup, contigs, low_meth_threshold, high_meth_threshold)
    motif_model = site_index.model(motif)
    motif_split = motif.split()

    parent_motifs = {}
    for i, base in enumerate(motif_split):
        if i == motif.mod_position:
            continue
        if (base != ".") and (base != "N"):
            new_motif = motif_split.copy()
            new_motif[i] = "."
            parent_motif = Motif("".join(new_motif), motif.mod_position)
            parent_model = site_index.model(parent_motif)
            score = predictive_evaluation_score(motif_model, parent_model)
            parent_motifs[parent_motif] = dict(
                motif_position=i,
                parent_model=parent_model,
                child_model=motif_model,
                score=score
            )
    return parent_motifs
    

def merge_motifs_in_df(motif_df, pileup, assembly, contig_bin, low_meth_threshold=0.3, high_meth_threshold=0.7, merge_threshold = 0.5,
                       site_index=None):
    """
    site_index: MotifSiteIndex of the bin's pileup, used when motif_df holds a single bin and mod type
    """
    new_df = []
    groups = motif_df.group_by("reference", "mod_type")
    single_group = motif_df.select("reference", "mod_type").n_unique() == 1
    for (bin_name, mod_type), df in groups:
        contigs = contig_bin[bin_name]
        log.debug(f"Merging motifs of bin: {bin_name}, modtype: {mod_type}")
        
        # Create dictionary of bin contigs to sequence
        bin_sequences = {contig: assembly[contig] for contig in contigs}
        if site_index is not None and single_group:
            group_index = site_index
        else:
            group_index = MotifSiteIndex(
                pileup.filter((col("contig").is_in(contigs)) & (col("mod_type") == mod_type)),
                bin_sequences, low_meth_threshold, high_meth_threshold
            )

          
        # Get list of motifs
        motif_seq = df["motif"].to_list()
        motif_pos = df["mod_position"].to_list()
        motifs = [nm.motif.Motif(seq, pos) for seq, pos in zip(motif_seq, motif_pos)]

        # Merge motifs
        merged_motifs = nm.motif.merge_motifs(motifs)
        all_merged_motifs = []
        all_premerge_motifs = []

        for _, (merged_motif, premerge_motifs, premerge_motif_variants, new_motif_variants) in merged_motifs.items():
            if len(new_motif_variants) == 0:
                log.debug(f"Merging motifs {premerge_motifs} to {merged_motif}, as no new motif variants were generated")
                all_merged_motifs.append(merged_motif)
                all_premerge_motifs.extend(premerge_motifs)
                continue
            merge_model = group_index.model(merged_motif)

            premerge_motif_variants_model = BetaBernoulliModel()
            for premerge_motif_variant in premerge_motif_variants:
                premerge_motif_variants_model = group_index.model(premerge_motif_variant, premerge_motif_variants_model)
            merge_score = predictive_evaluation_score(premerge_motif_variants_model, merge_model)
            if merge_score < merge_threshold:
                log.debug(f"Merging motifs {premerge_motifs} to {merged_motif}, merge score: {merge_score}")
                all_merged_motifs.append(merged_motif)
                all_premerge_motifs.extend(premerge_motifs)
            else:
                log.debug(f"Not merging motifs {premerge_motifs} to {merged_motif}, merge score: {merge_score}")
                continue

            
        log.info(f"Motif merge for bin {bin_name} modtype {mod_type} complete")
        # Create a new dataframe with the non merged motifs
        if  len(all_premerge_motifs) == 0:
            # No motifs were merged
            new_df.append(df)
            log.info(f"No motifs were merged for bin {bin_name} modtype {mod_type}")
            continue
        log.info(f"New merged motifs: {all_merged_motifs}")
        single_df = df.filter(col("motif").is_in(all_premerge_motifs).not_())
        new_df.append(single_df)
        merged_df = []
        for motif in all_merged_motifs:
            merged_model = group_index.model(motif)
            merged_model_parents = get_parent_scores(
                motif,
                pileup.filter((col("contig").is_in(contigs)) & (col("mod_type") == mod_type)),
                bin_sequences,
                low_meth_threshold,
                high_meth_threshold,
                site_index=group_index
            )
            if len(merged_model_parents) > 0:   
                merged_score = np.mean([d["score"] for d in merged_model_parents.values()])
            else:
                merged_score = -1
            merge_motif_df = nm.motif.MotifSearchResult(pl.DataFrame({
                "motif": motif.string,
                "score": merged_score,
                "mod_position": motif.mod_position,
                "reference": bin_name,
                "mod_type": mod_type,
                "model": merged_model
            }))
            merged_df.append(merge_motif_df)

        merged_df_con = pl.concat(merged_df)
        merged_df_con = nm.motif.MotifSearchResult(merged_df_con)
        new_df.append(merged_df_con)
    new_df = pl.concat(new_df)
    new_df = nm.motif.MotifSearchResult(new_df)
    return new_df
