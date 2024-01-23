import polars as pl
import multiprocessing
import os
import warnings
from multiprocessing import get_context
import logging as log
import nanomotif as nm
from nanomotif import candidate, evaluate
from nanomotif.parallel import update_progress_bar
from nanomotif.logger import configure_logger

def score_contig_with_motifs(contig, mod_type, pileup, sequence, motifs):
    log.info(f"Scoring motifs in contig {contig}, modtype {mod_type}")

    # Ensuring motifs are in regex for searching
    motifs = motifs.with_columns([
            pl.col("motif") \
                .map_elements(lambda x: nm.seq.regex_to_iupac(x)) \
                .map_elements(lambda x: nm.seq.iupac_to_regex(x)) \
                .alias("motif")
        ])
    
    # Getting motifs for scoring
    motifs = motifs.select("motif", "mod_position", "mod_type").unique()
    motifs = motifs.filter(pl.col("mod_type") == mod_type)
    motifs_str = motifs.get_column("motif").to_list()
    positions = motifs.get_column("mod_position").to_list()
    types = motifs.get_column("mod_type").to_list()

    # Scoring motifs
    models = []
    for motif, position, mod_type in zip(motifs_str, positions, types):
        log.debug(f"Scoring motif {motif}, {position}")
        model = evaluate.motif_model_contig(
            pileup,
            sequence, 
            candidate.Motif(motif, position)
        )
        models.append(model)

    # Creating a dataframe of the results
    result = pl.DataFrame({
            "motif": motifs_str,
            "mod_position": positions,
            "mod_type": types,
            "model": models
        })
    result = result.with_columns([
            pl.col("model").map_elements(lambda x: x._alpha).alias("n_mod"),
            pl.col("model").map_elements(lambda x: x._beta).alias("n_nomod"),
            pl.col("motif").map_elements(lambda x: nm.seq.regex_to_iupac(x)).alias("motif"),
            pl.lit(contig).alias("contig")
        ]).drop("model")
    return result

def worker_function(args, counter, lock, sequence, motifs, log_dir, verbose, seed):
    """
    Process a single subpileup for one contig and one modtype

    Parameters:
    - args (tuple): The arguments to the function: contig, modtype, subpileup
    - counter (multiprocessing.Value): The progress counter
    - lock (multiprocessing.Lock): The lock for the progress counter
    - assembly (Assembly): The assembly to be processed.
    - min_kl_divergence (float): Early stopping criteria, if max KL-divergence falls below, stops building motif.
    - padding (int): The padding to use for the motif.
    """
    warnings.filterwarnings("ignore")
    nm.seed.set_seed(seed)
    contig, modtype, subpileup = args
    
    process_id = os.getpid()
    log_file = f"{log_dir}/score-motifs.{process_id}.log"
    configure_logger(log_file, verbose=verbose)

    try:
        result = score_contig_with_motifs(contig, modtype, subpileup, sequence, motifs)
        with lock:
            counter.value += 1
        return result
    except:
        with lock:
            counter.value += 1
        return None

def score_sample_parallel(
        assembly, pileup, motifs,
        threads = 2,
        threshold_methylation_general = 0.5,
        threshold_valid_coverage = 1,
        verbose = False,
        log_dir = None,
        seed = None
    ):
    """
    Process a single sample
    
    Parameters:
    - assembly (Assembly): The assembly to be processed.
    - pileup (Pileup): The pileup to be processed.
    - max_candidate_size (int): The maximum size of the candidate motifs.
    - min_read_methylation_fraction (float): The minimum fraction of reads that must be methylated for a position to be considered methylated.
    - threshold_valid_coverage (int): The minimum number of reads that must cover a position for it to be considered valid.
    - min_kl_divergence (float): Early stopping criteria, if max KL-divergence falls below, stops building motif.
    - min_cdf_score (float): Minimum score of 1 - cdf(cdf_position) for a motif to be considered valid.
    - cdf_position (float): The position to evaluate the cdf at.
    - min_motif_frequency (int): Used to get minimum number of sequences to evaluate motif at.
    """
    assert pileup is not None, "Pileup is None"
    assert assembly is not None, "Assembly is None"
    assert threshold_methylation_general >= 0 and threshold_methylation_general <= 1, "threshold_methylation_general must be between 0 and 1"
    assert threshold_valid_coverage >= 0, "min_valid_coverage must be greater than 0"

    # Filter pileup
    pileup = pileup \
            .filter(pl.col("Nvalid_cov") > threshold_valid_coverage) \
            .filter(pl.col("fraction_mod") >= threshold_methylation_general) \
            .sort("contig") 
    
    # Create a list of tasks (TODO: not have a list of all data)
    tasks = [(contig, modtype, subpileup) for (contig, modtype), subpileup in pileup.groupby(["contig", "mod_type"])]

    # Create a progress manager
    manager = multiprocessing.Manager()
    counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Create a pool of workers
    pool = get_context("spawn").Pool(processes=threads)

    # Create a process for the progress bar
    progress_bar_process = multiprocessing.Process(target=update_progress_bar, args=(counter, len(tasks), True))
    progress_bar_process.start()

    # Put them workers to work
    results = pool.starmap(worker_function, [(
        task, 
        counter, lock, 
        assembly.assembly[task[0]].sequence, motifs,
        log_dir, verbose, seed
        ) for task in tasks])
    results = [result for result in results if result is not None]

    # Close the pool
    pool.close()
    pool.join()

    # Close the progress bar
    progress_bar_process.join()

    if len(results) == 0:
        return None
    motifs_all_scored = pl.concat(results)
    return motifs_all_scored



def score_all_contigs_all_motifs(pileup, assm, motifs, frac_mod=0.7):
    # Score all contigs with identified motifs
    all_motifs_scored = {
        "contig":[],
        "type":[],
        "motif":[],
        "mod_position":[],
        "model":[],
        "mean":[]
    }
    for contig in pileup.pileup.get_column("contig").unique().to_list():
        for row in motifs.select("motif", "mod_position", "mod_type").unique().rows():
            model = evaluate.motif_model_contig(
                pileup.pileup.filter(pl.col("Nvalid_cov") > 1).filter(pl.col("fraction_mod") > frac_mod),
                assm.assembly[contig].sequence, 
                candidate.Motif(row[0], row[1])
            )
            all_motifs_scored["contig"].append(contig)
            all_motifs_scored["motif"].append(row[0])
            all_motifs_scored["model"].append(model)
            all_motifs_scored["mean"].append(model.mean())
            all_motifs_scored["type"].append(row[2])
            all_motifs_scored["mod_position"].append(row[1])

    scored_all = pl.DataFrame(all_motifs_scored).with_columns([
        pl.col("model").map_elements(lambda x: x._alpha).alias("alpha"),
        pl.col("model").map_elements(lambda x: x._beta).alias("beta")
    ])
    return scored_all
