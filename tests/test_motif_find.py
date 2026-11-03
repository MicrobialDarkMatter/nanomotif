from nanomotif.constants import *
from nanomotif.motif import Motif
from nanomotif.seq import DNAsequence, DNAarray
import nanomotif as nm
import pytest
from hypothesis import given, strategies as st
import itertools
import pandas as pd
import polars as pl
import polars.testing
import numpy as np


class TestMethylatedMotifOccurrences:
    def test_known_motif_sequence(self):
        motif = nm.motif.Motif('ACG', 0)
        sequence = 'TACGGACGCCACG'
        methylated_positions = np.array([1, 5])
        non_methylated_positions = np.array([10])

        expected_methylated = np.array([1, 5])
        expected_nonmethylated = np.array([10])
        
        result = nm.find_motifs_bin.methylated_motif_occourances(motif, sequence, methylated_positions, non_methylated_positions)
        np.testing.assert_array_equal(result[0], expected_methylated)
        np.testing.assert_array_equal(result[1], expected_nonmethylated)

    def test_known_motif_sequence_no_methylated(self):
        motif = nm.motif.Motif('ACG', 0)
        sequence = 'TACGGACGCCACG'
        methylated_positions = np.array([])
        non_methylated_positions = np.array([1, 10])

        expected_methylated = np.array([])
        expected_nonmethylated = np.array([1, 10])
        
        result = nm.find_motifs_bin.methylated_motif_occourances(motif, sequence, methylated_positions, non_methylated_positions)
        np.testing.assert_array_equal(result[0], expected_methylated)
        np.testing.assert_array_equal(result[1], expected_nonmethylated)


def test_get_motif_parental_relationships():
    motif1 = Motif("ACGT", 0)
    motif2 = Motif("ACG", 0)
    motif3 = Motif("CGT", 1)
    motif4 = Motif("ACGTG", 0)
    motif5 = Motif("TGCA", 1)
    
    motifs = [motif1, motif2, motif3, motif4, motif5]
    expected_relationships = set([
        (motif2, motif1),
        (motif2, motif4),
        (motif3, motif1),
        (motif3, motif4),
        (motif1, motif4)
    ])
    relationships = set(nm.postprocess.get_motif_parental_relationship(motifs))
    assert relationships == expected_relationships

    
    motif6 = Motif("GGCA[AT]", 2)
    motif7 = Motif("GGCAAT", 2)
    motif8 = Motif("GGCAAT", 4)
    motif9 = Motif("AATTT", 0)
    motif10 = Motif("AATTT", 1)
    motif11 = Motif("AATTTT", 0)
    motifs = [motif6, motif7, motif8, motif9, motif10, motif11]
    expected_relationships = set([
        (motif6, motif7),
        (motif6, motif8),
        (motif9, motif11),
        (motif10, motif11)
    ])
    relationships = set(nm.postprocess.get_motif_parental_relationship(motifs))
    assert relationships == expected_relationships

    motif12 = Motif("TTAAGGAG", 6)
    motif13 = Motif("TTAA", 3)
    motifs = [motif12, motif13]
    expected_relationships = set([
        (motif13, motif12)
    ])
    relationships = set(nm.postprocess.get_motif_parental_relationship(motifs))
    assert relationships == expected_relationships

def _bin1_a_inputs():
    pileup = nm.datasets.geobacillus_plasmids_pileup()
    pileup = nm.dataload.filter_pileup(pileup)
    pileup = nm.dataload.filter_pileup_minimummod_frequency(pileup)
    pileup = nm.dataload.filter_pileup_adjacency_filter(pileup)
    pileup = pileup.filter(pl.col("mod_type") == "a")
    assembly = nm.datasets.geobacillus_plasmids_assembly()
    sequences = {c: assembly[c] for c in pileup.get_column("contig").unique().to_list()}
    return pileup, sequences


def test_contig_without_usable_windows_does_not_discard_bin(tmp_path):
    """A contig whose methylated positions are all too close to its ends to extract a
    search window must be skipped, not end the search for the whole bin."""
    pileup, sequences = _bin1_a_inputs()
    rng = np.random.default_rng(1)
    edge_seq = "".join(rng.choice(list("ACGT"), 400))
    sequences["edge_contig"] = DNAsequence(edge_seq)
    edge_positions = [i for i in range(15) if edge_seq[i] == "A"]
    edge_rows = pl.DataFrame({
        "contig": "edge_contig",
        "position": edge_positions,
        "mod_type": "a",
        "strand": "+",
        "fraction_mod": 0.95,
        "Nvalid_cov": 20,
    }).cast(pileup.schema)
    pileup = pl.concat([pileup, edge_rows])

    result = nm.find_motifs_bin.find_best_candidates(
        pileup, sequences, "a", "bin1", output_dir=str(tmp_path),
        low_meth_threshold=0.3, high_meth_threshold=0.7, padding=20, min_kl=0.05
    )
    assert result is not None
    _, best_candidates = result
    assert len(best_candidates) > 0


def test_process_subpileup_passes_thresholds_to_merge(tmp_path, monkeypatch):
    pileup, sequences = _bin1_a_inputs()
    seen = {}
    real_merge = nm.find_motifs_bin.merge_motifs_in_df

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return real_merge(*args, **kwargs)

    monkeypatch.setattr(nm.find_motifs_bin, "merge_motifs_in_df", spy)
    nm.find_motifs_bin.process_subpileup(
        {"bin1": list(sequences)}, "a", pileup, sequences,
        min_kl_divergence=0.05, padding=20,
        low_meth_threshold=0.2, high_meth_threshold=0.8, score_threshold=1.5,
        output_dir=str(tmp_path)
    )
    assert seen.get("low_meth_threshold") == 0.2
    assert seen.get("high_meth_threshold") == 0.8


def test_plain_strategy_task_count_matches_generator():
    """The progress bar total must equal the number of tasks the generator accounts for
    (yielded or skipped), otherwise it reports completion too early."""
    import threading
    import types

    class Counter:
        value = 0

    pileup, sequences = _bin1_a_inputs()
    strategy = nm.find_motifs_bin.PileupStrategyPlain(
        partitioned_assembly={"bin1": sequences, "bin2": sequences},
        partitioned_pileup={("bin1", "a"): pileup},
        bin_contig={"bin1": list(sequences), "bin2": list(sequences)},
        config=types.SimpleNamespace(assembly_path=None),
    )
    counter = Counter()
    n_yielded = len(list(strategy.task_generator(counter, threading.Lock())))
    assert n_yielded + counter.value == strategy.n_tasks()


def test_motif_site_index_matches_regex_counts():
    """MotifSiteIndex must give the counts of motif_model_bin (regex over the contigs)."""
    import random
    from nanomotif.model import BetaBernoulliModel
    pileup, sequences = _bin1_a_inputs()
    index = nm.find_motifs_bin.MotifSiteIndex(pileup, sequences, 0.3, 0.7, half_width=20)
    random.seed(0)
    motifs = [Motif("GATC", 1), Motif("..GAATTC..", 4), Motif("A", 0), Motif("[AG]GATC[CT]", 2),
              Motif("A" + "." * 25 + "GATC", 0)]  # the last does not fit in the window
    for _ in range(40):
        left = "".join(random.choice(["A", "C", "G", "T", ".", "[AG]", "N"]) for _ in range(random.randint(0, 5)))
        right = "".join(random.choice(["A", "C", "G", "T", ".", "[CT]"]) for _ in range(random.randint(0, 5)))
        motifs.append(Motif("." * 20 + left + "A" + right + "." * 20, 20 + len(Motif(left, 0).split())))
    for motif in motifs:
        expected = nm.find_motifs_bin.motif_model_bin(pileup, sequences, motif, BetaBernoulliModel(), 0.3, 0.7)
        got = index.model(motif)
        assert (got._alpha, got._beta) == (expected._alpha, expected._beta), motif


def test_popcount_fallback():
    words = np.array([0, 1, 2**63 + 5, 2**64 - 1], dtype=np.uint64)
    expected = sum(bin(int(w)).count("1") for w in words)
    assert nm.find_motifs_bin._popcount(words) == expected
    assert int(nm.find_motifs_bin._POPCOUNT_8[words.view(np.uint8)].sum()) == expected
