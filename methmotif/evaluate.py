import numpy as np
from itertools import product
import re
import polars as pl
from Bio import SeqIO
from scipy.stats import beta



class Pileup():
    """
    Class for loading pileup files
    """
    def __init__(self, path):
        self.path = path
        self.pileup = load_pileup(path)

    def get_contigs(self):
        """
        Get contigs in pileup
        """
        return self.pileup["contig"].unique()

    def get_mod_types(self):
        """
        Get modification types in pileup
        """
        return self.pileup["mod_type"].unique()
    
    def get_strands(self):
        """
        Get strands in pileup
        """
        return self.pileup["strand"].unique()
    
    def get_methylated_positions(self, threshold = 0.8):
        """
        Get methylated positions in pileup
        """
        return self.pileup.filter(pl.col("fraction_mod") >= threshold)


class BetaBernoulliModel():
    def __init__(self, alpha = 1, beta = 1):
        self.alpha = alpha
        self.beta = beta
        self._alpha = alpha
        self._beta = beta

    def update(self, n_positives, n_negatives):
        self._alpha += n_positives
        self._beta += n_negatives

    def reset(self):
        self._alpha = self.alpha
        self._beta = self.beta

    def sample(self):
        return np.random.beta(self._alpha, self._beta)

    def mean(self):
        return self._alpha / (self._alpha + self._beta)

    def variance(self):
        a = self._alpha
        b = self._beta
        return (a * b) / ((a + b)**2 * (a + b + 1))
    
    def standard_deviation(self):
        return np.sqrt(self.variance())
    
    def cdf(self, x):
        return beta.cdf(x, self._alpha, self._beta)

    def __str__(self):
        return f'BetaBernoulliModel(alpha={self._alpha}, beta={self._beta})'

    def __repr__(self):
        return str(self)


def generate_kmers(k, alphabet):
    # All possible k-mers of size 1 to k
    return [''.join(p) for i in range(1, k + 1) for p in product(alphabet, repeat=i)]
    


class MotifCandidates():
    def __init__(self, motifs, modification_position):
        self.candidates = motifs
        self.modification_position = modification_position

    def __iter__(self):
        return iter(zip(self.candidates, self.modification_position))

    def __len__(self):
        return len(self.candidates)

def load_pileup(path):
    """
    Load pileup file from path to pileup.bed output of modkit pileup
    """
    pileup = pl.read_csv(path, separator = "\t", has_header = False)
    pileup = pileup.select(["column_1", "column_2","column_4", "column_6", "column_11", "column_10"]) \
        .rename({"column_1":"contig", "column_2": "position", "column_4": "mod_type", "column_6": "strand", "column_11": "fraction_mod", "column_10":"Nvalid_cov"}) \
        .with_columns(pl.col("fraction_mod") / 100) \
        .sort("position")
    return pileup

def load_fasta(path, trim_names=False, trim_character=" ") -> dict:
    """
    Reads a fasta file and returns a dictionary with the contig names as 
    keys and the sequences as values
    """
    with open(path, 'r') as f:
        lines = f.readlines()
    data = {}
    active_sequence_name = "no_header"
    for line in lines:
        line = line.rstrip()
        if line[0] == '>':
            active_sequence_name = line[1:]
            if trim_names:
                active_sequence_name = active_sequence_name.split(trim_character)[0]
            if active_sequence_name not in data:
                data[active_sequence_name] = ''
        else:
            data[active_sequence_name] += line
    return data

class Assembly():
    def __init__(self, path):
        self.path = path
        self.assembly = load_fasta(path)
    
    def __getitem__(self, key):
        return self.assembly[key]
    def find_index_motif(self, contig, motif, offset=0):
        """
        Find occourance indeces of a motif of a single contig
        """
        index = [match.start() + offset for match in re.finditer(motif, self.assembly[contig])]
        return index


def methylated_motif_occourances(motif, assembly, contig_name, methylation_positions, offset):
    """
    Get occourances of a motif in a contig
    """
    motif_index = assembly.find_index_motif(contig_name, motif, offset)
    methylated_occurences = np.intersect1d(methylation_positions, motif_index)
    nonmethylated_occurences = np.setdiff1d(motif_index, methylation_positions)
    return methylated_occurences, nonmethylated_occurences

def score_candidates(methylation_pileup, assembly, contig, candidates):
    """
    Score candidates for a single contig
    """
    methylation_positions_fwd = methylation_pileup.filter(pl.col("strand") == "+")["position"].to_numpy()
    methylation_positions_rev = methylation_pileup.filter(pl.col("strand") == "-")["position"].to_numpy()

    score_canidates = {
        "motif": [],
        "motif_mod_index":[],
        "posterior": []
    }
    for candidate, motif_mod_index in candidates:
        index_meth_fwd, index_nonmeth_fwd = methylated_motif_occourances(candidate, assembly, contig, methylation_positions_fwd, offset = motif_mod_index)
        index_meth_rev, index_nonmeth_rev = methylated_motif_occourances(reverse_complement(candidate), assembly, contig, methylation_positions_rev, offset = len(candidate) - 1 - motif_mod_index)

        model = BetaBernoulliModel()
        model.update(len(index_meth_fwd) + len(index_meth_rev), len(index_nonmeth_fwd) + len(index_nonmeth_rev))

        score_canidates["motif"].append(candidate)
        score_canidates["motif_mod_index"].append(motif_mod_index)
        score_canidates["posterior"].append(model)
    return pl.DataFrame(score_canidates)



def reverse_complement(seq):
    """
    Reverse complement a sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    return ''.join([complement[base] for base in reversed(seq)])

if __name__ == "__main__":
    pass