from scipy.stats import entropy
import numpy as np
import warnings
import re
import random
random.seed(2403)
import nanomotif.utils as utils
from nanomotif.constants import *


class Assembly():
    def __init__(self, assembly: dict):
        self.assembly = {name:DNAsequence(seq) for name, seq in assembly.items()}
    
    def __getitem__(self, key):
        return self.assembly[key]
    
    def __repr__(self):
        return f"Assembly with {len(self.assembly)} contigs"

class DNAsequence:
    bases = ["A", "T", "G", "C"]
    iupac = ["A", "T", "G", "C", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
    iupac_dict = {
        "A": ["A"], "T": ["T"], "C": ["C"], "G": ["G"],
        "R": ["A", "G"], "Y": ["C", "T"], 
        "S": ["G", "C"], "W": ["A", "T"], 
        "K": ["G", "T"], "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "T", "C", "G"]
    }
    complement = {
        "A": "T", "T": "A", "G": "C", "C": "G", 
        "N": "N", "R": "Y", "Y": "R", "S": "S", 
        "W": "W", "K": "M", "M": "K", "B": "V", 
        "D": "H", "H": "D", "V": "B"
    }
    base_to_vector = {
        "A": [1, 0, 0, 0],
        "T": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "C": [0, 0, 0, 1],
        "N": [1, 1, 1, 1],
        ".": [1, 1, 1, 1]
    }

    base_int = {"A":1, "T":2, "G":3, "C":4, "N":5, "R":6, "Y":7, "S":8, "W":9, "K":10, "M":11, "B":12, "D":13, "H":14, "V":15}
    int_base = {i: n for n, i in base_int.items()}

    def __init__(self, sequence):
        self._check_sequence(sequence)
        sequence = sequence.upper()
        self.sequence = sequence
    
    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        self.sequence_alphabet = list(set(value))
        self._check_sequence(value)
        self._sequence = value

    def _check_sequence(self, sequence):
        assert isinstance(sequence, str), "DNA sequence must be a str"
        assert len(sequence) > 0, "DNA sequence must not be empty"
        assert all(letter in self.iupac for letter in list(set(sequence))), f"DNA sequence must be a nucleotide sequence of {''.join(self.iupac)}"

    
    def __getitem__(self, key):
        return self.sequence[key]

    def __len__(self):
        return len(self.sequence)

    def __iter__(self):
        return iter(self.sequence)
    
    def __repr__(self):
        return f"DNAsequence() | Number of sequence: {len(self)} | Unique Bases: {self.sequence_alphabet}"

    def __add__(self, other):
        if isinstance(other, DNAsequence):
            return DNAsequence(self.sequence + other.sequence)
        elif isinstance(other, str):
            return DNAsequence(self.sequence + other)
        else:
            return NotImplemented

    def __mul__(self, other):
        try:
            return DNAsequence(self.sequence * other)
        except:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, DNAsequence):
            return self.sequence == other.sequence
        elif isinstance(other, str):
            return self.sequence == other
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, DNAsequence):
            return self.sequence != other.sequence
        elif isinstance(other, str):
            return self.sequence != other
        else:
            return NotImplemented

    def __lt__(self, other):
        if isinstance(other, int):
            return len(self) < other
        elif isinstance(other, (DNAsequence, str)):
            return len(self) < len(other)
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, int):
            return len(self) <= other
        elif isinstance(other, (DNAsequence, str)):
            return len(self) <= len(other)
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, int):
            return len(self) > other
        elif isinstance(other, (DNAsequence, str)):
            return len(self) > len(other)
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, int):
            return len(self) >= other
        elif isinstance(other, (DNAsequence, str)):
            return len(self) >= len(other)
        else:
            return NotImplemented
    
    def sample_at_index(self, index: int, padding: int):
        """
        Randomly sample sequence around index in DNAsequence.

        Parameters
        ----------
        index : int
            Index to sample at.
        padding : int
            Number of position to include on each side

        Returns
        -------
        sequences : DNAsequence
            Sequence of length 2 * padding + 1
        """
        if index < padding:
            raise ValueError(f"Index {index} is too small for padding {padding}")
        if index > len(self) - padding:
            raise ValueError(f"Index {index} is too large for padding {padding}")
        return DNAsequence(self.sequence[(index - padding):(index + padding + 1)])

    def sample_at_indices(self, indices: list, padding: int):
        """
        Randomly sample sequence around indices in DNAsequence.

        Parameters
        ----------
        indices : list
            Indices to sample at.
        padding : int
            Number of position to include on each side

        Returns
        -------
        sequences : EqualLengthDNASet
            Sequences of length 2 * padding + 1
        """
        indices_checked = [index for index in indices if (index > padding) and (index < (len(self) - padding))]
        return EqualLengthDNASet([self.sample_at_index(index, padding) for index in indices_checked])

    def sample_subsequence(self, length):
        """Randomly sample a subsequence of a given length"""
        if length > len(self):
            raise ValueError(f"Cannot sample a subsequence of length {length} from a sequence of length {len(self)}")
        start_index = random.randint(0, len(self) - length)
        return DNAsequence(self.sequence[start_index:start_index + length])

    def sample_n_subsequences(self, length: int, n: int):
        """Randomly sample n sequences of a given length"""
        return EqualLengthDNASet([self.sample_subsequence(length=length) for _ in range(n)])

    def find_subsequence(self, subsequence: str):
        """
        Find all occourances of subseqeunce in DNAsequence
        """
        p = len(subsequence)
        t = len(self.sequence)
        subsequence_hash = hash(subsequence)
        sequence_hash = hash(self.sequence[:p])
        indices = []
        for i in range(t - p + 1):
            if subsequence_hash == sequence_hash:
                if subsequence == self.sequence[i:i+p]:
                    indices.append(i)
            if i < t - p:
                sequence_hash = hash(self.sequence[i+1:i+p+1])
        
        return np.array(indices)

    def sample_at_subsequence(self, padding: int, subsequence: str):
        """
        Randomly sample sequence around subsequence in DNAsequence.

        Parameters
        ----------
        padding : int
            Number of position to include on each side
        subsequence : str
            Subsequence to sample at.

        Returns
        -------
        sequences : DNAsequence
            Sequence of length 2 * padding + len(subsequence)
        """
        subsequence_index = self.find_subsequence(subsequence)
        subsequence_index = subsequence_index[subsequence_index >= padding]
        subsequence_index = subsequence_index[subsequence_index < len(self.sequence) - padding]
        if subsequence_index.shape[0] == 0:
            warnings.warn("Not occourances found " + subsequence + " found")
        else:
            sample_index = np.random.choice(subsequence_index, 1)[0]
            return DNAsequence(self.sequence[(sample_index - padding):(sample_index + padding + len(subsequence))])

    def sample_n_at_subsequence(self, padding: int, subsequence: str, n: int):
        """
        Randomly sample sequence around subsequence in DNAsequence.

        Parameters
        ----------
        padding : int
            Number of position to include on each side
        subsequence : str
            Subsequence to sample at.
        n: int
            number of sequences to sample (all unique)

        Returns
        -------
        sequences : EqualLengthDNASet
            Sequences of length 2 * padding + len(subsequence)
        """
        subsequence_index = self.find_subsequence(subsequence)
        subsequence_index = subsequence_index[subsequence_index >= padding]
        subsequence_index = subsequence_index[subsequence_index < (len(self.sequence) - padding)]
        if subsequence_index.shape[0] < n:
            warnings.warn("Only found " + str(subsequence_index.shape[0]) + " occourances, sampling " + str(subsequence_index.shape[0]) + " sequences")
            n = subsequence_index.shape[0]
        sample_index = np.random.choice(subsequence_index, n, replace=False)
        sampled_sequences = [DNAsequence(self.sequence[(i - padding):(i + padding + len(subsequence))]) for i in sample_index]
        return EqualLengthDNASet(sampled_sequences)

    def reverse_complement(self):
        """Return the reverse complement of the sequence."""
        return DNAsequence("".join([self.complement[base] for base in reversed(self.sequence)]))

    def count(self, subsequence):
        """Count the number of occurrences of a subsequence in the sequence."""
        return self.sequence.count(subsequence)
    
    def gc(self) -> list:
        '''
        Returns the GC content of a sequence

        Returns
        -------
        Float
            GC content of sequence

        Examples
        --------
        >>> DNAsequence("ATCG").gc()
        0.5
        '''
        return (self.count("G") + self.count("C")) / len(self)





class EqualLengthDNASet:
    def __init__(self, sequences: list):
        self._check_sequences(sequences)
        self.sequences = sequences
        self.sequence_length = len(sequences[0])

    @property
    def sequences(self):
        return self._sequences

    @sequences.setter
    def sequences(self, value):
        self._check_sequences(value)
        self._sequences = value
    
    def __add__(self, other):
        if isinstance(other, EqualLengthDNASet):
            return EqualLengthDNASet(self.sequences + other.sequences)
        elif isinstance(other, list):
            return EqualLengthDNASet(self.sequences + other)
        else:
            return NotImplemented

    def _check_sequences(self, sequences):
        assert isinstance(sequences, list), "DNA sequences must be a list"
        assert len(sequences) > 0, "DNA sequences must not be empty"
        assert all(isinstance(seq, DNAsequence) for seq in sequences), "DNA sequences must be a list of strings"
        assert all(len(seq) > 0 for seq in sequences), "DNA sequences must not contain empty strings"
        assert utils.all_lengths_equal(sequences), "All sequences must be of equal length"

    def __getitem__(self, key):
        return EqualLengthDNASet(self.sequences[key])

    def levenshtein_distances(self) -> np.ndarray:
        """
        All vs. all levenshtein distances between sequences.

        Returns
        -------
        np.ndarray
            A matrix of edit distances between all sequences

        Examples
        --------
        >>> EqualLengthDNASet(["ATCG", "GCTA", "TACT", "AGCT"]).levenshtein_distances()
        array([[0, 4, 3, 2],
               [4, 0, 3, 2],
               [3, 3, 0, 2],
               [2, 2, 2, 0]], dtype=int16)
        """
        lengths = [len(s) for s in self]
        assert utils.all_equal(lengths), "All sequences must be of equal length"

        n = len(self)
        distances = np.empty(shape = (n, n), dtype=np.int16)
        for i in range(0, n):
            for j in range(i, n):
                d = editdistance.eval(self[i], self[j])
                distances[i, j] = d
                distances[j, i] = d
        return distances

    def reverse_compliment(self):
        """Return the reverse complement of the sequence."""
        return EqualLengthDNASet([seq.reverse_complement() for seq in self.sequences])

    def pssm(self, pseudocount=0, only_canonical: bool = True):
        """
        Calculate the positional frequency of each nucleotide for all sequences.

        Parameters
        ----------
        pseudocount : float, optional
            Pseudocount to be added to each nucleotide count. Default is 0.5.

        Returns
        -------
        np.ndarray
            A numpy array of positional frequencies.
        """
        if only_canonical:
            bases = self.sequences[0].bases
        else:
            bases = self.sequences[0].iupac

        n = len(self.sequences)

        pssm = []
        for nuc in bases:
            pssm_nuc = []
            for i in range(self.sequence_length):
                count = 0
                for seq in self.sequences:
                    if seq[i] == nuc:
                        count += 1
                pssm_nuc.append((count + pseudocount)/n)
            pssm.append(pssm_nuc)
        return np.array(pssm)

    def kl_divergence(self, pssm):
        """
        Calculate positional Kullback-Liebler divergence from sequences to provided PSSM.
        Parameters
        ----------
        pssm : Position-specific scoring matrix
            PSSM to calculate distances from
        
        Returns
        -------
        np.ndarray
            A numpy array of Kullback-Leibler divergence.
        """
        assert pssm.shape == (4, self.sequence_length), "Shape does not match (should be [4, sequence_length])"
        return entropy(self.pssm(), pssm)

    def get_sequences_with_match(self, pattern: str):
        """
        Remove all sequences that contain a given pattern.

        Parameters
        ----------
        pattern : str
            pattern to remove

        Returns
        -------
        EqualLengthDNASet with sequences containing pattern
        """
        return EqualLengthDNASet([seq for seq in self.sequences if bool(re.match(pattern, seq.sequence))])

    def get_sequences_without_match(self, pattern: str):
        """
        Remove all sequences that contain a given pattern.

        Parameters
        ----------
        pattern : str
            pattern to remove

        Returns
        -------
        EqualLengthDNASet with sequences containing pattern
        """
        sequences = [seq for seq in self.sequences if not bool(re.match(pattern, seq.sequence))]
        if len(sequences) == 0:
            warnings.warn("No sequences left after filtering")
            return None
        return EqualLengthDNASet(sequences)

    def convert_to_DNAarray(self):
        """
        Convert EqualLengthDNASet to DNAarray
        """
        return DNAarray([[DNAsequence.base_to_vector[base] for base in sequence] for sequence in self.sequences])

class DNAarray(np.ndarray):
    """
    A numpy array of DNA sequences
    """
    base_to_vector = {
        "A": [1, 0, 0, 0],
        "T": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "C": [0, 0, 0, 1],
        "N": [1, 1, 1, 1],
        ".": [1, 1, 1, 1]
    }
    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None: return
    
    def filter_sequence_matches(self, sequence, keep_matches=True):
        """
        Filter the array based on sequence matches

        Parameters
        ----------
        sequence : np.ndarray
        keep_matches : bool, optional
            Whether to keep matches or non-matches, by default True

        Returns
        -------
        DNAarray
            The filtered array
        """
        assert isinstance(sequence, np.ndarray), "Sequence must be a numpy array"
        assert sequence.shape == self.shape[1:3], "Sequence must have the same length as sequences in the array"

        if keep_matches:
            result = self[np.all(self <= sequence, axis=(1,2)), :]
        else:
            result = self[np.any(self > sequence, axis=(1,2)), :]
        if result.shape[0] == 0:
            warnings.warn("No sequences left after filtering")
            return None
        return DNAarray(result)
    
    def pssm(self):
        """
        Calculate the position specific scoring matrix for the array. 
        The PSSM is the sum of the array along the first axis (sequences) divided by the number of sequences, 
        e.g. the positional frequency of bases in all of the sequences in the array.

        Returns
        -------
        np.ndarray
            The PSSM
        """
        return self.sum(axis = 0).transpose() / self.shape[0]

def regex_to_iupac(regex):
    # Dictionary to map sorted nucleotide combinations to IUPAC codes
    mapping = {
        "A": "A",
        "T": "T",
        "C": "C",
        "G": "G",
        "AG": "R",
        "CT": "Y",
        "CG": "S",
        "AT": "W",
        "GT": "K",
        "AC": "M",
        "CGT": "B",
        "AGT": "D",
        "ACT": "H",
        "ACG": "V"
    }

    # Function to sort characters inside brackets and replace with IUPAC code
    def replace_with_iupac(match):
        sorted_seq = "".join(sorted(match.group(1)))
        return mapping.get(sorted_seq, "")

    # Replace all [xy] patterns with their corresponding IUPAC code
    iupac_seq = re.sub(r'\[([A|T|C|G]+)\]', replace_with_iupac, regex)

    # Replace remaining dots with N for any nucleotide
    iupac_seq = iupac_seq.replace(".", "N")
    
    return iupac_seq

def iupac_to_regex(iupac):
    # Dictionary to map sorted nucleotide combinations to IUPAC codes
    mapping = {
        "A": "A",
        "T": "T",
        "C": "C",
        "G": "G",
        "R": "AG",
        "Y": "CT",
        "S": "CG",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
        "N": "."
    }

    # Function to sort characters inside brackets and replace with IUPAC code
    def replace_with_iupac(match):
        regex_seq = mapping.get(match.group(0))
        if len(regex_seq) == 1:
            return regex_seq
        else:
            return "[" + mapping.get(match.group(0)) + "]"

    # Replace all [xy] patterns with their corresponding IUPAC code
    regex_seq = re.sub(r'.', replace_with_iupac, iupac)
    return regex_seq



def regex_to_iupac_n(regex):
    # Dictionary to map sorted nucleotide combinations to IUPAC codes
    mapping = {
        "A": "A",
        "T": "T",
        "C": "C",
        "G": "G",
        "AG": "R",
        "CT": "Y",
        "CG": "S",
        "AT": "W",
        "GT": "K",
        "AC": "M",
        "CGT": "B",
        "AGT": "D",
        "ACT": "H",
        "ACG": "V"
    }

    # Function to sort characters inside brackets and replace with IUPAC code
    def replace_with_iupac(match):
        sorted_seq = "".join(sorted(match.group(1)))
        return mapping.get(sorted_seq, "")

    # Replace all [xy] patterns with their corresponding IUPAC code
    iupac_seq = re.sub(r'\[([A|T|C|G]+)\]', replace_with_iupac, regex)

    # Replace remaining dots with N for any nucleotide
    iupac_seq = iupac_seq.replace(".", "N")

    # Replace consecutive N's with N followed by the count (only if 2 or more)
    def replace_ns(match):
        length = len(match.group(0))
        return r'$\mathregular{N_' + str(length) + '}$' if length > 4 else 'N'

    iupac_seq = re.sub(r'N+', replace_ns, iupac_seq)
    
    return iupac_seq


def reverse_compliment(seq):
    """Return the reverse complement of the sequence."""
    return "".join([COMPLEMENT[base] for base in reversed(seq)])
