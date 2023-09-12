from scipy.stats import entropy
import numpy as np
import warnings
import re
import nanomotif.utils
import matplotlib.pyplot as plt
from matplotlib import cm

class Assembly():
    def __init__(self, assembly: dict):
        self.assembly = {name:DNAsequence(seq) for name, seq in assembly.items()}
    
    def __getitem__(self, key):
        return self.assembly[key]
    
    def __repr__(self):
        return f"Assembly with {len(self.assembly)} contigs"



def reverse_complement(seq):
    """
    Reverse complement a sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "[":"]", "]":"[", ".": "."} 
    return ''.join([complement[base] for base in reversed(seq)])



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
    
    def sample_subsequence(self, length):
        """Randomly sample a subsequence of a given length"""
        if length > len(self):
            raise ValueError(f"Cannot sample a subsequence of length {length} from a sequence of length {len(self)}")
        start_index = random.randint(0, len(self) - length)
        return self.sequence[start_index:start_index + length]

    def sample_n_subsequences(self, length: int, n: int):
        """Randomly sample n sequences of a given length"""
        return [self.sample_subsequence(length=length) for _ in range(n)]

    def find_subsequence(self, subsequence: str):
        """
        Find all occourances of subseqeunce in DNAsequence
        """
        index = [match.start() for match in re.finditer(subsequence, self.sequence)]
        return np.array(index)

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
            warnings.warn("Not enough occourances of " + subsequence + " found (n: " + str(n) + ")")
        else:
            sample_index = np.random.choice(subsequence_index, n)
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
        assert all_lengths_equal(sequences), "All sequences must be of equal length"

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
        return EqualLengthDNASet([seq for seq in self.sequences if not bool(re.match(pattern, seq.sequence))])
    def plot_pssm(self, ax=None, center_x_axis=True):
        """
        Method to plot the positional conservation of DNA sequences in the object.

        Parameters:
        ax (matplotlib.axes.Axes, optional): An existing Axes object to draw the plot onto, 
                                            default is None in which case a new figure will be created.
        Returns:
        None

        This method generates a plot showing the conservation of sequences in the object.
        It plots the positional frequencies of each nucleotide.
        """
        
        # If no Axes object provided, create a new figure and axes
        if ax is None:
            _, ax = plt.subplots()

        # Plot positional frequencies for each nucleotide
        int_base = self.sequences[0].int_base
        for i in range(1,5):
            ax.plot(self.pssm()[i-1], label=f"{int_base[i]}", linewidth=2, c=cm.cubehelix(i/5))

        # Configure x-axis labels and ticks
        if center_x_axis:
            n_labs = min(4, (self.sequence_length-1)//2) 
            window = (self.sequence_length-1)//2
            x_axis_step = max(window//n_labs, 1)
            labs_positive = np.arange(0, x_axis_step * n_labs + 1, x_axis_step)
            tick_label = np.concatenate((-labs_positive, labs_positive))
            tick_position = tick_label + window
            ax.set_xticks(tick_position)
            ax.set_xticklabels(tick_label)
            ax.set_xlabel("Relative Position")
        else:
            ax.set_xticks(np.arange(0, self.sequence_length, 1))
            ax.set_xlabel("Position")
        # Set title and labels
        ax.set_title(f"Number of sequences: {len(self.sequences)}")
        ax.set_ylabel("Frequency")

        return ax

def all_lengths_equal(iterator):
    """
    Checks whether the lengths of all elements in an iterable are equal.

    The function will return True even if the iterable is empty. It requires that the elements of the iterable
    also be iterable, such as strings, lists, and tuples.

    Args:
        iterator (Iterable): An iterable object containing other iterable elements.

    Returns:
        bool: True if all lengths are equal or if iterable is empty, False otherwise.

    Examples:
    >>> all_lengths_equal(['abc', 'def', 'ghi'])
    True
    >>> all_lengths_equal([[1, 2, 3], [4, 5, 6]])
    True
    >>> all_lengths_equal([])
    True
    >>> all_lengths_equal(['abc', 'de'])
    False
    """
    iterator = iter(iterator)
    try:
        first = len(next(iterator))
    except StopIteration:
        return True
    return all(first == len(x) for x in iterator)