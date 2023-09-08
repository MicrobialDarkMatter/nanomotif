

class Assembly():
    def __init__(self, assembly: dict):
        self.assembly = assembly
    
    def __getitem__(self, key):
        return self.assembly[key]
    
    def __repr__(self):
        return f"Assembly with {len(self.assembly)} contigs"



def reverse_complement(seq):
    """
    Reverse complement a sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    return ''.join([complement[base] for base in reversed(seq)])
