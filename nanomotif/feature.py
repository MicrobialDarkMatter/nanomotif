


class Pileup():
    """
    Class for loading pileup files
    """
    def __init__(self, pileup):
        self.pileup = pileup

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
    
