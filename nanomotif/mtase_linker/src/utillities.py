#%% Utility functions

#%% Reverse complement
def reverse_complement(seq):
    """
    Generate the reverse complement of a DNA sequence.

    Given a string representing a DNA sequence, this function returns its 
    reverse complement. The reverse complement is formed by first reversing 
    the original sequence and then taking the complement of each nucleotide. 
    This function recognizes standard DNA bases (A, T, G, C), ambiguous 
    nucleotides (R, Y, S, W, K, M, B, D, H, V, N), and special characters 
    (., [, and ]).

    Parameters:
    seq (str): A string representing the DNA sequence.

    Returns:
    str: The reverse complement of the given DNA sequence.

    Raises:
    KeyError: If the sequence contains characters not defined in the 
    complement mapping.

    Examples:
    >>> reverse_complement("ATGC")
    "GCAT"

    >>> reverse_complement("RYSWKM")
    "MKSWYR"

    Note:
    The function assumes that the input sequence is in uppercase. If the input
    contains lowercase characters, they will not be recognized and a KeyError
    will be raised.
    """

    complement = {
        "A": "T", "T": "A", "G": "C", "C": "G", 
        "N": "N", "R": "Y", "Y": "R", "S": "S", 
        "W": "W", "K": "M", "M": "K", "B": "V", 
        "D": "H", "H": "D", "V": "B", ".": ".",
        "[": "]", "]": "[", "N": "N"
    }

    return ''.join(complement[base] for base in reversed(seq))

#%% Function to predict motif type based on RM-type
def motif_type_predictions(df, subtype_column):
    """
    Adds a new column 'motif type' to the DataFrame based on predefined match conditions.

    Parameters:
    - df: MTase table.
    - subtype_column (str): A column which state the RM type of the MTase.

    Returns:
    - pandas.DataFrame: The modified DataFrame with a new 'motif type' column added.

    The function matches each row in the 'subtype_column' against predefined strings
    for four motif types: 'bipartite', 'palindromic', 'non-palindromic', and 'unknown'.
    It then assigns a corresponding label to each row in the new 'motif type' column.
    """

   
    match_string1 = ['Type_I']  
    match_string2 = ['Type_II'] 
    match_string3 = ['Type_III']  
    match_string4 = ['Type_IIG'] 
    output_string1 = 'bipartite'
    output_string2 = 'palindrome'
    output_string3 = 'non-palindrome'
    output_string4 = 'bipartite, non-palindrome'

    def classify_motif(x):
        if x in match_string1:
            return output_string1
        elif x in match_string2:
            return output_string2
        elif x in match_string3:
            return output_string3
        elif x in match_string4:
            return output_string4
        else:
            return 'undefined' 

    df['motif_type'] = df[subtype_column].apply(classify_motif)

    return df

#%% Function to convert HMM-hit string to RM subtype type
def RM_type_converter(df, gene_name_column):
    """
    Adds a new column 'sub_type' to the DataFrame based on predefined match conditions.

    Parameters:
    - df: MTase table.
    - subtype_column (str): A column which state the name of the HMM, which identified the MTase.

    Returns:
    - pandas.DataFrame: The modified DataFrame with a new 'sub_type' column added.

    The function matches each entry in the 'gene name' column against predefined strings
    for the HMM names: 'RM__Type_I_MTases', 'RM_Type_II__Type_II_MTases', 'RM_Type_IIG__Type_IIG', and 'RM_Type_III__Type_III_MTases'.
    It then assigns a corresponding RM type to each entry in the new 'sub_type' column.
    """

 
    match_string1 = 'RM__Type_I_MTases' 
    match_string2 = 'RM_Type_II__Type_II_MTases'
    match_string3 = 'RM_Type_IIG__Type_IIG'
    match_string4 = 'RM_Type_III__Type_III_MTases'
    output_string1 = 'Type_I'
    output_string2 = 'Type_II'
    output_string3 = 'Type_IIG'
    output_string4 = 'Type_III'

    def classify_motif(x):
        if match_string1 in x:
            return output_string1
        elif match_string2 in x:
            return output_string2
        elif match_string3 in x:
            return output_string3
        elif match_string4 in x:
            return output_string4
        else:
            return 'undefined'  

    df['sub_type'] = df[gene_name_column].apply(classify_motif)

    return df

#%% Functions to predict mod type based on HMM hits

def mod_predictions_hmm(df, pfam_hmm_acc_column):
    """
    Creates a new column in the pfam hit table predicting mod type based on pfam family hits.

    Parameters:
    - df (pandas.DataFrame): pfam_hit_df.
    - pfam_hmm_acc_column: the pfam hmm accession number column.

    Returns:
    - Return the df with the new mod type column added.
    """

    a_interpro_acc = ["PF01555", "PF02384", "PF12161", "PF05869", "PF02086", "PF07669", "PF13651"]
    m_interpro_acc = ["PF00145"]

    def classify_mod(x):
        if x.split('.')[0] in a_interpro_acc:
            return "ac"
        elif x.split('.')[0] in m_interpro_acc:
            return "m"
        else:
            return None

    df['mod_type'] = df[pfam_hmm_acc_column].apply(classify_mod)

    return df
# %%
#Function to recode nanomotif bin-motifs modtype values.
def recode_mod_type(value):
    if value == 'm':
        return 'm'
    elif value == 'a' or value == '21839':
        return 'ac'
    else:
        return value