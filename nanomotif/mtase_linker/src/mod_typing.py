# pyright: reportUndefinedVariable=false
#%% load libraries
import pandas as pd
import os
from utillities import mod_predictions_hmm


#%% Function to decide the mod type of the MTase protein sequences.
def select_and_filter_pfam_entries(pfam_hit_mod_df):
    """
    Based on a user defined hierachy, this function decides whether the protein belong to modtype "a" or "m" or None

    Parameters:
    pfam_hit_mod_df: pfam_hit_table with mod type prediction for each entry
    """

    # Define hierarchies
    hierarchy_a = ["PF02086", "PF05869", "PF02384", "PF01555", "PF1216", "PF07669", "PF13651"]
    hierarchy_m = ["PF0014"]

    # Function to select modtype based on hierarchy
    def select_row_based_on_hierarchy(name, group, hierarchy_a, hierarchy_m):
        group_a = group[group['mod_type'] == 'a'].copy()
        group_m = group[group['mod_type'] == 'm'].copy()

        if not group_a.empty and not group_m.empty:
            return pd.DataFrame({'gene_id': [name], 'mod_type': ['ambiguous']})
        elif not group_a.empty:
            group = group_a
            hierarchy = hierarchy_a
        elif not group_m.empty:
            group = group_m
            hierarchy = hierarchy_m
        else:
            return group.head(1)

        group['rank'] = group['HMM_acc'].apply(lambda x: hierarchy.index(x.split('.')[0]) if x.split('.')[0] in hierarchy else len(hierarchy))
        return group.sort_values('rank').head(1)

    # Iterate over the groups in the dataframe
    return pd.concat([select_row_based_on_hierarchy(name, group, hierarchy_a, hierarchy_m) for name, group in pfam_hit_mod_df.groupby('gene_id')])

#%% Defining pfam hit table ath
pfam_hit_tsv_path = snakemake.input[0]

if os.path.getsize(pfam_hit_tsv_path) > 0: 
    
    # Importing pfam hit table
    pfam_hit_df = pd.read_csv(pfam_hit_tsv_path, sep = '\t', header = None)
    pfam_hit_df = pfam_hit_df.loc[:, 0:3]
    pfam_hit_df.columns = ["gene_id", "acc", "HMM_name", "HMM_acc"]

    #Adding modtype prediction to pfam hits
    pfam_hit_mod_df = mod_predictions_hmm(pfam_hit_df, 'HMM_acc')

    # Running above function to extract gene name and modtype table
    protein_acc_mod_df = select_and_filter_pfam_entries(pfam_hit_mod_df)[["gene_id", "mod_type"]]
    protein_acc_mod_df.to_csv(snakemake.output[0] , sep='\t', index=False)

else:
    open(snakemake.output[0], 'w').close()
