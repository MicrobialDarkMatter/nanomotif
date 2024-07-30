# pyright: reportUndefinedVariable=false
#%% load libraries
import pandas as pd
import os
from utillities import mod_predictions_hmm

#%% Defining pfam hit table ath
pfam_hit_tsv_path = snakemake.input[0]

if os.path.getsize(pfam_hit_tsv_path) > 0: 
    
    # Importing pfam hit table
    pfam_hit_df = pd.read_csv(pfam_hit_tsv_path, sep = '\t', header = None)
    pfam_hit_df = pfam_hit_df.loc[:, 0:3]
    pfam_hit_df.columns = ["gene_id", "acc", "HMM_name", "HMM_acc"]
#%%
    #Adding modtype prediction to pfam hits
    pfam_hit_mod_df = mod_predictions_hmm(pfam_hit_df, 'HMM_acc')

    def process_mod_type(group):

        mod_types_concatenated = ', '.join(group['mod_type'])
    
        mod_types_condensed = ', '.join(sorted(set(mod_types_concatenated.split(', '))))
    
    # Check if both 'm' and 'a' are present, if so, choose 'm' only
        if 'm' in mod_types_condensed and 'ac' in mod_types_condensed:
            mod_types_condensed = mod_types_condensed.replace('ac, ', '').replace(', ac', '').replace('ac', '')
    

        return mod_types_condensed

    protein_acc_mod_df = pfam_hit_mod_df.groupby('gene_id').apply(process_mod_type).reset_index()
    protein_acc_mod_df.columns = ['gene_id', 'mod_type']
#%%
    protein_acc_mod_df.to_csv(snakemake.output[0] , sep='\t', index=False)

else:
    open(snakemake.output[0], 'w').close()

