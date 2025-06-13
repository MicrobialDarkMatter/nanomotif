# pyright: reportUndefinedVariable=false
#The above line removes warning messages coming from the snakemake.input and .output parameters

#%% Importing packages and utility functions
import pandas as pd
import numpy as np
from utillities import reverse_complement
from utillities import motif_type_predictions
from utillities import RM_type_converter
from utillities import recode_mod_type

#%% Importing modtype table
gene_mod_tsv_path = snakemake.input['mod_table_path']
gene_mod_df = pd.read_csv(gene_mod_tsv_path, sep = '\t', header = 0)


#%% Importing Defensefinder table
DefenseFinder_df_tsv_path = snakemake.input['DefenseFinder_path']
DefenseFinder_df = pd.read_csv(DefenseFinder_df_tsv_path, sep = "\t", header = 0)
DefenseFinder_df = RM_type_converter(DefenseFinder_df, "gene_name") # Converting 'gene name' column to 'sub_type'
DefenseFinder_df = DefenseFinder_df[["hit_id", "replicon", "sub_type", "RM_system", "sys_id"]] 
DefenseFinder_df.columns = ["gene_id", "contig", "sub_type", "RM_system", "sys_id"]


#%% Importing BLASTP table
BLASTP_df_tsv_path = snakemake.input['BLASTP_path']
BLASTP_df = pd.read_csv(BLASTP_df_tsv_path, sep = "\t", header = 0, index_col=False)
BLASTP_df = BLASTP_df[["qseqid", "REbase_id", "motif"]]
BLASTP_df.columns = ['gene_id', 'REbase_ID', 'motif_guess']

DefenseFinder_blastp_df = DefenseFinder_df.merge(BLASTP_df, on = 'gene_id', how = 'left')

#%% Importing nanomotif table
nanomotif_table_path = snakemake.input['nanomotif_table_path']
nanomotif_table = pd.read_csv(nanomotif_table_path, sep = "\t", header = 0)
nanomotif_table['assign_mod_type'] = nanomotif_table['mod_type'].apply(recode_mod_type)

#Set motif acceptance threshold as >=0.5 mean_methylation and remove ambiguous motifs in nanomotif_table
mean_methylation = nanomotif_table['n_mod'] / (nanomotif_table['n_mod'] + nanomotif_table['n_nomod'])
nanomotif_table_mm50 = nanomotif_table[mean_methylation >= snakemake.params['MINIMUM_METHYLATION']]



#%% Importing bin_contig table from mmlong2 
contig_bin_tsv_path = snakemake.input['contig_bin']
contig_bin_df = pd.read_csv(contig_bin_tsv_path , sep = "\t", header = None)
contig_bin_df.columns = ['contig', 'bin name']


#%% Merging dataframes
# Merging bin origin information into DefenseFinder based on contig number.
bin_DF_blastp_df = DefenseFinder_blastp_df.merge(contig_bin_df[['bin name', 'contig']], on = 'contig', how = 'left')
bin_DF_blastp_df_filtered = bin_DF_blastp_df.dropna(subset=['bin name'])

# Merging modtype table into Defensefinder table
MTase_table = bin_DF_blastp_df_filtered.merge(gene_mod_df[['gene_id', 'mod_type']], on = 'gene_id', how = 'left')


#%% adding methylation motif type column based on MTase subtype
MTase_table = motif_type_predictions(MTase_table, 'sub_type')
MTase_table_assigned = MTase_table.copy()


#%% Add "a" mod type known for RM type I, type IIG and type III
for idx, row in MTase_table_assigned.iterrows():
    if row['sub_type'] == 'Type_I' or row['sub_type'] == 'Type_III' or row['sub_type'] == 'Type_IIG':
        row['mod_type'] = "ac"


#%% Preparing dataframes
        
#Extracting all entries with missing 'mod type' prediction
nan_genes = MTase_table_assigned[MTase_table_assigned['mod_type'].isna()]

# Generate entries with pseudo mod types for all entries with missing 'mod type' prediction
extra_entries = []

for index, entry in nan_genes.iterrows():
        # Create a copy of the current nan gene for 'a' and 'm'
        new_entry_ac = entry.copy()
        new_entry_m = entry.copy()

        # Set 'mod_type' to 'ac' and 'm' respectively
        new_entry_ac['mod_type'] = 'ac'
        new_entry_m['mod_type'] = 'm'


        extra_entries.append(new_entry_ac)
        extra_entries.append(new_entry_m)

extra_entries_df = pd.DataFrame(extra_entries)

# Add entries with pseudo mod types to MTase_table_assigned resulting in combined_MTase
combined_MTase_df = pd.concat([MTase_table_assigned, extra_entries_df]) #Indicies stay the same for extra entries and original entry


# Initialize a column to track whether an entry has been assigned
nanomotif_table_mm50 = nanomotif_table_mm50.copy()
nanomotif_table_mm50.loc[:,'linked'] = False
MTase_table_assigned = MTase_table_assigned.copy()
MTase_table_assigned.loc[:,'linked'] = False


#Group by 'bin name' and 'mod type'
Cgrouped_MTase = combined_MTase_df.groupby(['bin name', 'mod_type'], dropna=True)
grouped_nanomotif = nanomotif_table_mm50.groupby(['bin', 'assign_mod_type'])


#%% Assigment rule system

# Step 1: Priority 1 Assignment based on similarity to motif guess
for (bin_name, assign_mod_type), nanomotif_group in grouped_nanomotif:
    if (bin_name, assign_mod_type) in Cgrouped_MTase.groups:
        mtase_group = Cgrouped_MTase.get_group((bin_name, assign_mod_type))
        for idx, row in nanomotif_group.iterrows():
            # Calculate reverse complement of the motif
            rev_comp_motif = reverse_complement(row['motif'])

            # Find matching motif guess in MTase_table and assign
            matching_indices = mtase_group.loc[mtase_group['motif_guess'] == row['motif']].index
            matching_rev_comp_indices = mtase_group.loc[mtase_group['motif_guess'] == rev_comp_motif].index
            

            # Expand matching_indices to include sys_id matches
            if not matching_indices.empty:
                sys_id_value = mtase_group.loc[matching_indices[0], 'sys_id']  # Get sys_id of the first match
                sys_id_matching_indices = mtase_group.loc[mtase_group['sys_id'] == sys_id_value].index
                matching_indices = matching_indices.union(sys_id_matching_indices)

            # Similarly expand for reverse complement matches
            if not matching_rev_comp_indices.empty:
                sys_id_value_rev_comp = mtase_group.loc[matching_rev_comp_indices[0], 'sys_id']
                sys_id_matching_rev_comp_indices = mtase_group.loc[mtase_group['sys_id'] == sys_id_value_rev_comp].index
                matching_rev_comp_indices = matching_rev_comp_indices.union(sys_id_matching_rev_comp_indices)

            # Update detected motif
            MTase_table_assigned.loc[matching_indices, 'detected_motif'] = row['motif']
            MTase_table_assigned.loc[matching_rev_comp_indices, 'detected_motif'] = rev_comp_motif


            # Mark as linked
            MTase_table_assigned.loc[matching_indices, 'linked'] = True
            MTase_table_assigned.loc[matching_rev_comp_indices, 'linked'] = True

            # If there are matching indices, mark as linked in nanomotif_table
            if not matching_indices.empty:
                nanomotif_group.loc[idx, 'linked'] = True
                genes_str = ', '.join(MTase_table_assigned.loc[matching_indices, 'gene_id'].astype(str).unique())
                nanomotif_table_mm50.at[idx, 'candidate_genes'] = genes_str
            if not matching_rev_comp_indices.empty:
                nanomotif_group.loc[idx, 'linked'] = True
                genes_str = ', '.join(MTase_table_assigned.loc[matching_rev_comp_indices, 'gene_id'].astype(str).unique())
                nanomotif_table_mm50.at[idx, 'candidate_genes'] = genes_str

        # After processing each group, update the 'assigned' status back to the main dataframe
        nanomotif_table_mm50.update(nanomotif_group['linked'])

MTase_table_assigned['detected_motif'] = MTase_table_assigned['detected_motif'].replace('nan', np.nan)
#%%
#%% Re-group to exclude assigned entries
grouped_nanomotif = nanomotif_table_mm50[nanomotif_table_mm50['linked'] != True].groupby(['bin', 'assign_mod_type'])
idx_MTase_table_assigned = MTase_table_assigned[MTase_table_assigned['linked'] != True].index #Extract indicies of non assigned entries 
#%%
combined_MTase_df = combined_MTase_df.loc[idx_MTase_table_assigned] #FIlter combined_MTase based on extracted indicies
Cgrouped_MTase = combined_MTase_df.groupby(['bin name', 'mod_type'], dropna=True)
#%%

# Step 2: Priority 2 Assignment: based single subtypes in group
for (bin_name, assign_mod_type), nanomotif_group in grouped_nanomotif:
    if (bin_name, assign_mod_type) in Cgrouped_MTase.groups:
        mtase_group = Cgrouped_MTase.get_group((bin_name, assign_mod_type))

        # Count unique motif types in nanomotif and MTase table
        for idx, row in nanomotif_group.iterrows():
            motif_type = row['motif_type']
            candidate_genes = mtase_group[mtase_group['motif_type'].apply(lambda x: motif_type in x.split(', '))][['gene_id', 'sys_id']]
            
            # Count unique sys_ids
            unique_sys_ids = candidate_genes['sys_id'].unique()
            if len(unique_sys_ids) == 1 and not pd.isna(unique_sys_ids).any():
                count_in_mtase = 1
            elif len(candidate_genes) == 1 and pd.isna(unique_sys_ids).any():
                count_in_mtase = 1
            else:
                count_in_mtase = len(candidate_genes)

            count_in_nanomotif = nanomotif_group['motif_type'].eq(motif_type).sum()

            # Check if the motif_type is unique in both groups
            if count_in_mtase == 1 and count_in_nanomotif == 1:
                
                #Assigning single candidate genes
                genes_list = candidate_genes['gene_id'].tolist()
                genes_str = ', '.join(genes_list)
                nanomotif_table_mm50.loc[idx, 'candidate_genes'] = genes_str

                if motif_type in ['non-palindrome', 'bipartite']:
                    # Explicitly handle 'bipartite, non-palindrome' if motif_type is 'non-palindrome' or 'bipartite'
                    mtase_index = mtase_group[(mtase_group['motif_type'].eq(motif_type)) | (mtase_group['motif_type'] == 'bipartite, non-palindrome')].index
                else:
                    # Find the corresponding index in mtase_group
                    mtase_index = mtase_group[mtase_group['motif_type'].eq(motif_type)].index
                
                # Assign the motif
                MTase_table_assigned.loc[mtase_index, 'detected_motif'] = row['motif']
                MTase_table_assigned.loc[mtase_index, 'linked'] = True
                nanomotif_table_mm50.loc[idx, 'linked'] = True

MTase_table_assigned['detected_motif'] = MTase_table_assigned['detected_motif'].replace('nan', np.nan)

#%% Re-group to exclude assigned entries
grouped_nanomotif = nanomotif_table_mm50[nanomotif_table_mm50['linked'] != True].groupby(['bin', 'assign_mod_type'])
idx_MTase_table_assigned = MTase_table_assigned[MTase_table_assigned['linked'] != True].index #Extract indicies of non assigned entries
combined_MTase_df = combined_MTase_df.loc[idx_MTase_table_assigned] #Filter combined_MTase based on extracted indicies
Cgrouped_MTase = combined_MTase_df.groupby(['bin name', 'mod_type'], dropna=True)

#%% Step 3: Assign candidate genes
for (bin_name, assign_mod_type), nanomotif_group in grouped_nanomotif:
    if (bin_name, assign_mod_type) in Cgrouped_MTase.groups:
        mtase_group = Cgrouped_MTase.get_group((bin_name, assign_mod_type))

        # Count unique motif types in nanomotif and MTase table
        for idx, row in nanomotif_group.iterrows():
            motif_type = row['motif_type']            
            candidate_genes = mtase_group[mtase_group['motif_type'].apply(lambda x: motif_type in x.split(', '))][['gene_id', 'sys_id']]

            # Count unique sys_ids
            unique_sys_ids = candidate_genes['sys_id'].unique()
            if len(unique_sys_ids) == 1 and not pd.isna(unique_sys_ids).any():
                count_in_mtase = 1
            elif len(candidate_genes) == 1 and pd.isna(unique_sys_ids).any():
                count_in_mtase = 1
            else:
                count_in_mtase = len(candidate_genes)

            count_in_nanomotif = nanomotif_group['motif_type'].eq(motif_type).sum()
            
            # Check if the motif_type is unique in both groups
            if count_in_mtase == 1 and count_in_nanomotif == 1:
                
                #Assigning single candidate genes
                genes_list = candidate_genes['gene_id'].tolist()
                genes_str = ', '.join(genes_list)
                nanomotif_table_mm50.loc[idx, 'candidate_genes'] = genes_str

                if motif_type in ['non-palindrome', 'bipartite']:
                    # Explicitly handle 'bipartite, non-palindrome' if motif_type is 'non-palindrome' or 'bipartite'
                    mtase_index = mtase_group[(mtase_group['motif_type'].eq(motif_type)) | (mtase_group['motif_type'] == 'bipartite, non-palindrome')].index
                else:
                    # Find the corresponding index in mtase_group
                    mtase_index = mtase_group[mtase_group['motif_type'].eq(motif_type)].index
                
                # Assign the motif
                MTase_table_assigned.loc[mtase_index, 'detected_motif'] = row['motif']
                MTase_table_assigned.loc[mtase_index, 'linked'] = True
                nanomotif_table_mm50.loc[idx, 'linked'] = True

            #Assigning all candidate genes
            if count_in_mtase >= 1:
                genes_list = candidate_genes['gene_id'].tolist()
                genes_str = ', '.join(genes_list)

                nanomotif_table_mm50.loc[idx, 'candidate_genes'] = genes_str


# %%
MTase_table_assigned = MTase_table_assigned[['bin name', 'gene_id', 'contig', 'mod_type', 'sub_type', 'RM_system', 'sys_id', 'motif_type', 'REbase_ID', 'motif_guess', 'linked', 'detected_motif']]
MTase_table_assigned.columns = ['bin', 'gene_id', 'contig', 'mod_type_pred', 'sub_type_pred', 'RM_system', 'DF_system_ID', 'motif_type_pred', 'REbase_ID', 'motif_pred', 'linked', 'detected_motif']
MTase_table_assigned_cl = MTase_table_assigned.dropna(subset=['bin'])

nanomotif_table_mm50 = nanomotif_table_mm50[['bin', 'mod_type', 'motif', 'mod_position', 'n_mod', 'n_nomod', 'motif_type', 'motif_complement', 'mod_position_complement', 'n_mod_complement', 'n_nomod_complement', 'linked', 'candidate_genes']]
#%%
MTase_table_assigned_cl.to_csv(snakemake.output['MTase_assignment_table'] , sep='\t', index=False)
nanomotif_table_mm50.to_csv(snakemake.output['nanomotif_assignment_table'] , sep='\t', index=False)

