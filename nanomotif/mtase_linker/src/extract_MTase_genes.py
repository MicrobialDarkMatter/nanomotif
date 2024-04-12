# pyright: reportUndefinedVariable=false
#%% Load modules
import pandas as pd
from Bio import SeqIO

#%% Importing defenseFinder hmmer and systems output
df_hmmer_path = snakemake.input['hmmer_file']
df_hmmer = pd.read_csv(df_hmmer_path, sep='\t')

df_systems_path = snakemake.input["system_file"]
df_systems = pd.read_csv(df_systems_path, sep='\t')

#%% Filter defensefinder data for RM systems only.
df_hmmer_RM = df_hmmer[df_hmmer['gene_name'].str.contains('_Type_', na=False)]

# Group by gene_id (hit id), then retrieve the entry with the maximum 'hit_score' in each group
idx = df_hmmer_RM.groupby('hit_id')['hit_score'].idxmax()
df_hmmer_RM_filtered = df_hmmer_RM.loc[idx]

#%% Filter defensefinder data for MTase systems only
df_hmmer_MTase = df_hmmer_RM_filtered[(df_hmmer_RM_filtered['gene_name'].str.contains("MTases")) | (df_hmmer_RM_filtered['gene_name'].str.contains("Type_IIG"))].copy()

#%% Determining whether the MTase is part of an RM system
df_systems_RM = df_systems[df_systems['type'] == "RM"]

# Defining a function to check whether a 'qseqid' is in any 'protein_in_syst' row
def check_RM_system(hit_id, protein_in_syst_column):
    for proteins in protein_in_syst_column:
        protein_list = proteins.split(',')
        if hit_id in protein_list:
            return True
    return False

# Checking whether MTase is part of any system
df_hmmer_MTase.loc[:, 'system'] = df_hmmer_MTase["hit_id"].apply(lambda x: check_RM_system(x, df_systems['protein_in_syst']))

# Checking wether MTase is part of an RM-system
df_hmmer_MTase.loc[:, 'RM_system'] = df_hmmer_MTase["hit_id"].apply(lambda x: check_RM_system(x, df_systems_RM['protein_in_syst']))

# Removing MTase hits which are part of another defense system 
df_hmmer_MTase_filtered = df_hmmer_MTase[~((df_hmmer_MTase['system'] == True) & (df_hmmer_MTase['RM_system'] == False))]
df_hmmer_MTase_filtered.drop('system', axis = 1, inplace = True) 

#%% Generating processed defensefinder.tsv
df_hmmer_MTase_filtered.to_csv(snakemake.output['DF_MTase_table'], sep = "\t", index = False)


#%% Generating MTase sequence files
MTase_sequence_file = snakemake.output['DF_MTase_AA']
prodigal_sequence_file = snakemake.input['prodigal_AAs']

# Load identifiers from the DataFrame
identifiers = df_hmmer_MTase_filtered['hit_id'].tolist()

with open(MTase_sequence_file, 'w') as output_file:
    for record in SeqIO.parse(prodigal_sequence_file, 'fasta'):
        # Check if the identifier is in the list
        if record.id.split(' ')[0] in identifiers:
            # Write the sequence to the new FASTA file
            SeqIO.write(record, output_file, 'fasta')
