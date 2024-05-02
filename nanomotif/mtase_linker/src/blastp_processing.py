# pyright: reportUndefinedVariable=false
#%% Load modules
import pandas as pd

#%% Importing blastp results
blastp_table_path = snakemake.input[0]
blastp_df = pd.read_csv(blastp_table_path, sep='\t')

#%%
if blastp_df.empty == False:

    #Extract motif information from subject information
    blastp_df[['REbase_id', 'motif_length']] = blastp_df['salltitles'].str.split('   ', expand = True)
    blastp_df[['motif', 'aa_length']] = blastp_df['motif_length'].str.split('  ', expand = True)
    blastp_df.drop(['motif_length', 'aa_length'], axis=1, inplace=True)


    # Filter blastp results based on identity and coverage parameters.
    qseq_groups = blastp_df.groupby('qseqid')
    idx_to_keep = []
    for name, group in qseq_groups:
        qseqs_sign = group.loc[(group['pident'] >= snakemake.params['identity']) & (group['qcovs'] >= snakemake.params['qcovs'])]
        if len(qseqs_sign) > 1:
            idx = qseqs_sign['bitscore'].idxmax()
            idx_to_keep.append(idx)
        elif len(qseqs_sign) == 1:
            idx_to_keep.append(qseqs_sign.index[0])
        else:
            pass

    blastp_df_sign = blastp_df.loc[idx_to_keep]

    # Export filterede blastp results
    blastp_df_sign.to_csv(snakemake.output[0], sep = "\t", index = False)

else:
    blastp_df.to_csv(snakemake.output[0], sep = "\t", index = False)
