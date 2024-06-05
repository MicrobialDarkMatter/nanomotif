import snakemake
import os
import sys
import glob

def run_MTase_linker(args):
    
    thisdir = os.path.abspath(os.path.dirname(__file__))
    cwd = os.getcwd()

    snakefile_path = os.path.join(thisdir,"MTase_linker.smk")
    if os.path.exists(snakefile_path):
        snakefile = snakefile_path
        print(snakefile_path)
    else:
        msg = 'Error: cannot find Snakefile at the following location:\n'
        msg += '{}\n'.format(snakefile_path)
        sys.stderr.write(msg)
        sys.exit(1)

    # Bins directory
    a1 = os.path.join(cwd, args.assembly)
    a2 = os.path.join(args.assembly)
    if os.path.exists(a1):
        assembly_path = a1
    elif os.path.exists(a2):
        assembly_path = a2
    else:
        msg = 'Error: cannot find assembly files in paths {} or {} '.format(a1,a2)
        sys.stderr.write(msg)
        sys.exit(1)


    # contig_bin.tsv
    b1 = os.path.join(cwd, args.contig_bin)
    b2 = os.path.join(args.contig_bin)
    if os.path.exists(b1) and not os.path.isdir(b1) and os.path.splitext(b1)[1] == '.tsv':
        contig_bin = b1
    elif os.path.exists(b2) and not os.path.isdir(b2) and os.path.splitext(b2)[1] == '.tsv':
        contig_bin = b2
    elif os.path.exists(b1) and not os.path.isdir(b1):
        msg = 'Error: contig_bin file {} does not have .tsv. extension'.format(b1)
        sys.stderr.write(msg)
        sys.exit(1)
    elif os.path.exists(b2) and not os.path.isdir(b2):
        msg = 'Error: contig_bin file {} does not have .tsv. extension'.format(b2)
        sys.stderr.write(msg)
        sys.exit(1)
    else:
        msg = 'Error: cannot find bins directory in paths {} or {} '.format(b1,b2)
        sys.stderr.write(msg)
        sys.exit(1)


    c1 = os.path.join(cwd, args.dependency_dir)
    c2 = os.path.join(args.dependency_dir)
    if os.path.exists(c1) and os.path.isdir(c1):
        dependency_dir = c1
    elif os.path.exists(c2) and os.path.isdir(c2):
        dependency_dir = c2
    else:
        msg = 'Dependency dir not found at {} or {}, run MTase-linker install'.format(c1, c2)
        sys.stderr.write(msg)
        sys.exit(1)
    
    if os.path.exists(os.path.join(c1, "df_models", "defense-finder_update.done")):
        print("Defensefinder models found!")
    elif os.path.exists(os.path.join(c2, "df_models", "defense-finder_update.done")):
            print("Defensefinder models found!")
    else:
        msg = 'Defensefinder models not found at {} or {}, run MTase-linker install'.format(c1, c2)
        sys.stderr.write(msg)
        sys.exit(1)

    if os.path.exists(os.path.join(c1, "ML_envs")):
        print("Environments found!")
    elif os.path.exists(os.path.join(c2, "ML_envs")):
            print("Environments found!")
    else:
        msg = 'Environments not found at {} or {}, run MTase-linker install'.format(c1, c2)
        sys.stderr.write(msg)
        sys.exit(1)

    if os.path.exists(os.path.join(c1, "REbase_RecSeq_all")):
        print("REbase database found!")
    elif os.path.exists(os.path.join(c2, "REbase_RecSeq_all")):
            print("REbase database found!")
    else:
        msg = 'REbase database not found at {} or {}, run MTase-linker install'.format(c1, c2)
        sys.stderr.write(msg)
        sys.exit(1)

    if os.path.exists(os.path.join(c1, "PFAM_MTase_profiles.hmm")):
        print("PFAM models found!")
    elif os.path.exists(os.path.join(c2, "PFAM_MTase_profiles.hmm")):
            print("PFAM models found!")
    else:
        msg = 'PFAM models not found at {} or {}, run MTase-linker install'.format(c1, c2)
        sys.stderr.write(msg)
        sys.exit(1)


    workflow = None
    workflow = {"THREADS": args.threads,
                "ASSEMBLY": assembly_path,
                "CONTIG_BIN": contig_bin,
                "OUTPUTDIRECTORY": args.out,
                "DEPENDENCY_PATH": dependency_dir,
                "IDENTITY": args.identity,
                "QCOVS": args.qcovs,
                "NANOMOTIF": args.bin_motifs}


    status = snakemake.snakemake(snakefile, 
                                config=workflow,
                                targets = ["all"], 
                                use_conda = True, 
                                forceall=args.forceall, 
                                cores = args.threads,
                                dryrun=args.dryrun, 
                                workdir = cwd,
                                conda_prefix = os.path.join(dependency_dir, "ML_envs"))

    if status:
        print("MTase-Linker done!")
    else:
        print("MTase-linker failed.")

