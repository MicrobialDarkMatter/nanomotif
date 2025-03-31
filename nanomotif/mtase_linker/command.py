import subprocess
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

    try:
        identity = float(args.identity)
    except ValueError:
        sys.stderr.write("Error: identity must be a number.\n")
        sys.exit(1)

    if identity < 0 or identity > 100:
        sys.stderr.write("Error: identity must be between 0 and 100.\n")
        sys.exit(1)

    try:
        QCOVS = float(args.qcovs)
    except ValueError:
        sys.stderr.write("Error: query coverage must be a number.\n")
        sys.exit(1)

    if QCOVS < 0 or QCOVS > 100:
        sys.stderr.write("Error: query coverage must be between 0 and 100.\n")
        sys.exit(1)

    try:
        MINIMUM_METHYLATION = float(args.minimum_motif_methylation)
    except ValueError:
        sys.stderr.write("Error: minimum motif methylation must be a number.\n")
        sys.exit(1)

    if MINIMUM_METHYLATION < 0 or MINIMUM_METHYLATION > 1:
        sys.stderr.write("Error: minimum motif methylation must be between 0 and 1.\n")
        sys.exit(1)

    command = [
        "snakemake",
        "--snakefile", snakefile,
        "--cores", str(args.threads),
        "--config",
        f"THREADS={args.threads}",
        f"ASSEMBLY={assembly_path}",
        f"CONTIG_BIN={contig_bin}",
        f"OUTPUTDIRECTORY={args.out}",
        f"DEPENDENCY_PATH={dependency_dir}",
        f"IDENTITY={args.identity}",
        f"QCOVS={args.qcovs}",
        f"NANOMOTIF={args.bin_motifs}",
        f"MINIMUM_METHYLATION={args.minimum_motif_methylation}",
        "--use-conda",
        "--conda-prefix", os.path.join(dependency_dir, "ML_envs")
        
        ]
    
    if args.forceall:
        command.append("--forceall")
    if args.dryrun:
        command.append("--dryrun")
        
    try:
        result = subprocess.run(command, cwd=cwd, check=True)
        if result.returncode == 0:
            print("MTase-Linker done!")
    except subprocess.CalledProcessError as e:
        print("MTase-linker failed with error:", e)
        sys.exit(1)
