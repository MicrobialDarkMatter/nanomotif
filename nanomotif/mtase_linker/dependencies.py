import os
import subprocess
import sys
import shutil
import snakemake

def snakemake_create_environments(args):

    thisdir = os.path.abspath(os.path.dirname(__file__))
    cwd = os.getcwd()

    snakefile_path = os.path.join(thisdir, "setup.smk")
    if os.path.exists(snakefile_path):
        snakefile = snakefile_path
        print(snakefile_path)
    else:
        msg = 'Error: cannot find Snakefile at the following location:\n'
        msg += '{}\n'.format(snakefile_path)
        sys.stderr.write(msg)
        sys.exit(-1)

    dependency_dir = os.path.join(args.dependency_dir, "ML_dependencies")


    command = [
        "snakemake",
        "--snakefile", snakefile,
        "--config", f"DEPENDENCYDIR={dependency_dir}",
        "--use-conda",
        "--conda-prefix", os.path.join(dependency_dir, "ML_envs"),
        "--conda-create-envs-only",
        "--cores", "1",
        "all"
        ]
        
    try:
        result = subprocess.run(command, cwd=cwd, check=True)
        if result.returncode == 0:
            print("Conda environments were successfully created.")
    except subprocess.CalledProcessError as e:
        print("Creation of conda environments failed.", e)
        sys.exit(1)


def get_models(args):

    thisdir = os.path.abspath(os.path.dirname(__file__))
    cwd = os.getcwd()

    snakefile_path = os.path.join(thisdir, "setup.smk")
    if os.path.exists(snakefile_path):
        snakefile = snakefile_path
        print(snakefile_path)
    else:
        msg = 'Error: cannot find Snakefile at the following location:\n'
        msg += '{}\n'.format(snakefile_path)
        sys.stderr.write(msg)
        sys.exit(-1)

    dependency_dir = os.path.join(args.dependency_dir, "ML_dependencies")


    command = [
        "snakemake",
        "--snakefile", snakefile,
        "--config", f"DEPENDENCYDIR={dependency_dir}",
        "--use-conda",
        "--conda-prefix", os.path.join(dependency_dir, "ML_envs"),
        "--cores", "1",
        "make_REbase_db"
        ]
        
    try:
        result = subprocess.run(command, cwd=cwd, check=True)
        if result.returncode == 0:
            print("pfam models and REbase database were successfully retreived.")
    except subprocess.CalledProcessError as e:
        print("pfam models and REbase database could not be retreived.", e)
        sys.exit(1)



def defensefinder_update(args):

    thisdir = os.path.abspath(os.path.dirname(__file__))
    cwd = os.getcwd()

    snakefile_path = os.path.join(thisdir, "setup.smk")
    if os.path.exists(snakefile_path):
        snakefile = snakefile_path
        print(snakefile_path)
    else:
        msg = 'Error: cannot find Snakefile at the following location:\n'
        msg += '{}\n'.format(snakefile_path)
        sys.stderr.write(msg)
        sys.exit(-1)

    dependency_dir = os.path.join(args.dependency_dir, "ML_dependencies")


    command = [
        "snakemake",
        "--snakefile", snakefile,
        "--config", f"DEPENDENCYDIR={dependency_dir}",
        "--use-conda",
        "--conda-prefix", os.path.join(dependency_dir, "ML_envs"),
        "--cores", "1",
        "defense_finder_update"
        ]
        
    try:
        result = subprocess.run(command, cwd=cwd, check=True)
        if result.returncode == 0:
            print("Defensefinder models updated succesfully")
    except subprocess.CalledProcessError as e:
        print("Update of Defensefinder models failed.", e)
        sys.exit(1)


def check_installation_MTase_linker(args):

    thisdir = os.path.abspath(os.path.dirname(__file__))
    cwd = os.getcwd()

    dependency_dir = os.path.join(args.dependency_dir, "ML_dependencies")

    command = [
        "nanomotif", 
        "MTase-linker", 
        "run", 
        "--threads", "1", 
        "--assembly", os.path.join(thisdir, "..", "datasets", "e_coli_assembly.polished.fasta"), 
        "--contig_bin", os.path.join(thisdir, "..", "datasets", "e_coli_contig-bin.tsv"), 
        "--bin_motifs", os.path.join(thisdir, "..", "datasets", "e_coli_bin-motifs.tsv"), 
        "--out", os.path.join(dependency_dir, "check_installation"),
        "-d", dependency_dir, 
        "--forceall", "True"
        ]

    print("Testing installation of MTase-linker.")

    try:
        result = subprocess.run(command, cwd=cwd, check=True)
        if result.returncode == 0:
            print("Installation of MTase-linker was successful.")
    except subprocess.CalledProcessError as e:
        print("Installation of MTase-linker failed.", e)
        sys.exit(1)