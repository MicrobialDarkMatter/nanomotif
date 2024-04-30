import os
import sys
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
    workflow = None
    workflow = {"DEPENDENCYDIR": dependency_dir}
    
    status = snakemake.snakemake(snakefile,
                                config = workflow,
                                targets = ["all"],
                                conda_prefix = os.path.join(args.dependency_dir, "ML_dependencies", "ML_envs"),
                                use_conda = True,
                                conda_create_envs_only = True,
                                cores = 1, workdir = cwd)

     # Check if the workflow executed successfully
    if status:
        print("Conda environments were successfully created.")
    else:
        print("Creation of conda environments failed.")

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
    workflow = None
    workflow = {"DEPENDENCYDIR": dependency_dir}
    
    status = snakemake.snakemake(snakefile,
                                config = workflow,
                                targets = ["make_REbase_db"],
                                conda_prefix = os.path.join(args.dependency_dir, "ML_dependencies", "ML_envs"),
                                use_conda = True,
                                cores = 1, workdir = cwd)

     # Check if the workflow executed successfully
    if status:
        print("pfam models and REbase database were successfully retreived.")
    else:
        print("pfam models and REbase database could not be retreived.")


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
    workflow = None
    workflow = {"DEPENDENCYDIR": dependency_dir}


    status = snakemake.snakemake(snakefile, config=workflow,
                                targets = ["defense_finder_update"], 
                                use_conda = True,
                                conda_prefix = os.path.join(args.dependency_dir, "ML_dependencies", 'ML_envs'),
                                cores = 1, workdir = cwd)

     # Check if the workflow executed successfully
    if status:
        print("Defensefinder models updated succesfully")
    else:
        print("Update of Defensefinder models failed.")