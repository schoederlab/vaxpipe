# vaxpipe

**vaxpipe** is a Snakemake-based pipeline designed to automate structural modeling of symmetric protein complexes using the [Rosetta](https://www.rosettacommons.org/) molecular modeling suite. 

It streamlines the process of the prediction of possible amino acid mutations using deep learning based methods, such as ProteinMPNN and ESM and database driven approaches, such as Rosetta FastDesign, with a focus on shape complementarity of protein interfaces.

---

## Purpose

The goal of `vaxpipe` is to support **vaccine design** workflows by automating the structural modeling of symmetric protein complexes using the Rosetta molecular modeling suite. It focuses on streamlining the process of generating symmetry definition files, performing symmetric relaxation, predicting amino acid mutations using deep learning (ProteinMPNN, ESM), and designing novel sequences with an emphasis on shape complementarity at protein interfaces. This pipeline ensures reproducibility and modularity through Snakemake.

- Generating symmetry definitions using Rosetta's `make_symmdef_file.pl` to speed up the calculations
- predicting amino acid probabilities using deep learning based methods (P-MPNN, ESM) followed by designing novel sequences which are further evaluated using Rosetta
- Ensuring reproducibility and modularity using [Snakemake](https://snakemake.readthedocs.io/). Snakemake is a workflow management system that enables the creation of reproducible and scalable data analysis pipelines. It allows to define rules for how data is processed, automatically handling dependencies and parallelization.

## Operational Guidelines

- **Input PDB Structure**: The pipeline expects an input PDB file (e.g., `test/3ft7.pdb`) which represents a symmetric protein complex.
- **Rosetta Installation**: A functional Rosetta installation compiled with PyTorch and TensorFlow libraries is required for the deep learning-based predictions. For cluster usage, a singularity image can be downloaded from Docker Hub (see details below). Ensure that the `ROSETTA` path is correctly specified in the `config.yaml`.
- **Snakemake**: The pipeline is orchestrated using Snakemake. While users should ideally be familiar with Snakemake, this documentation aims to provide sufficient guidance for basic local and HPC execution. Detailed instructions are provided within the relevant sections.
- **Configuration**: The `config.yaml` file must be updated with appropriate paths to Rosetta binaries, the name of the pdb file, and input/output directories. This file is critical for the pipeline's correct operation.
- **Output Files**: Intermediate and final output files will be generated in a structured directory format, as defined in the `snakefile`. If the pipeline run successfully, several .png files will be present in the output directory, which represent both frequency and energy evaluation of tested mutations.
- **Resource Management**: For HPC execution, users are responsible for configuring SLURM or other cluster submission parameters within the `snakefile-hpc` to match their system's requirements and resource availability.

## Local Execution Requirements

- `Python ≥ 3.7, biopython, tqdm, matplotlib, pandas`
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- Rosetta compiled with pytorch and tensorflow libraries. A detailed information on how to compile Rosetta with pytorch and tensorflow support can be found [here](https://docs.rosettacommons.org/docs/latest/build_documentation/Building-Rosetta-with-TensorFlow-and-Torch)
- ESM model (will be downloaded automatically)

## HPC Requirements
- `Python ≥ 3.7, biopython, tqdm, matplotlib, pandas`
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- **singularity**
- We provide a snakemake file that relies on a Rosetta docker/singularity image. The image is available on Docker Hub and can be pulled using singularity: `singularity pull docker://rosettacommons/rosetta:ml-387` (used for cluster execution).
- ESM model [download](https://git.iwe-lab.de/moritzertelt/ML_graphs/-/tree/main/tensorflow_graphs/ESM/esm2_t33_650M_UR50D). Currently, the pipeline just accepts this ESM model. **The model needs to be downloaded and copied into the repository path.**

## Repository Structure

```
vaxpipe/
├── snakefile # the main workflow
├── snakefile-hpc # the main workflow optimized for the HPC cluster
├── config.yaml # containes path information
├── input_files/ # containes the necessary input .xml files to run the rosetta jobs
└── test/ # containes a test input pdb for validation/benchmarking

```
---

## Installation

generate snakemake conda environment

```
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
pip install biopython tqdm matplotlib
```

**change the `config.yaml` file with the corresponding paths**

---

## Understanding and Executing the Pipeline with Snakemake

`vaxpipe` is built with Snakemake, a powerful workflow management system that helps you create reproducible and scalable data analysis pipelines. Here's a brief guide to get you started:

### Core Concepts:
- **Rules**: Define how output files are generated from input files. Each rule specifies a command to be run.
- **Wildcards**: Allow rules to be applied to many files (e.g., `data/{sample}.txt` where `{sample}` is a wildcard).
- **DAG (Directed Acyclic Graph)**: Snakemake automatically builds a DAG of jobs, ensuring that all dependencies are met before a rule is executed.
- **Execution**: Snakemake determines which rules need to be run based on the desired output files and the current state of your input files.

### Essential Commands:

1.  **Dry Run**: Always start with a dry run to see what Snakemake plans to do without actually executing any commands:
    ```bash
    snakemake --dry-run
    # or a shorter version
    snakemake -n
    ```
2.  **Execute Locally**: Run the pipeline on your local machine, utilizing `X` CPU cores:
    ```bash
    snakemake -j X
    ```
3.  **Execute on HPC Cluster**: For cluster execution, the `--snakefile snakefile-hpc` option is used to specify the HPC-optimized workflow. Cluster parameters are passed via the `--cluster` flag (refer to the `hpc execution` section below for a detailed example):
    ```bash
    snakemake --snakefile snakefile-hpc --cluster "sbatch ..."
    ```
For more in-depth information, please refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

## execution

```
snakemake -j X # executing vaxpipe using X cpu cores
```

## hpc execution

For High-Performance Computing (HPC) environments, `vaxpipe` leverages Snakemake's cluster submission capabilities, typically via workload managers like SLURM. The `snakefile-hpc` is specifically optimized for this purpose.

Here's a breakdown of the command-line options:

```
snakemake --jobs 500 \
 --cores 1000 \
 --local-cores 1 \
 --snakefile snakefile-hpc \
 --latency-wait 120 \
 --cluster "sbatch --ntasks=1 \
                   --cpus-per-task=1 \
                   --job-name=vaxpipe_{rule}_{wildcards} \
                   --mem=4G \
                   --time=10:00:00 \
                   --error=logs/{rule}_{wildcards}.err \
                   --output=logs/{rule}_{wildcards}.out"
```

**Explanation of Options:**

-   `--jobs 500`: This flag tells Snakemake to submit up to 500 jobs to the cluster simultaneously. This is the maximum number of jobs that Snakemake will attempt to run in parallel on the HPC system.
-   `--cores 1000`: This specifies the total number of CPU cores available across all cluster jobs. Snakemake uses this to manage resource allocation, ensuring that the total core requests from submitted jobs do not exceed this limit. Note that individual jobs will request resources as defined in the `--cluster` command.
-   `--local-cores 1`: This reserves 1 CPU core for local tasks (e.g., Snakemake's internal processing or very small, quick jobs that are not submitted to the cluster). This ensures Snakemake itself has sufficient resources to manage the workflow.
-   `--snakefile snakefile-hpc`: This explicitly tells Snakemake to use the `snakefile-hpc` workflow definition, which is tailored for cluster execution and may include specific resource requests or directives for job submission.
-   `--latency-wait 120`: This sets a grace period (in seconds) that Snakemake waits for output files to appear after a job completes. This is particularly useful in networked file systems where there might be a delay in file synchronization, preventing Snakemake from prematurely marking a job as failed if its output isn't immediately visible.
-   `--cluster "sbatch ..."`: This is the most critical option for HPC execution. It passes the specified `sbatch` command directly to the cluster's workload manager (SLURM in this example) for each job submission. The parameters within the quotes define the resources and properties of each individual job:
    *   `--ntasks=1`: Requests 1 task per job.
    *   `--cpus-per-task=1`: Requests 1 CPU core per task. Thus, each job will run on a single CPU core.
    *   `--job-name=vaxpipe_{rule}_{wildcards}`: Assigns a dynamic name to each job, incorporating the Snakemake rule name and any wildcards, which helps in tracking jobs on the cluster.
    *   `--mem=4G`: Requests 4 Gigabytes of memory for each job.
    *   `--time=10:00:00`: Sets a time limit of 10 hours for each job's execution.
    *   `--error=logs/{rule}_{wildcards}.err`: Redirects standard error output to a specific error log file, dynamically named based on the rule and wildcards.
    *   `--output=logs/{rule}_{wildcards}.out`: Redirects standard output to a specific output log file, also dynamically named.

Users should adjust the values for `--jobs`, `--cores`, and particularly the `--cluster` parameters (e.g., memory, time, cpus-per-task, job queues) to match their specific HPC environment and the requirements of their protein design tasks.