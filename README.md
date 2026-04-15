# vaxpipe

**vaxpipe** is a Snakemake-based pipeline designed to automate structural modeling of symmetric protein complexes using the [Rosetta](https://www.rosettacommons.org/) molecular modeling suite. 

It streamlines the process of the prediction of possible amino acid mutations using deep learning based methods, such as ProteinMPNN and ESM and database driven approaches, such as Rosetta FastDesign, with a focus on shape complementarity of protein interfaces.


### Operational Guidelines

- **Input PDB Structure**: The pipeline expects an input PDB file (e.g., `test/3ft7.pdb`) which represents a symmetric protein complex.
- **Rosetta Installation**: A functional Rosetta installation compiled with PyTorch and TensorFlow libraries is required for the deep learning-based predictions. For cluster usage, a singularity image can be downloaded from Docker Hub (see details below). Ensure that the `ROSETTA` path is correctly specified in the `config.yaml`.
- **Snakemake**: The pipeline is orchestrated using Snakemake. While users should ideally be familiar with Snakemake, this documentation aims to provide sufficient guidance for basic local and HPC execution. Detailed instructions are provided within the relevant sections.
- **Configuration**: The `config.yaml` file must be updated with appropriate paths to Rosetta binaries, the name of the pdb file, and input/output directories. This file is critical for the pipeline's correct operation.
- **Output Files**: Intermediate and final output files will be generated in a structured directory format, as defined in the `snakefile`. If the pipeline run successfully, several .png files will be present in the output directory, which represent both frequency and energy evaluation of tested mutations.
- **Resource Management**: For HPC execution, users are responsible for configuring SLURM or other cluster submission parameters within the `snakefile-hpc` to match their system's requirements and resource availability.

### Local Execution Requirements

- Python ≥ 3.11, biopython, matplotlib, numpy, pandas (installed via `uv sync`)
- [Snakemake ≥ 9.0](https://snakemake.readthedocs.io/en/v9.3.0/)
- Rosetta compiled with pytorch and tensorflow libraries. A detailed information on how to compile Rosetta with pytorch and tensorflow support can be found [here](https://docs.rosettacommons.org/docs/latest/build_documentation/Building-Rosetta-with-TensorFlow-and-Torch)
- ESM model (will be downloaded automatically)

### HPC Requirements
- Python ≥ 3.11, biopython, matplotlib, numpy, pandas (installed via `uv sync`)
- Snakemake ≥ 9.0
- **singularity**
- We provide a snakemake file that relies on a Rosetta docker/singularity image. The image is available on Docker Hub and can be pulled using singularity: `singularity pull docker://rosettacommons/rosetta:ml-387` (used for cluster execution).
- ESM model [download](https://git.iwe-lab.de/moritzertelt/ML_graphs/-/tree/main/tensorflow_graphs/ESM/esm2_t33_650M_UR50D). Currently, the pipeline just accepts this ESM model. **The model needs to be downloaded and copied into the repository path.**

### Repository Structure

```
vaxpipe/
├── snakefile             # local workflow
├── snakefile-hpc         # HPC workflow (SLURM)
├── config.yaml           # paths and execution settings
├── pyproject.toml        # Python dependencies (install with uv sync)
├── profiles/
│   └── slurm/
│       └── config.yaml   # Snakemake SLURM profile
├── input_files/          # Rosetta XML protocols and helper scripts
└── test/                 # test input PDB for validation
```
---

### Installation

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) if you haven't already, then create the environment and install all dependencies from the included `pyproject.toml`:

```bash
uv sync
```

This creates a `.venv` in the repository root. Set `python` in `config.yaml` to the resulting binary:

```yaml
python: /path/to/vaxpipe/.venv/bin/python
```

**Change the remaining paths in `config.yaml` to match your setup.**

---

### Understanding and Executing the Pipeline with Snakemake

`vaxpipe` is built with Snakemake, a powerful workflow management system that helps you create reproducible and scalable data analysis pipelines. Here's a brief guide to get you started:

#### Core Concepts:
- **Rules**: Define how output files are generated from input files. Each rule specifies a command to be run.
- **Wildcards**: Allow rules to be applied to many files (e.g., `data/{sample}.txt` where `{sample}` is a wildcard).
- **DAG (Directed Acyclic Graph)**: Snakemake automatically builds a DAG of jobs, ensuring that all dependencies are met before a rule is executed.
- **Execution**: Snakemake determines which rules need to be run based on the desired output files and the current state of your input files.

#### Essential Commands:

1.  **Dry Run**: Always start with a dry run to see what Snakemake plans to do without actually executing any commands (locally):
    ```bash
    snakemake --dry-run 
    # or a shorter version
    snakemake -n
    ```
2.  **Execute Locally**: Run the pipeline on your local machine, utilizing `X` CPU cores and opting into Conda-based software deployment (replacement for the deprecated `--use-conda` flag):
    ```bash
    snakemake --cores X --software-deployment-method conda
    ```
3.  **Execute on HPC Cluster**: For cluster execution, use the HPC-optimized snakefile together with Snakemake’s executor interface (replacement for the deprecated `--cluster` flag). See the `hpc execution` section below for a detailed SLURM example:
    ```bash
    snakemake --snakefile snakefile-hpc --cluster "sbatch ..."
    ```
For more in-depth information, please refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

## execution

```
snakemake --cores X --software-deployment-method conda  # executing vaxpipe using X CPU cores
```

### HPC execution

For High-Performance Computing (HPC) environments, `vaxpipe` leverages Snakemake's cluster submission capabilities, typically via workload managers like SLURM. The `snakefile-hpc` is specifically optimized for this purpose.

All cluster options are kept in `profiles/slurm/config.yaml` instead of the command line. To run the pipeline on HPC simply execute:

```bash
snakemake --profile profiles/slurm
```

Adjust `profiles/slurm/config.yaml` to match your cluster (partition, memory, time limits, etc.) before submitting. Per-rule resource overrides (e.g. higher memory for ESM) are set via the `resources` block in `snakefile-hpc` and picked up automatically by the SLURM executor.