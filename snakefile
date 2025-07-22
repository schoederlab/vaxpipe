# Snakefile

import os

# Configuration
ROSETTA_DIR = os.environ.get("ROSETTA", "/home/iwe14/ML_Rosetta/rosetta.source.release-362")
WORKDIR = "/home/iwe14/Documents/test"

SAMPLES = ["MERS"]
TAG = "test1"

rule all:
    input:
        expand(f"{WORKDIR}/{{sample}}_relax_{TAG}.sc", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_2.symm", sample=SAMPLES)

rule make_symmdef_file1:
    input:
        pdb = f"{WORKDIR}/{{sample}}.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}.symm",
        pdb = f"{WORKDIR}/{{sample}}_INPUT.pdb"
    shell:
        """
        {ROSETTA_DIR}/main/source/src/apps/public/symmetry/make_symmdef_file.pl \
        -p {input.pdb} -a A -i B > {output.symm}
        """

rule relax:
    input:
        pdb = f"{WORKDIR}/{{sample}}_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}.symm"
    output:
        scorefile = f"{WORKDIR}/{{sample}}_relax_{TAG}.sc",
        relaxed_pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001.pdb"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/relax.linuxgccrelease \
        -s {input.pdb} \
        -constrain_relax_to_start_coords \
        -nstruct 1 \
        -multiple_processes_writing_to_one_directory \
        -out:prefix relax_ \
        -out:path:all {WORKDIR} \
        -scorefile {output.scorefile} \
        -symmetry_definition {input.symm}
        """

rule make_symmdef_file2:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_INPUT.pdb"
    shell:
        """
        {ROSETTA_DIR}/main/source/src/apps/public/symmetry/make_symmdef_file.pl \
        -p {input.pdb} -a A -i B > {output.symm}
        """
