# Snakefile

import os
import yaml

configfile: "config.yaml"

# Configuration
ROSETTA_DIR = config["rosettadir"]
WORKDIR = config["workdir"]
SAMPLES = config["samples"]
TAG = config["tag"]

#No. of designs
DESIGN_IDS =  list(range(1, 1001))

rule all:
    input:
        expand(f"{WORKDIR}/{{sample}}_relax_{TAG}.sc", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_2.symm", sample=SAMPLES),
        #interface design
        expand(f"{WORKDIR}/{{sample}}_in-des_{{design_id}}.pdb", sample=SAMPLES, design_id=range(1, 1001))

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
        {ROSETTA_DIR}/main/source/bin/relax.pytorchtensorflow.linuxgccrelease \
        -s {input.pdb} \
        -constrain_relax_to_start_coords \
        -beta \
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
        pdb = f"{WORKDIR}/{{sample}}_relaxed_symm.pdb"
    shell:
        """
        {ROSETTA_DIR}/main/source/src/apps/public/symmetry/make_symmdef_file.pl \
        -p {input.pdb} -a A -i B > {output.symm}
        """

rule run_esm:
    input:
        pdb = f"{WORKDIR}/{{sample}}_relaxed_symm.pdb"
    output:
        weights = f"{WORKDIR}/esm_probs.weights"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol ./input_files/esm/run_esm_and_save.xml \
        -s {input.pdb} \
        -beta \
        -auto_download \
        """

rule esm_sampling:
    input:
        pdb = f"{WORKDIR}/{{sample}}_relaxed_symm.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm"
    output:
        pdbs = expand(f"{WORKDIR}/{{sample}}_esm_{{i}}.pdb", sample=SAMPLES, i=range(1, 201))
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol ./input_files/esm/sample_mutations.xml \
        -s {input.pdb} \
        -parser:script_vars sym={input.symm} \
        -parser:script_vars resfile=./input_files/esm/resfile.resfile \
        -out:prefix esm_ \
        -nstruct 200 \
        -beta \
        -overwrite > esm_sampling.log
        """

rule run_pmpnn:
    input:
        pdb = f"{WORKDIR}/{{sample}}_relaxed_symm.pdb"
    output:
        weights = f"{WORKDIR}/mpnn_probs.weights"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol ./input_files/pmpnn/run_mpnn_and_save.xml \
        -s {input.pdb} \
        -beta \
        """

rule pmpnn_sampling:
    input:
        pdb = f"{WORKDIR}/{{sample}}_relaxed_symm.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm"
    output:
        pdbs = expand(f"{WORKDIR}/{{sample}}_pmpnn_{{i}}.pdb", sample=SAMPLES, i=range(1, 201))
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol ./input_files/pmpnn/sample_mutations.xml \
        -s {input.pdb} \
        -parser:script_vars sym={input.symm} \
        -parser:script_vars resfile=./input_files/esm/resfile.resfile \
        -out:prefix pmpnn_ \
        -nstruct 200 \
        -beta \
        -overwrite > pmpnn_sampling.log
        """

rule interface_design:
    input:
        pdb = f"{WORKDIR}/{{sample}}_relaxed_symm.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm"
    output:
        pdb = f"{WORKDIR}/{{sample}}_in-des_{{design_id}}.pdb"
    params:
        protocol = f"./input_files/interface_desing/sym_design.xml"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -s {input.pdb} \
        -parser:script_vars sym={input.symm} \
        -parser:script_vars resfile=./input_files/esm/resfile.resfile \
        -beta \
        -out:suffix {wildcards.design_id} \
        -overwrite > interface_design.log
        """    