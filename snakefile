#!/bin/python
# Snakefile

import os
import yaml

min_version("9.0")

configfile: "config.yaml"

# Configuration
ROSETTA_DIR = config["rosettadir"]
INPUTDIR = config["inputdir"]
WORKDIR = config["workdir"]
SAMPLES = config["samples"]
ITERATIONS = config["iterations"]
MUTATIONS = config["mutations"]

#wildcards
variants = ["indes","pmpnn", "esm"]
modes = ["design", "control"]
wildcard_constraints:
    i = r"\d{4}"

rule all:
    input:
        #preprocessing
        expand(f"{WORKDIR}/{{sample}}_clean_0001.pdb", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}.symm", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_new.pdb", sample=SAMPLES),
        expand(f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_2.symm", sample=SAMPLES),
        #esm
        expand(f"{WORKDIR}/{{sample}}_esm_probs.weights", sample=SAMPLES),
        expand(f"{WORKDIR}/esm_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb", sample=SAMPLES, i=ITERATIONS),
        #mpnn
        expand(f"{WORKDIR}/{{sample}}_mpnn_probs.weights", sample=SAMPLES),
        expand(f"{WORKDIR}/pmpnn_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb", sample=SAMPLES, i=ITERATIONS),
        #interface design
        expand(f"{WORKDIR}/indes_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb", sample=SAMPLES, i=ITERATIONS),
        #expand(f"{WORKDIR}/{{variant}}_{{sample}}_{{i}}.log", variant=variants, sample=SAMPLES, i=ITERATIONS),
        #analysis
        expand(f"{WORKDIR}/{{sample}}_WT.fasta", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_{{variant}}.fasta", sample=SAMPLES, variant=variants),
        expand(f"{WORKDIR}/{{sample}}_{{variant}}_frequency.png", sample=SAMPLES, variant=variants),
        expand(f"{WORKDIR}/{{sample}}_{{variant}}_frequency.csv", sample=SAMPLES, variant=variants),
        expand(f"{WORKDIR}/{{sample}}_{{variant}}/{{m}}.txt", sample=SAMPLES, variant=variants, m=MUTATIONS),
        #design
        expand(f"{WORKDIR}/{{sample}}_{{variant}}/{{mode}}/{{m}}.sc", m=MUTATIONS, sample=SAMPLES, mode=modes, variant=variants),
        #plotting
        expand(f"{WORKDIR}/{{sample}}_energydifference_{{variant}}.png", sample=SAMPLES, variant=variants)

ruleorder: clean_pdb > make_symmdef_file1 > rename_file > relax > make_symmdef_file2 > get_wt_fasta > run_pmpnn > pmpnn_sampling > interface_design \
  > get_fasta_from_pdbs > plot_frequencies > get_mutation_list > run_design_or_control > plot_energy

rule clean_pdb:
    localrule: True
    input:
        pdb = f"{WORKDIR}/{{sample}}.pdb"
    output:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001.pdb"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/score_jd2.pytorchtensorflow.linuxgccrelease \
        -renumber_pdb -ignore_unrecognized_res -s {input.pdb} \
        -out:pdb -out:suffix _clean -out:path:all {WORKDIR}
        """

rule make_symmdef_file1:
    localrule: True
    input:
        pdb = lambda wc: f"{WORKDIR}/{wc.sample}_clean_0001.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}.symm",
        pdb = f"{WORKDIR}/{{sample}}_clean_0001_INPUT.pdb"
    shell:
        """
        {INPUTDIR}/symmetry/make_symmdef_file.pl \
        -p {input.pdb} -a A -i B > {output.symm}
        """

rule rename_file:
    localrule: True
    input:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001_INPUT.pdb"
    output:
        pdbs = f"{WORKDIR}/{{sample}}_new.pdb"
    shell:
        """
        mv {input.pdb} {output.pdbs}
        """

rule relax:
    localrule: True
    input:
        pdb = f"{WORKDIR}/{{sample}}_new.pdb",
        symm = f"{WORKDIR}/{{sample}}.symm"
    output:
        relaxed_pdb = f"{WORKDIR}/relax_{{sample}}_new_0001.pdb"
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
        -symmetry_definition {input.symm}
        """

rule make_symmdef_file2:
    localrule: True
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    shell:
        """
        {INPUTDIR}/symmetry/make_symmdef_file.pl \
        -p {input.pdb} -a A -i B > {output.symm}
        """

rule run_esm:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        weights = f"{WORKDIR}/{{sample}}_esm_probs.weights"
    params:
        protocol = f"{INPUTDIR}/esm/run_esm_and_save.xml"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -parser:script_vars weights={output.weights} \
        -s {input.pdb} \
        -beta \
        -auto_download
        """

rule esm_sampling:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        weights = f"{WORKDIR}/{{sample}}_esm_probs.weights"
    output:
        pdb = f"{WORKDIR}/esm_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb",
    params:
        protocol = f"{INPUTDIR}/esm/sample_mutations.xml",
        resfile = f"{INPUTDIR}/esm/resfile.resfile",
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -parser:script_vars weights={input.weights} \
            -parser:script_vars resfile={params.resfile} \
            -out:pdb true \
            -out:path:all {WORKDIR} \
            -out:prefix esm_ \
            -out:suffix _{wildcards.i} \
            -beta \
            -overwrite > esm.log
        """

rule run_pmpnn:
    localrule: True
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        weights = f"{WORKDIR}/{{sample}}_mpnn_probs.weights"
    params:
        protocol = f"{INPUTDIR}/pmpnn/run_mpnn_and_save.xml"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -parser:script_vars weights={output.weights} \
        -s {input.pdb} \
        -beta \
        """

rule pmpnn_sampling:
    localrule: True
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        weights = f"{WORKDIR}/{{sample}}_mpnn_probs.weights"
    output:
        pdb = f"{WORKDIR}/pmpnn_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb",
    params:
        protocol = f"{INPUTDIR}/pmpnn/sample_mutations.xml",
        resfile = f"{INPUTDIR}/pmpnn/resfile.resfile",
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -parser:script_vars weights={input.weights} \
            -parser:script_vars resfile={params.resfile} \
            -out:pdb true \
            -out:path:all {WORKDIR} \
            -out:prefix pmpnn_ \
            -out:suffix _{wildcards.i} \
            -beta \
            -overwrite > pmpnn.log
        """

rule interface_design:
    localrule: True
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm"
    output:
        pdb = f"{WORKDIR}/indes_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb"
    params:
        protocol = f"{INPUTDIR}/interface-design/sym_design.xml",
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -beta \
            -out:pdb true \
            -out:path:all {WORKDIR} \
            -out:prefix indes_ \
            -out:suffix _{wildcards.i} \
            -overwrite > indes.log
        """

rule get_fasta_from_pdbs:
    localrule: True
    input:
        pdbs = lambda wildcards: expand(
            f"{WORKDIR}/{wildcards.variant}_relax_{wildcards.sample}_new_0001_INPUT_{{i}}_0001.pdb",
            i=ITERATIONS
        )
    output:
        fastafile = f"{WORKDIR}/{{sample}}_{{variant}}.fasta"
    params:
        script = f"{INPUTDIR}/get_fasta/get_multifasta_from_pdb_path.py"
    shell:
        """
        python {params.script} -p {input.pdbs} -c A -o {output.fastafile}
        """

rule get_wt_fasta:
    localrule: True
    input:
        pdb = f"{WORKDIR}/{{sample}}.pdb"
    output:
        fastafile = f"{WORKDIR}/{{sample}}_WT.fasta"
    params:
        script = f"{INPUTDIR}/get_fasta/get_multifasta_from_pdb_path.py"
    shell:
        """
        python {params.script} \
        -p {input.pdb} \
        -c A \
        -o {output.fastafile}
        """

rule plot_frequencies:
    localrule: True
    input:
        fastafile = f"{WORKDIR}/{{sample}}_{{variant}}.fasta",
        wtfile = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        figure = f"{WORKDIR}/{{sample}}_{{variant}}_frequency.png",
        csv = f"{WORKDIR}/{{sample}}_{{variant}}_frequency.csv"
    params:
        script = f"{INPUTDIR}/validate/plot_frequencies.py",
        mutations = len(MUTATIONS)
    shell:
        """
        python {params.script} -i {input.fastafile} -r {input.wtfile} -m {params.mutations} -o {output.figure}
        """

rule get_mutation_list:
    localrule: True
    input:
        csv = f"{WORKDIR}/{{sample}}_{{variant}}_frequency.csv"
    output:
        out = f"{WORKDIR}/{{sample}}_{{variant}}/{{m}}.txt"
    params:
        script = f"{INPUTDIR}/validate/design-mutations.py"
    shell:
        """
        mkdir -p $(dirname {output.out})
        python {params.script} -i {input.csv} -o {output.out} \
        """

rule run_design_or_control:
    localrule: True
    input:
        txt = f"{WORKDIR}/{{sample}}_{{variant}}/{{m}}.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        sc =  f"{WORKDIR}/{{sample}}_{{variant}}/{{mode}}/{{m}}.sc"
    params:
        xml = f"{INPUTDIR}/validate/design.v02.xml",
        outdir = f"{WORKDIR}/{{sample}}_{{variant}}/{{mode}}"
    shell:
        """
        MUTATION_LINE=$(cat {input.txt})
        
        mutpos=$(echo $MUTATION_LINE | cut -d'_' -f1)
        mutaa=$(echo $MUTATION_LINE | cut -d'_' -f2)

        mkdir -p {params.outdir}
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
            -parser:protocol {params.xml} \
            -parser:script_vars mutpos=$mutpos mut_aa=$mutaa protocol={wildcards.mode} symfile={input.symfile} \
            -in:file:s {input.pdb} \
            -corrections:beta_nov16 \
            -out:file:scorefile {output.sc} \
            -overwrite \
            -nstruct 20 \
            -beta
        """

rule plot_energy:
    localrule: True
    input:
        control = lambda wildcards: expand(f"{WORKDIR}/{wildcards.sample}_{wildcards.variant}/control/{{m}}.sc", m=MUTATIONS),
        design = lambda wildcards: expand(f"{WORKDIR}/{wildcards.sample}_{wildcards.variant}/design/{{m}}.sc", m=MUTATIONS),
    output:
        image = f"{WORKDIR}/{{sample}}_energydifference_{{variant}}.png"
    params:
        script = f"{INPUTDIR}/validate/plot_energies.py"
    shell:
        """
        python {params.script} \
            -i1 $(dirname {input.design[0]}) \
            -i2 $(dirname {input.control[0]}) \
            -o {output.image}
        """