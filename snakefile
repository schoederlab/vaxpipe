#!/bin/python
# Snakefile

import os
import yaml

configfile: "config.yaml"

# Configuration
ROSETTA_DIR = config["rosettadir"]
INPUTDIR = config["inputdir"]
WORKDIR = config["workdir"]
SAMPLES = config["samples"]
ITERATIONS = [f"{i:04d}" for i in range(1, 6)]
MUTATIONS = [f"{m}" for m in range(1, 6)]
PROSS_TEMPS = config["pross_temps"]
PROSS_TEMPS_STR = [str(t) for t in PROSS_TEMPS]

#wildcards
variants = ["pmpnn","indes","esm","pross"]
analysis_variants = ["pmpnn","indes","esm"]
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
        expand(f"{WORKDIR}/esm/{{sample}}_esm_probs.weights", sample=SAMPLES),
        expand(f"{WORKDIR}/esm/esm_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb", sample=SAMPLES, i=ITERATIONS),
        #mpnn
        expand(f"{WORKDIR}/pmpnn/{{sample}}_mpnn_probs.weights", sample=SAMPLES),
        expand(f"{WORKDIR}/pmpnn/pmpnn_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb", sample=SAMPLES, i=ITERATIONS),
        #interface design
        expand(f"{WORKDIR}/indes/indes_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb", sample=SAMPLES, i=ITERATIONS),
        #expand(f"{WORKDIR}/{{variant}}_{{sample}}_{{i}}.log", variant=variants, sample=SAMPLES, i=ITERATIONS),
        #pross
        expand(f"{WORKDIR}/pross/{{sample}}.pssm", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}.cst", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}.hhr", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}.a3m", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}.psi", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}_resfiles_pross/designable_aa_resfile.{{t}}", sample=SAMPLES, t=PROSS_TEMPS),
        expand(f"{WORKDIR}/pross/{{sample}}_pross_design_{{t}}.sc", sample=SAMPLES, t=PROSS_TEMPS),
        expand(f"{WORKDIR}/pross/{{sample}}_pross_wt_{{t}}.sc", sample=SAMPLES, t=PROSS_TEMPS),
        #analysis
        expand(f"{WORKDIR}/{{sample}}_WT.fasta", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_{{variant}}.fasta", sample=SAMPLES, variant=analysis_variants),
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.png", sample=SAMPLES, variant=analysis_variants),
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.csv", sample=SAMPLES, variant=analysis_variants),
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{m}}.txt", sample=SAMPLES, variant=analysis_variants, m=MUTATIONS),
        #design
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{mode}}/{{m}}.sc", m=MUTATIONS, sample=SAMPLES, mode=modes, variant=analysis_variants),
        #plotting
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_energydifference_{{variant}}.png", sample=SAMPLES, variant=analysis_variants)

ruleorder: clean_pdb > make_symmdef_file1 > rename_file > relax > make_symmdef_file2 > get_wt_fasta > run_esm > run_pmpnn > esm_sampling > pmpnn_sampling > interface_design \
  > get_fasta_from_pdbs > plot_frequencies > get_mutation_list > run_design_or_control > plot_energy

rule clean_pdb:
    localrule: True
    input:
        pdb = f"{WORKDIR}/{{sample}}.pdb"
    output:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001.pdb"
    shell:
        """
        mkdir -p {WORKDIR}/esm
        mkdir -p {WORKDIR}/pmpnn
        mkdir -p {WORKDIR}/indes
        mkdir -p {WORKDIR}/pross
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
        weights = f"{WORKDIR}/esm/{{sample}}_esm_probs.weights"
    params:
        protocol = f"{INPUTDIR}/esm/run_esm_and_save.xml"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -parser:script_vars weights={output.weights} \
        -s {input.pdb} \
        -beta \
        -auto_download \
        -out:path:all {WORKDIR}/esm \
        -overwrite
        """

rule esm_sampling:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        weights = f"{WORKDIR}/esm/{{sample}}_esm_probs.weights"
    output:
        pdb = f"{WORKDIR}/esm/esm_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb",
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
            -out:path:all {WORKDIR}/esm \
            -overwrite > {WORKDIR}/esm/esm.log
        """

rule run_pmpnn:
    localrule: True
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        weights = f"{WORKDIR}/pmpnn/{{sample}}_mpnn_probs.weights"
    params:
        protocol = f"{INPUTDIR}/pmpnn/run_mpnn_and_save.xml"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -parser:script_vars weights={output.weights} \
        -s {input.pdb} \
        -beta \
        -out:path:all {WORKDIR}/pmpnn \
        -overwrite
        """

rule pmpnn_sampling:
    localrule: True
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        weights = f"{WORKDIR}/pmpnn/{{sample}}_mpnn_probs.weights"
    output:
        pdb = f"{WORKDIR}/pmpnn/pmpnn_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb",
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
            -out:path:all {WORKDIR}/pmpnn \
            -overwrite > {WORKDIR}/pmpnn/pmpnn.log
        """

rule interface_design:
    localrule: True
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm"
    output:
        pdb = f"{WORKDIR}/indes/indes_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb"
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
            -out:path:all {WORKDIR}/indes \
            -out:prefix indes_ \
            -out:suffix _{wildcards.i} \
            -out:path:all {WORKDIR}/indes \
            -overwrite > {WORKDIR}/indes/indes.log
        """

rule get_fasta_from_pdbs:
    localrule: True
    input:
        pdbs = lambda wildcards: expand(
            f"{WORKDIR}/{wildcards.analysis_variant}/{wildcards.analysis_variant}_relax_{wildcards.sample}_new_0001_INPUT_{{i}}_0001.pdb",
            i=ITERATIONS
        )
    output:
        fastafile = f"{WORKDIR}/{{sample}}_{{analysis_variant}}.fasta"
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

rule generate_PSSM_and_constraints:
    input:
        fastafile = f"{WORKDIR}/{{sample}}_WT.fasta",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        hhr = f"{WORKDIR}/pross/{{sample}}.hhr",
        a3m = f"{WORKDIR}/pross/{{sample}}.a3m",
        psi = f"{WORKDIR}/pross/{{sample}}.psi",
        pssm = f"{WORKDIR}/pross/{{sample}}.pssm",
        cst = f"{WORKDIR}/pross/{{sample}}.cst",
        hhr_log = f"{WORKDIR}/pross/{{sample}}_hhr.log"
    shell:
        """
        bash {INPUTDIR}/pross/pssm/generate_pssm.file \
        {input.fastafile} {output.hhr} {output.a3m} {output.psi} {output.hhr_log} {output.pssm}
        bash {INPUTDIR}/pross/filter/make_cst.sh {input.pdb} > {output.cst}
        """

rule filterscan:
    localrule: True
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        cst = f"{WORKDIR}/pross/{{sample}}.cst",
        pssm = f"{WORKDIR}/pross/{{sample}}.pssm",
        fasta = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        resfiles_path = f"{WORKDIR}/pross/{{sample}}_resfiles_pross/designable_aa_resfile.{{t}}",
    params:
        protocol = f"{INPUTDIR}/pross/filter/filterscan.xml",
        path = f"{WORKDIR}/pross/{{sample}}_resfiles_pross/designable_aa_resfile",
    shell:
        """
        for res in $(seq 1 $(grep -v '^>' {input.fasta} | tr -d '\n' | wc -c)); do
            {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
                -parser:protocol {params.protocol} \
                -s {input.pdb} \
                -parser:script_vars sym={input.symm} \
                -parser:script_vars pdb_reference={input.pdb} \
                -parser:script_vars cst_full_path={input.cst} \
                -parser:script_vars cst_value=0.4 \
                -parser:script_vars pssm_full_path={input.pssm} \
                -parser:script_vars resfiles_path={params.path} \
                -parser:script_vars current_res=$res \
                -out:path:all {WORKDIR}/pross \
                -beta \
                -overwrite > {WORKDIR}/pross/filterscan.log
        done
        """

rule pross_design:
    input:
        resfile = f"{WORKDIR}/pross/{{sample}}_resfiles_pross/designable_aa_resfile.{{t}}",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        cst = f"{WORKDIR}/pross/{{sample}}.cst",
        pssm = f"{WORKDIR}/pross/{{sample}}.pssm"
    output:
        sc = f"{WORKDIR}/pross/{{sample}}_pross_design_{{t}}.sc"
    params:
        protocol = f"{INPUTDIR}/pross/design/design.xml",
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -parser:script_vars pdb_reference={input.pdb} \
            -parser:script_vars cst_full_path={input.cst} \
            -parser:script_vars cst_value=0.4 \
            -parser:script_vars pssm_full_path={input.pssm} \
            -parser:script_vars in_resfile={input.resfile} \
            -overwrite > {WORKDIR}/pross/pross_design.log \
            -ignore_unrecognized_res \
            -use_input_sc \
            -use_occurrence_data \
            -out:file:scorefile {output.sc} \
            -out:path:all {WORKDIR}/pross \
            -out:prefix pross_design_{wildcards.t} \
            -beta \
        """

rule pross_design_wt:
    localrule: True
    input:
        resfile = f"{WORKDIR}/pross/{{sample}}_resfiles_pross/designable_aa_resfile.{{t}}",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        cst = f"{WORKDIR}/pross/{{sample}}.cst",
        pssm = f"{WORKDIR}/pross/{{sample}}.pssm"
    output:
        sc = f"{WORKDIR}/pross/{{sample}}_pross_wt_{{t}}.sc"
    params:
        protocol = f"{INPUTDIR}/pross/design/design_WT.xml",
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -parser:script_vars pdb_reference={input.pdb} \
            -parser:script_vars cst_full_path={input.cst} \
            -parser:script_vars cst_value=0.4 \
            -parser:script_vars pssm_full_path={input.pssm} \
            -parser:script_vars in_resfile={input.resfile} \
            -overwrite > {WORKDIR}/pross/pross_wt.log \
            -out:path:all {WORKDIR}/pross \
            -ignore_unrecognized_res \
            -use_input_sc \
            -use_occurrence_data \
            -out:file:scorefile {output.sc} \
            -out:prefix pross_wt_{wildcards.t} \
            -out:path:all {WORKDIR}/pross \
            -beta
        """

rule plot_frequencies:
    localrule: True
    input:
        fastafile = f"{WORKDIR}/{{sample}}_{{variant}}.fasta",
        wtfile = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        figure = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.png",
        csv = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.csv"
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
        csv = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.csv"
    output:
        out = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{m}}.txt"
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
        txt = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{m}}.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        sc =  f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{mode}}/{{m}}.sc"
    params:
        xml = f"{INPUTDIR}/validate/design.v02.xml",
        outdir = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{mode}}"
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
            -out:path:all {WORKDIR}/{wildcards.variant}/ \
            -overwrite \
            -nstruct 20 \
            -beta
        """

rule plot_energy:
    localrule: True
    input:
        control = lambda wildcards: expand(f"{WORKDIR}/{wildcards.variant}/{wildcards.sample}_{wildcards.variant}/control/{{m}}.sc", m=MUTATIONS),
        design = lambda wildcards: expand(f"{WORKDIR}/{wildcards.variant}/{wildcards.sample}_{wildcards.variant}/design/{{m}}.sc", m=MUTATIONS),
    output:
        image = f"{WORKDIR}/{{variant}}/{{sample}}_energydifference_{{variant}}.png"
    params:
        script = f"{INPUTDIR}/validate/plot_energies.py"
    shell:
        """
        python {params.script} \
            -i1 $(dirname {input.design[0]}) \
            -i2 $(dirname {input.control[0]}) \
            -o {output.image}
        """