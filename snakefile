# Snakefile

import os
import yaml

configfile: "config.yaml"

# Configuration
ROSETTA_DIR = config["rosettadir"]
INPUTDIR = config["inputdir"]
WORKDIR = config["workdir"]
SAMPLES = config["samples"]

#wildcards
variants = ["esm", "pmpnn", "indes"]
modes = ["design", "control"]
    
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
        expand(f"{WORKDIR}/esm_relax_{{sample}}_new_0001_INPUT_{{i}}.pdb", sample=SAMPLES, i=[f"{i:04d}" for i in range(1, 51)]),
        expand(f"{WORKDIR}/{{sample}}_esm.fasta", sample=SAMPLES),
        #mpnn
        expand(f"{WORKDIR}/{{sample}}_mpnn_probs.weights", sample=SAMPLES),
        expand(f"{WORKDIR}/pmpnn_relax_{{sample}}_new_0001_INPUT_{{i}}.pdb", sample=SAMPLES, i=[f"{i:04d}" for i in range(1, 51)]),
        expand(f"{WORKDIR}/{{sample}}_pmpnn.fasta", sample=SAMPLES),
        #interface design
        expand(f"{WORKDIR}/indes_relax_{{sample}}_new_0001_INPUT_{{i}}.pdb", sample=SAMPLES, i=[f"{i:04d}" for i in range(1, 51)]),
        expand(f"{WORKDIR}/{{sample}}_indes.fasta", sample=SAMPLES),
        #analysis
        expand(f"{WORKDIR}/{{sample}}_WT.fasta", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_{{variant}}.fasta", sample=SAMPLES, variant=variants),
        expand(f"{WORKDIR}/{{sample}}_{{variant}}_frequency.png", sample=SAMPLES, variant=variants),
        expand(f"{WORKDIR}/{{sample}}_{{variant}}_frequency.csv", sample=SAMPLES, variant=variants),
        expand(f"{WORKDIR}/{{sample}}_{{variant}}_mutations.txt", sample=SAMPLES, variant=variants),
        #design
        expand(f"{WORKDIR}/{{sample}}_{{mode}}_{{variant}}/DONE.txt", sample=SAMPLES, mode=modes, variant=variants),
        #plotting
        expand(f"{WORKDIR}/{{sample}}_energydifference_{{variant}}.png", sample=SAMPLES, variant=variants)

ruleorder: clean_pdb > make_symmdef_file1 > rename_file > relax > make_symmdef_file2 > run_esm > run_pmpnn > esm_sampling \
> pmpnn_sampling > interface_design > get_wt_fasta > get_fasta_from_pdbs > plot_frequencies > get_mutation_list \
> run_design_or_control > plot_energy

rule clean_pdb:
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
    input:
        pdb = lambda wc: f"{WORKDIR}/{wc.sample}_clean_0001.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}.symm",
        pdb = f"{WORKDIR}/{{sample}}_clean_0001_INPUT.pdb"
    shell:
        """
        {ROSETTA_DIR}/main/source/src/apps/public/symmetry/make_symmdef_file.pl \
        -p {input.pdb} -a A -i B > {output.symm}
        """

rule rename_file:
    input:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001_INPUT.pdb"
    output:
        pdbs = f"{WORKDIR}/{{sample}}_new.pdb"
    shell:
        """
        mv {input.pdb} {output.pdbs}
        """

rule relax:
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
        -nstruct 50 \
        -multiple_processes_writing_to_one_directory \
        -out:prefix relax_ \
        -out:path:all {WORKDIR} \
        -symmetry_definition {input.symm}
        """

rule make_symmdef_file2:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    shell:
        """
        {ROSETTA_DIR}/main/source/src/apps/public/symmetry/make_symmdef_file.pl \
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
        -s {input.pdb} \
        -parser:script_vars weights={output.weights} \
        -beta \
        -auto_download \
        -out:path:all {WORKDIR}
        """

rule esm_sampling:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        weights = f"{WORKDIR}/{{sample}}_esm_probs.weights"
    output:
        pdbs = f"{WORKDIR}/esm_relax_{{sample}}_new_0001_INPUT_{{i}}.pdb"
    params:
        protocol = f"{INPUTDIR}/esm/sample_mutations.xml",
        resfile = f"{INPUTDIR}/esm/resfile.resfile"
    threads: 1
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -s {input.pdb} \
        -parser:script_vars sym={input.symm} \
        -parser:script_vars weights={input.weights} \
        -parser:script_vars resfile={params.resfile} \
        -nstruct 50 \
        -out:path:all {WORKDIR} \
        -out:prefix esm_ \
        -beta \
        -overwrite > {WORKDIR}/esm_sampling.log
        """

rule run_pmpnn:
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
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        weights = f"{WORKDIR}/{{sample}}_mpnn_probs.weights"
    output:
        pdbs = f"{WORKDIR}/pmpnn_relax_{{sample}}_new_0001_INPUT_{{i}}.pdb"
    params:
        protocol = f"{INPUTDIR}/pmpnn/sample_mutations.xml",
        resfile = f"{INPUTDIR}/pmpnn/resfile.resfile"
    threads: 1
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -s {input.pdb} \
        -parser:script_vars sym={input.symm} \
        -parser:script_vars weights={input.weights} \
        -parser:script_vars resfile={params.resfile} \
        -out:path:all {WORKDIR} \
        -out:prefix pmpnn_ \
        -nstruct 50 \
        -beta \
        -overwrite > {WORKDIR}/pmpnn_sampling.log
        """

rule interface_design:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm"
    output:
        pdb = f"{WORKDIR}/indes_relax_{{sample}}_new_0001_INPUT_{{i}}.pdb"
    params:
        protocol = f"{INPUTDIR}/interface-design/sym_design.xml"
    threads: 1 
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -s {input.pdb} \
        -parser:script_vars sym={input.symm} \
        -beta \
        -nstruct 50 \
        -out:path:all {WORKDIR} \
        -out:prefix indes_ \
        -overwrite > {WORKDIR}/interface_design.log
        """

rule get_fasta_from_pdbs:
    input:
        pdbs = lambda wc: [ 
            f"{WORKDIR}/{wc.variant}_relax_{wc.sample}_new_0001_INPUT_{i:04d}.pdb"
            for i in range(1, 51)
        ]
    output:
        fastafile = f"{WORKDIR}/{{sample}}_{{variant}}.fasta"
    params:
        script = f"{INPUTDIR}/get_fasta/get_multifasta_from_pdb_path.py"
    shell:
        """
        python {params.script} -p {input.pdbs} -c A -o {output.fastafile}
        """

rule get_wt_fasta:
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
    input:
        fastafile = f"{WORKDIR}/{{sample}}_{{variant}}.fasta",
        wtfile = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        figure = f"{WORKDIR}/{{sample}}_{{variant}}_frequency.png",
        csv = f"{WORKDIR}/{{sample}}_{{variant}}_frequency.csv"
    params:
        script = f"{INPUTDIR}/validate/plot_frequencies.py"
    shell:
        """
        python {params.script} -i {input.fastafile} -r {input.wtfile} -o {output.figure}
        """

rule get_mutation_list:
    input:
        csv = f"{WORKDIR}/{{sample}}_{{variant}}_frequency.csv"
    output:
        txt = f"{WORKDIR}/{{sample}}_{{variant}}_mutations.txt"
    params:
        script = f"{INPUTDIR}/validate/design-mutations.py"
    shell:
        """
        python {params.script} -i {input.csv} -o {output.txt} \
        """

rule run_design_or_control:
    input:
        txt = f"{WORKDIR}/{{sample}}_{{variant}}_mutations.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        designdir = directory(f"{WORKDIR}/{{sample}}_{{mode}}_{{variant}}/"),
        donefile = f"{WORKDIR}/{{sample}}_{{mode}}_{{variant}}/DONE.txt"
    wildcard_constraints:
        sample = r"[a-zA-Z0-9]+",
        mode = r"(design|control)",
        variant = r"[a-zA-Z0-9]+"
    params:
        script = lambda wildcards: (
            f"{INPUTDIR}/validate/command_multiple_design.sh"
            if wildcards.mode == "design"
            else f"{INPUTDIR}/validate/command_multiple_control.sh"
        ),
        xml = f"{INPUTDIR}/validate/design.v02.xml"
    shell:
        """
        bash {params.script} {input.txt} {output.designdir} {input.symfile} {input.pdb} {params.xml} {ROSETTA_DIR}
        touch {output.donefile}
        """

rule plot_energy:
    input:
        directory_control = f"{WORKDIR}/{{sample}}_control_{{variant}}/",
        directory_design = f"{WORKDIR}/{{sample}}_design_{{variant}}/"
    output:
        image = f"{WORKDIR}/{{sample}}_energydifference_{{variant}}.png"
    params:
        script = f"{INPUTDIR}/validate/plot_energies.py"
    shell:
        """
        python {params.script} -i1 {input.directory_design} -i2 {input.directory_control} -o {output.image}
        """