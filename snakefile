# Snakefile

import os
import yaml

configfile: "config.yaml"

# Configuration
ROSETTA_DIR = config["rosettadir"]
INPUTDIR = config["inputdir"]
WORKDIR = config["workdir"]
SAMPLES = config["samples"]
TAG = config["tag"]

rule all:
    input:
        #preprocessing
        expand(f"{WORKDIR}/{{sample}}_clean_0001.pdb", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}.symm", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_INPUT.pdb", sample=SAMPLES),
        expand(f"{WORKDIR}/relax_{{sample}}_INPUT_0001.pdb", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_2.symm", sample=SAMPLES),
        #esm
        expand(f"{WORKDIR}/{{sample}}_esm_probs.weights", sample=SAMPLES),
        expand(f"{WORKDIR}/esm_relax_{{sample}}_INPUT_0001_symm_{{i}}.pdb", sample=SAMPLES, i=[f"{i:04d}" for i in range(1, 12)]),
        expand(f"{WORKDIR}/{{sample}}_esm.fasta", sample=SAMPLES),
        #mpnn
        expand(f"{WORKDIR}/{{sample}}_mpnn_probs.weights", sample=SAMPLES),
        expand(f"{WORKDIR}/mpnn_relax_{{sample}}_INPUT_0001_symm_{{i}}.pdb", sample=SAMPLES, i=[f"{i:04d}" for i in range(1, 12)]),
        expand(f"{WORKDIR}/{{sample}}_pmpnn.fasta", sample=SAMPLES),
        #interface design
        expand(f"{WORKDIR}/indes_relax_{{sample}}_INPUT_0001_symm_{{i}}.pdb", sample=SAMPLES, i=[f"{i:04d}" for i in range(1, 12)]),
        expand(f"{WORKDIR}/{{sample}}_in-des.fasta", sample=SAMPLES),
        #analysis
        expand(f"{WORKDIR}/{{sample}}_WT.fasta", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_esm_frequency.png", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_mpnn_frequency.png", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_indes_frequency.png", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_esm_frequency.csv", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_mpnn_frequency.csv", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_indes_frequency.csv", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_esm_mutations.txt", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_mpnn_mutations.txt", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_indes_mutations.txt", sample=SAMPLES),
        #design
        expand(f"{WORKDIR}/{{sample}}_design_esm/", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_control_esm/", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_design_mpnn/", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_control_mpnn/", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_design_indes/", sample=SAMPLES),
        expand(f"{WORKDIR}/{{sample}}_control_indes/", sample=SAMPLES)

rule clean_pdb:
    input:
        pdb = f"{WORKDIR}/{{sample}}.pdb"
    output:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001.pdb"
    shell:
        """
        {ROSETTA_DIR}//main/source/bin/score_jd2.pytorchtensorflow.linuxgccrelease \
        -renumber_pdb -ignore_unrecognized_res -s {input.pdb} \
        -out:pdb -out:suffix _clean -out:path:all {WORKDIR}
        """

rule make_symmdef_file1:
    input:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001.pdb"
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
        pdb = f"{WORKDIR}/{{sample}}_INPUT.pdb"
    wildcard_constraints:
        sample = "(?!.*(_clean_0001)).+"
    shell:
        """
        mv {input.pdb} {output.pdb}
        """

rule relax:
    input:
        pdb = f"{WORKDIR}/{{sample}}_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}.symm"
    output:
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
        -symmetry_definition {input.symm}
        """

rule make_symmdef_file2:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb"
    shell:
        """
        {ROSETTA_DIR}/main/source/src/apps/public/symmetry/make_symmdef_file.pl \
        -p {input.pdb} -a A -i B > {output.symm}
        """

rule run_esm:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb"
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
        """

rule esm_sampling:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        weights=f"{WORKDIR}/{{sample}}_esm_probs.weights"
    output:
        pdbs = f"{WORKDIR}/esm_relax_{{sample}}_INPUT_0001_symm_{{i}}.pdb"
    params:
        protocol = f"{INPUTDIR}/esm/sample_mutations.xml",
        resfile = f"{INPUTDIR}/esm/resfile.resfile"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -s {input.pdb} \
        -parser:script_vars sym={input.symm} \
        -parser:script_vars weights={input.weights} \
        -parser:script_vars resfile={params.resfile} \
        -nstruct 11 \
        -out:path:all {WORKDIR} \
        -out:prefix esm_ \
        -beta \
        -overwrite > {WORKDIR}/esm_sampling.log
        """

rule run_pmpnn:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb"
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
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        weights=f"{WORKDIR}/{{sample}}_mpnn_probs.weights"
    output:
        pdbs = f"{WORKDIR}/mpnn_relax_{{sample}}_INPUT_0001_symm_{{i}}.pdb"
    params:
        protocol = f"{INPUTDIR}/pmpnn/sample_mutations.xml",
        resfile = f"{INPUTDIR}/pmpnn/resfile.resfile"
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -s {input.pdb} \
        -parser:script_vars sym={input.symm} \
        -parser:script_vars weights={input.weights} \
        -parser:script_vars resfile={params.resfile} \
        -out:path:all {WORKDIR} \
        -out:prefix mpnn_ \
        -nstruct 11 \
        -beta \
        -overwrite > {WORKDIR}/pmpnn_sampling.log
        """

rule interface_design:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm"
    output:
        pdb = f"{WORKDIR}/indes_relax_{{sample}}_INPUT_0001_symm_{{i}}.pdb"
    params:
        protocol = f"{INPUTDIR}/interface-design/sym_design.xml" 
    shell:
        """
        {ROSETTA_DIR}/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol {params.protocol} \
        -s {input.pdb} \
        -parser:script_vars sym={input.symm} \
        -beta \
        -nstruct 11 \
        -out:path:all {WORKDIR} \
        -out:prefix indes_ \
        -overwrite > {WORKDIR}/interface_design.log
        """

rule get_fasta_from_esm:
    input:
        pdbs = expand(f"{WORKDIR}/esm_relax_{{sample}}_INPUT_0001_symm_{{i}}.pdb", sample=SAMPLES, i=[f"{i:04d}" for i in range(1, 12)])
    output:
        fastafile = f"{WORKDIR}/{{sample}}_esm.fasta"
    params:
        script = f"{INPUTDIR}/get_fasta/get_multifasta_from_pdb_path.py"
    shell:
        """
        python {params.script} \
        -p {input.pdbs} \
        -c A \
        -o {output.fastafile}
        """

rule get_fasta_from_pmpnn:
    input:
        pdbs = expand(f"{WORKDIR}/mpnn_relax_{{sample}}_INPUT_0001_symm_{{i}}.pdb", sample=SAMPLES, i=[f"{i:04d}" for i in range(1, 12)])
    output:
        fastafile = f"{WORKDIR}/{{sample}}_pmpnn.fasta"
    params:
        script = f"{INPUTDIR}/get_fasta/get_multifasta_from_pdb_path.py"
    shell:
        """
        python {params.script} \
        -p {input.pdbs} \
        -c A \
        -o {output.fastafile}
        """

rule get_fasta_in_des:
    input:
        pdbs = expand(f"{WORKDIR}/indes_relax_{{sample}}_INPUT_0001_symm_{{i}}.pdb", sample=SAMPLES, i=[f"{i:04d}" for i in range(1, 12)])
    output:
        fastafile = f"{WORKDIR}/{{sample}}_in-des.fasta"
    params:
        script = f"{INPUTDIR}/get_fasta/get_multifasta_from_pdb_path.py"
    shell:
        """
        python {params.script} \
        -p {input.pdbs} \
        -c A \
        -o {output.fastafile}
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

rule plot_esm_frequencies:
    input:
        fastafile = f"{WORKDIR}/{{sample}}_esm.fasta",
        wtfile = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        figure = f"{WORKDIR}/{{sample}}_esm_frequency.png",
        csv = f"{WORKDIR}/{{sample}}_esm_frequency.csv"
    params:
        script = f"{INPUTDIR}/validate/plot_frequencies.py"
    shell:
        """
        python {params.script} \
        -i {input.fastafile} \
        -r {input.wtfile} \
        -o {output.figure}
        """

rule plot_mpnn_frequencies:
    input:
        fastafile = f"{WORKDIR}/{{sample}}_pmpnn.fasta",
        wtfile = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        figure = f"{WORKDIR}/{{sample}}_mpnn_frequency.png",
        csv = f"{WORKDIR}/{{sample}}_mpnn_frequency.csv"
    params:
        script = f"{INPUTDIR}/validate/plot_frequencies.py"
    shell:
        """
        python {params.script} \
        -i {input.fastafile} \
        -r {input.wtfile} \
        -o {output.figure}
        """

rule plot_indes_frequencies:
    input:
        fastafile = f"{WORKDIR}/{{sample}}_in-des.fasta",
        wtfile = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        figure = f"{WORKDIR}/{{sample}}_indes_frequency.png",
        csv = f"{WORKDIR}/{{sample}}_indes_frequency.csv"
    params:
        script = f"{INPUTDIR}/validate/plot_frequencies.py"
    shell:
        """
        python {params.script} \
        -i {input.fastafile} \
        -r {input.wtfile} \
        -o {output.figure}
        """

rule get_esm_mutation_list:
    input:
        csv = f"{WORKDIR}/{{sample}}_esm_frequency.csv"
    output:
        txt = f"{WORKDIR}/{{sample}}_esm_mutations.txt"
    params:
        script = f"{INPUTDIR}/validate/design-mutations.py"
    shell:
        """
        python {params.script} \
        -i {input.csv} \
        -o {output.txt} \
        """

rule get_pmpnn_mutation_list:
    input:
        csv = f"{WORKDIR}/{{sample}}_mpnn_frequency.csv"
    output:
        txt = f"{WORKDIR}/{{sample}}_mpnn_mutations.txt"
    params:
        script = f"{INPUTDIR}/validate/design-mutations.py"
    shell:
        """
        python {params.script} \
        -i {input.csv} \
        -o {output.txt} \
        """

rule get_indes_mutation_list:
    input:
        csv = f"{WORKDIR}/{{sample}}_indes_frequency.csv"
    output:
        txt = f"{WORKDIR}/{{sample}}_indes_mutations.txt"
    params:
        script = f"{INPUTDIR}/validate/design-mutations.py"
    shell:
        """
        python {params.script} \
        -i {input.csv} \
        -o {output.txt} \
        """

rule design_esm:
    input:
        txt = f"{WORKDIR}/{{sample}}_esm_mutations.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb"
    output:
        designdir = f"{WORKDIR}/{{sample}}_design_esm/"
    params:
        script = f"{INPUTDIR}/validate/command_multiple_design.sh",
        xml = f"{INPUTDIR}/validate/design.v02.xml"
        
    shell:
        """
        bash {params.script} {input.txt} {output.designdir} {input.symfile} {input.pdb} {params.xml}
        """

rule control_esm:
    input:
        txt = f"{WORKDIR}/{{sample}}_esm_mutations.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb"
    output:
        designdir = f"{WORKDIR}/{{sample}}_control_esm/"
    params:
        script = f"{INPUTDIR}/validate/command_multiple_control.sh",
        xml = f"{INPUTDIR}/validate/design.v02.xml"
    shell:
        """
        bash {params.script} {input.txt} {output.designdir} {input.symfile} {input.pdb} {params.xml}
        """

rule design_mpnn:
    input:
        txt = f"{WORKDIR}/{{sample}}_mpnn_mutations.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb"
    output:
        designdir = f"{WORKDIR}/{{sample}}_design_mpnn/"
    params:
        script = f"{INPUTDIR}/validate/command_multiple_design.sh",
        xml = f"{INPUTDIR}/validate/design.v02.xml"
        
    shell:
        """
        bash {params.script} {input.txt} {output.designdir} {input.symfile} {input.pdb} {params.xml}
        """

rule control_mpnn:
    input:
        txt = f"{WORKDIR}/{{sample}}_mpnn_mutations.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb"
    output:
        designdir = f"{WORKDIR}/{{sample}}_control_mpnn/"
    params:
        script = f"{INPUTDIR}/validate/command_multiple_control.sh",
        xml = f"{INPUTDIR}/validate/design.v02.xml"
    shell:
        """
        bash {params.script} {input.txt} {output.designdir} {input.symfile} {input.pdb} {params.xml}
        """

rule design_indes:
    input:
        txt = f"{WORKDIR}/{{sample}}_indes_mutations.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb"
    output:
        designdir = f"{WORKDIR}/{{sample}}_design_indes/"
    params:
        script = f"{INPUTDIR}/validate/command_multiple_design.sh",
        xml = f"{INPUTDIR}/validate/design.v02.xml"
        
    shell:
        """
        bash {params.script} {input.txt} {output.designdir} {input.symfile} {input.pdb} {params.xml}
        """

rule control_indes:
    input:
        txt = f"{WORKDIR}/{{sample}}_indes_mutations.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_INPUT_0001_symm.pdb"
    output:
        designdir = f"{WORKDIR}/{{sample}}_control_indes/"
    params:
        script = f"{INPUTDIR}/validate/command_multiple_control.sh",
        xml = f"{INPUTDIR}/validate/design.v02.xml"
    shell:
        """
        bash {params.script} {input.txt} {output.designdir} {input.symfile} {input.pdb} {params.xml}
        """