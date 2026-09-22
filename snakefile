# Snakefile

import os
import yaml

configfile: "config.yaml"

# Configuration
ROSETTA_DIR = config["rosettadir"]
INPUTDIR = config["inputdir"]
WORKDIR = config["workdir"]
SAMPLES = config["samples"]
UNIREF_DB = config["uniref_db"]
PROSS_TEMPS = config["pross_temps"]

#wildcards
analysis_variants = ["esm","indes","pmpnn"]
modes = ["design", "control"]
# Number of designs per branch and number of mutations carried into validation.
# Keep these small for a local test run, large for a production run on the cluster.
ITERATIONS = [f"{i:04d}" for i in range(1, config["iterations"] + 1)]
MUTATIONS = [f"{m}" for m in range(1, config["mutations"] + 1)]

# Wildcards must not swallow path separators, otherwise e.g. {variant}/{sample}_{variant}
# can be matched in several ways and rules become ambiguous.
wildcard_constraints:
    i = r"\d{4}",
    m = r"\d+",
    sample = r"[^/]+",
    variant = r"[^/]+",
    mode = r"design|control",
    t = r"[^/]+"

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
        #pross
        #expand(f"{WORKDIR}/pross/{{sample}}.pssm", sample=SAMPLES),
        #expand(f"{WORKDIR}/pross/{{sample}}.cst", sample=SAMPLES),
        #expand(f"{WORKDIR}/pross/{{sample}}.hhr", sample=SAMPLES),
        #expand(f"{WORKDIR}/pross/{{sample}}.a3m", sample=SAMPLES),
        #expand(f"{WORKDIR}/pross/{{sample}}.psi", sample=SAMPLES),
        #expand(f"{WORKDIR}/pross/{{sample}}_resfiles_pross/designable_aa_resfile.{{t}}", sample=SAMPLES, t=PROSS_TEMPS),
        #expand(f"{WORKDIR}/pross/{{sample}}_pross_design_{{t}}.sc", sample=SAMPLES, t=PROSS_TEMPS),
        #expand(f"{WORKDIR}/pross/{{sample}}_pross_wt_{{t}}.sc", sample=SAMPLES, t=PROSS_TEMPS),
        #analysis
        expand(f"{WORKDIR}/{{sample}}_WT.fasta", sample=SAMPLES),
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}.fasta", sample=SAMPLES, variant=analysis_variants),
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.png", sample=SAMPLES, variant=analysis_variants),
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.csv", sample=SAMPLES, variant=analysis_variants),
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{m}}.txt", sample=SAMPLES, variant=analysis_variants, m=MUTATIONS),
        #design
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{mode}}/{{m}}.sc", m=MUTATIONS, sample=SAMPLES, mode=modes, variant=analysis_variants),
        #plotting
        expand(f"{WORKDIR}/{{variant}}/{{sample}}_energydifference_{{variant}}.png", sample=SAMPLES, variant=analysis_variants)
        #add proliNNator and disulfiNNate
        #expand(f"{WORKDIR}/{{sample}}_prolinnator.csv", sample=SAMPLES),
        #expand(f"{WORKDIR}/{{sample}}_disulfinnate.csv", sample=SAMPLES)

rule clean_pdb:
    input:
        pdb = f"{WORKDIR}/{{sample}}.pdb"
    output:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001.pdb"
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/clean_pdb_{{sample}}.log"
    shell:
        """
        mkdir -p {WORKDIR}/esm
        mkdir -p {WORKDIR}/pmpnn
        mkdir -p {WORKDIR}/indes
        mkdir -p {WORKDIR}/pross
        singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif score_jd2 \
        -renumber_pdb -ignore_unrecognized_res -s {input.pdb} \
        -out:pdb -out:suffix _clean -out:path:all {WORKDIR} > {log} 2>&1
        """

rule make_symmdef_file1:
    input:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}.symm",
        pdb = f"{WORKDIR}/{{sample}}_clean_0001_INPUT.pdb"
    resources:
        cpus=1
    log:
        f"{WORKDIR}/logs/make_symmdef_file1_{{sample}}.log"
    shell:
        """
        perl {INPUTDIR}/symmetry/make_symmdef_file.pl \
        -p {input.pdb} -a A -i B 2> {log} > {output.symm}
        """

# cp rather than mv: _clean_0001_INPUT.pdb is a declared output of make_symmdef_file1,
# and moving it leaves that rule's outputs incomplete on disk.
rule rename_file:
    input:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001_INPUT.pdb"
    output:
        pdbs = f"{WORKDIR}/{{sample}}_new.pdb"
    resources:
        cpus=1
    shell:
        """
        cp {input.pdb} {output.pdbs}
        """

rule relax:
    input:
        pdb = f"{WORKDIR}/{{sample}}_new.pdb",
        symm = f"{WORKDIR}/{{sample}}.symm"
    output:
        relaxed_pdb = f"{WORKDIR}/relax_{{sample}}_new_0001.pdb"
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/relax_{{sample}}.log"
    shell:
        """
        singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif relax \
        -s {input.pdb} \
        -constrain_relax_to_start_coords \
        -beta \
        -nstruct 1 \
        -multiple_processes_writing_to_one_directory \
        -out:prefix relax_ \
        -out:path:all {WORKDIR} \
        -symmetry_definition {input.symm} > {log} 2>&1
        """

rule make_symmdef_file2:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    resources:
        cpus=1
    log:
        f"{WORKDIR}/logs/make_symmdef_file2_{{sample}}.log"
    shell:
        """
        perl {INPUTDIR}/symmetry/make_symmdef_file.pl \
        -p {input.pdb} -a A -i B 2> {log} > {output.symm}
        """

rule run_esm:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        weights = f"{WORKDIR}/esm/{{sample}}_esm_probs.weights"
    params:
        protocol = f"{INPUTDIR}/esm/run_esm_and_save.xml"
    resources:
        mem_mb=16000,
        cpus=1
    log:
        f"{WORKDIR}/logs/run_esm_{{sample}}.log"
    shell:
        """
        singularity run -B {WORKDIR} -B {WORKDIR}/../esm2_t33_650M_UR50D:/usr/local/database/protocol_data/tensorflow_graphs/tensorflow_graph_repo_submodule/ESM/esm2_t33_650M_UR50D \
        {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
        -parser:protocol {params.protocol} \
        -parser:script_vars weights={output.weights} \
        -s {input.pdb} \
        -beta \
        -overwrite > {log} 2>&1
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
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/esm_sampling_{{sample}}_{{i}}.log"
    shell:
        """
        singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -parser:script_vars weights={input.weights} \
            -parser:script_vars resfile={params.resfile} \
            -out:pdb true \
            -out:path:all {WORKDIR}/esm/ \
            -out:prefix esm_ \
            -out:suffix _{wildcards.i} \
            -beta \
            -overwrite > {log} 2>&1
        """

rule run_pmpnn:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        weights = f"{WORKDIR}/pmpnn/{{sample}}_mpnn_probs.weights"
    params:
        protocol = f"{INPUTDIR}/pmpnn/run_mpnn_and_save.xml"
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/run_pmpnn_{{sample}}.log"
    shell:
        """
        singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
        -parser:protocol {params.protocol} \
        -parser:script_vars weights={output.weights} \
        -s {input.pdb} \
        -beta \
        -overwrite > {log} 2>&1
        """

rule pmpnn_sampling:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        weights = f"{WORKDIR}/pmpnn/{{sample}}_mpnn_probs.weights"
    output:
        pdb = f"{WORKDIR}/pmpnn/pmpnn_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb",
    params:
        protocol = f"{INPUTDIR}/pmpnn/sample_mutations.xml",
        resfile = f"{INPUTDIR}/pmpnn/resfile.resfile",
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/pmpnn_sampling_{{sample}}_{{i}}.log"
    shell:
        """
        singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -parser:script_vars weights={input.weights} \
            -parser:script_vars resfile={params.resfile} \
            -out:pdb true \
            -out:path:all {WORKDIR}/pmpnn/ \
            -out:prefix pmpnn_ \
            -out:suffix _{wildcards.i} \
            -beta \
            -overwrite > {log} 2>&1
        """

rule interface_design:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm"
    output:
        pdb = f"{WORKDIR}/indes/indes_relax_{{sample}}_new_0001_INPUT_{{i}}_0001.pdb"
    params:
        protocol = f"{INPUTDIR}/interface-design/sym_design.xml",
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/interface_design_{{sample}}_{{i}}.log"
    shell:
        """
        singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -beta \
            -out:pdb true \
            -out:path:all {WORKDIR}/indes/ \
            -out:prefix indes_ \
            -out:suffix _{wildcards.i} \
            -overwrite > {log} 2>&1
        """

rule get_fasta_from_pdbs:
    input:
        pdbs = lambda wildcards: expand(
            f"{WORKDIR}/{wildcards.variant}/{wildcards.variant}_relax_{wildcards.sample}_new_0001_INPUT_{{i}}_0001.pdb",
            i=ITERATIONS
        )
    output:
        fastafile = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}.fasta"
    params:
        script = f"{INPUTDIR}/get_fasta/get_multifasta_from_pdb_path.py"
    resources:
        cpus=1
    shell:
        """
        python {params.script} -p {input.pdbs} -c A -o {output.fastafile}
        """

# The reference sequence has to come from the same structure the designs are built
# from, otherwise cleaning/renumbering shifts every reported mutation position.
rule get_wt_fasta:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        fastafile = f"{WORKDIR}/{{sample}}_WT.fasta"
    params:
        script = f"{INPUTDIR}/get_fasta/get_multifasta_from_pdb_path.py"
    resources:
        cpus=1
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
    resources:
        cpus=1
    shell:
        """
        mkdir -p {WORKDIR}/pross
        bash {INPUTDIR}/pross/pssm/generate_pssm.file \
        {input.fastafile} {output.hhr} {output.a3m} {output.psi} {output.hhr_log} {output.pssm} {UNIREF_DB}
        bash {INPUTDIR}/pross/filter/make_cst.sh {input.pdb} > {output.cst}
        """

# One FilterScan pass writes the resfiles for *all* delta thresholds at once, so this
# rule has to declare all of them. Parameterising it by {t} instead would run the
# whole per-residue scan once per temperature, with every job writing the same files.
rule filterscan:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        cst = f"{WORKDIR}/pross/{{sample}}.cst",
        pssm = f"{WORKDIR}/pross/{{sample}}.pssm",
        fasta = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        resfiles = expand(
            WORKDIR + "/pross/{sample}_resfiles_pross/designable_aa_resfile.{t}",
            t=PROSS_TEMPS,
            allow_missing=True
        )
    params:
        protocol = f"{INPUTDIR}/pross/filter/filterscan.xml",
        path = f"{WORKDIR}/pross/{{sample}}_resfiles_pross/designable_aa_resfile",
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/filterscan_{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {params.path})
        : > {log}
        nres=$(grep -v '^>' {input.fasta} | tr -d '\n' | wc -c)
        for res in $(seq 1 $nres); do
            singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
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
                -overwrite >> {log} 2>&1
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
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/pross_design_{{sample}}_{{t}}.log"
    shell:
        """
        singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -parser:script_vars pdb_reference={input.pdb} \
            -parser:script_vars cst_full_path={input.cst} \
            -parser:script_vars cst_value=0.4 \
            -parser:script_vars pssm_full_path={input.pssm} \
            -parser:script_vars in_resfile={input.resfile} \
            -overwrite \
            -ignore_unrecognized_res \
            -use_input_sc \
            -use_occurrence_data \
            -out:file:scorefile {output.sc} \
            -out:path:all {WORKDIR}/pross \
            -out:prefix pross_design_{wildcards.t} \
            -beta > {log} 2>&1
        """

rule pross_design_wt:
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
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/pross_design_wt_{{sample}}_{{t}}.log"
    shell:
        """
        singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -parser:script_vars pdb_reference={input.pdb} \
            -parser:script_vars cst_full_path={input.cst} \
            -parser:script_vars cst_value=0.4 \
            -parser:script_vars pssm_full_path={input.pssm} \
            -parser:script_vars in_resfile={input.resfile} \
            -overwrite \
            -ignore_unrecognized_res \
            -use_input_sc \
            -use_occurrence_data \
            -out:file:scorefile {output.sc} \
            -out:prefix pross_wt_{wildcards.t} \
            -out:path:all {WORKDIR}/pross \
            -beta > {log} 2>&1
        """

rule plot_frequencies:
    input:
        fastafile = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}.fasta",
        wtfile = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        figure = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.png",
        csv = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.csv"
    params:
        script = f"{INPUTDIR}/validate/plot_frequencies.py",
        mutations = len(MUTATIONS)
    resources:
        cpus=1
    shell:
        """
        python {params.script} -i {input.fastafile} -r {input.wtfile} -m {params.mutations} -o {output.figure}
        """

rule get_mutation_list:
    input:
        csv = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}_frequency.csv"
    output:
        out = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{m}}.txt"
    params:
        script = f"{INPUTDIR}/validate/design-mutations.py"
    resources:
        cpus=1
    shell:
        """
        mkdir -p $(dirname {output.out})
        python {params.script} -i {input.csv} -o {output.out}
        """

rule run_design_or_control:
    input:
        txt = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{m}}.txt",
        symfile = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
    output:
        sc =  f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{mode}}/{{m}}.sc"
    params:
        xml = f"{INPUTDIR}/validate/design.v02.xml",
        outdir = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{mode}}",
        # Every job reads the same input structure, so Rosetta names all nstruct decoys
        # relax_{sample}_new_0001_INPUT_00XX.pdb. Without a private -out:path:pdb the
        # concurrent jobs of all mutations overwrite each other's structures.
        pdbdir = f"{WORKDIR}/{{variant}}/{{sample}}_{{variant}}/{{mode}}/{{m}}_decoys"
    resources:
        mem_mb=4000,
        cpus=1
    log:
        f"{WORKDIR}/logs/run_design_or_control_{{variant}}_{{sample}}_{{mode}}_{{m}}.log"
    shell:
        """
        MUTATION_LINE=$(cat {input.txt})

        mutpos=$(echo $MUTATION_LINE | cut -d'_' -f1)
        mutaa=$(echo $MUTATION_LINE | cut -d'_' -f2)

        mkdir -p {params.outdir} {params.pdbdir}
        singularity run -B {WORKDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
            -parser:protocol {params.xml} \
            -parser:script_vars mutpos=$mutpos mut_aa=$mutaa protocol={wildcards.mode} symfile={input.symfile} \
            -in:file:s {input.pdb} \
            -corrections:beta_nov16 \
            -out:file:scorefile {output.sc} \
            -out:path:pdb {params.pdbdir} \
            -overwrite \
            -nstruct 20 \
            -beta > {log} 2>&1
        """

rule plot_energy:
    input:
        control = lambda wildcards: expand(f"{WORKDIR}/{wildcards.variant}/{wildcards.sample}_{wildcards.variant}/control/{{m}}.sc", m=MUTATIONS),
        design = lambda wildcards: expand(f"{WORKDIR}/{wildcards.variant}/{wildcards.sample}_{wildcards.variant}/design/{{m}}.sc", m=MUTATIONS),
    output:
        image = f"{WORKDIR}/{{variant}}/{{sample}}_energydifference_{{variant}}.png"
    params:
        script = f"{INPUTDIR}/validate/plot_energies.py"
    resources:
        cpus=1
    shell:
        """
        python {params.script} \
            -i1 $(dirname {input.design[0]}) \
            -i2 $(dirname {input.control[0]}) \
            -o {output.image}
        """
