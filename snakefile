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
# Which profile drives the PROSS design: the hhblits/psiblast alignment ("msa"),
# or the ESM / ProteinMPNN probabilities converted to log-odds.
PSSM_SOURCES = config["pssm_sources"]

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
    p = r"msa|esm|pmpnn",
    r = r"\d+",
    sample = r"[^/]+",
    variant = r"[^/]+",
    mode = r"design|control",
    t = r"[^/]+"

# "<sample>_2.symm" matches make_symmdef_file1 as well, with sample="<sample>_2".
# Both rules can therefore claim that file, so state which one actually produces it.
ruleorder: make_symmdef_file2 > make_symmdef_file1

# Same trap for "<sample>_<source>.pssm", which generate_PSSM_and_constraints matches
# with sample="<sample>_<source>". The converted neural profiles own that name.
ruleorder: pssm_from_weights > generate_PSSM_and_constraints

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
        #pross (needs hhblits, BLAST+ and the UniRef30 database on the host)
        expand(f"{WORKDIR}/pross/{{sample}}.pssm", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}.cst", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}.hhr", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}.a3m", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}.psi", sample=SAMPLES),
        expand(f"{WORKDIR}/pross/{{sample}}_resfiles_{{p}}/designable_aa_resfile.{{t}}", sample=SAMPLES, t=PROSS_TEMPS, p=PSSM_SOURCES),
        expand(f"{WORKDIR}/pross/{{sample}}_pross_design_{{p}}_{{t}}.sc", sample=SAMPLES, t=PROSS_TEMPS, p=PSSM_SOURCES),
        expand(f"{WORKDIR}/pross/{{sample}}_pross_wt_{{p}}_{{t}}.sc", sample=SAMPLES, t=PROSS_TEMPS, p=PSSM_SOURCES),
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
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/clean_pdb_{{sample}}.log"
    shell:
        """
        mkdir -p {WORKDIR}/esm
        mkdir -p {WORKDIR}/pmpnn
        mkdir -p {WORKDIR}/indes
        mkdir -p {WORKDIR}/pross
        singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif score_jd2 \
        -renumber_pdb -ignore_unrecognized_res -s {input.pdb} \
        -out:pdb -out:suffix _clean -out:path:all {WORKDIR} > {log} 2>&1
        """

rule make_symmdef_file1:
    input:
        pdb = f"{WORKDIR}/{{sample}}_clean_0001.pdb"
    output:
        symm = f"{WORKDIR}/{{sample}}.symm",
        pdb = f"{WORKDIR}/{{sample}}_clean_0001_INPUT.pdb"
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
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/relax_{{sample}}.log"
    shell:
        """
        singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif relax \
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
        mem_mb=16000
    log:
        f"{WORKDIR}/logs/run_esm_{{sample}}.log"
    shell:
        """
        singularity run -B {WORKDIR} -B {INPUTDIR} -B {WORKDIR}/../esm2_t33_650M_UR50D:/usr/local/database/protocol_data/tensorflow_graphs/tensorflow_graph_repo_submodule/ESM/esm2_t33_650M_UR50D \
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
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/esm_sampling_{{sample}}_{{i}}.log"
    shell:
        """
        singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
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
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/run_pmpnn_{{sample}}.log"
    shell:
        """
        singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
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
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/pmpnn_sampling_{{sample}}_{{i}}.log"
    shell:
        """
        singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
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
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/interface_design_{{sample}}_{{i}}.log"
    shell:
        """
        singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
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
    shell:
        """
        python {params.script} -p {input.pdbs} -c A -o {output.fastafile}
        """

# The reference sequence has to come from the same structure the designs are built
# from, otherwise cleaning/renumbering shifts every reported mutation position.
# A checkpoint rather than a plain rule: filterscan fans out into one job per residue,
# and the residue count is only known once this sequence exists, which is after the DAG
# has been built. The checkpoint makes Snakemake re-evaluate the DAG at that point.
checkpoint get_wt_fasta:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb"
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
    # hhblits holds the UniRef30 prefilter in memory: ~12 GB for a 44-residue query,
    # more for longer ones. The 4 GB default would get the job OOM-killed.
    threads: 8
    resources:
        mem_mb=32000
    log:
        f"{WORKDIR}/logs/generate_PSSM_and_constraints_{{sample}}.log"
    shell:
        """
        mkdir -p {WORKDIR}/pross
        bash {INPUTDIR}/pross/pssm/generate_pssm.file \
        {input.fastafile} {output.hhr} {output.a3m} {output.psi} {output.hhr_log} {output.pssm} {UNIREF_DB} {threads} > {log} 2>&1
        bash {INPUTDIR}/pross/filter/make_cst.sh {input.pdb} > {output.cst} 2>> {log}
        """

def wt_sequence_length(sample):
    """Residue count of the scanned chain, read after get_wt_fasta has run."""
    fasta = checkpoints.get_wt_fasta.get(sample=sample).output.fastafile
    with open(fasta) as handle:
        seq = "".join(line.strip() for line in handle if not line.startswith(">"))
    return len(seq)


def filterscan_residue_resfiles(wildcards):
    return expand(
        WORKDIR + "/pross/{sample}_resfiles_{p}/res{r}/designable_aa_resfile.{t}",
        sample=wildcards.sample,
        p=wildcards.p,
        t=wildcards.t,
        r=range(1, wt_sequence_length(wildcards.sample) + 1),
    )


# "msa" is the classic PROSS profile built by hhblits/psiblast; "esm" and "pmpnn" are
# the same design protocol driven by the neural probabilities the sampling branches
# already produce, converted to log-odds.
def pssm_for_source(wildcards):
    if wildcards.p == "msa":
        return f"{WORKDIR}/pross/{wildcards.sample}.pssm"
    return f"{WORKDIR}/pross/{wildcards.sample}_{wildcards.p}.pssm"


# The two neural branches name their probability tables differently.
PROBS_FILE = {
    "esm": WORKDIR + "/esm/{sample}_esm_probs.weights",
    "pmpnn": WORKDIR + "/pmpnn/{sample}_mpnn_probs.weights",
}


rule pssm_from_weights:
    input:
        weights = lambda wc: PROBS_FILE[wc.p].format(sample=wc.sample),
        fasta = f"{WORKDIR}/{{sample}}_WT.fasta"
    output:
        pssm = f"{WORKDIR}/pross/{{sample}}_{{p}}.pssm"
    params:
        script = f"{INPUTDIR}/pross/pssm/weights_to_pssm.py"
    wildcard_constraints:
        p = "esm|pmpnn"
    log:
        f"{WORKDIR}/logs/pssm_from_weights_{{sample}}_{{p}}.log"
    shell:
        """
        python {params.script} -i {input.weights} -s {input.fasta} -o {output.pssm} > {log} 2>&1
        """


# One FilterScan call scans a single residue but writes it to the resfiles of *all*
# delta thresholds at once, so this rule has to declare all of them. The residues are
# independent, hence one job each; they cannot share an output file, because FilterScan
# appends and concurrent jobs would interleave their lines.
rule filterscan_residue:
    input:
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        cst = f"{WORKDIR}/pross/{{sample}}.cst",
        pssm = pssm_for_source
    output:
        resfiles = expand(
            WORKDIR + "/pross/{sample}_resfiles_{p}/res{r}/designable_aa_resfile.{t}",
            t=PROSS_TEMPS,
            allow_missing=True
        )
    params:
        protocol = f"{INPUTDIR}/pross/filter/filterscan.xml",
        path = f"{WORKDIR}/pross/{{sample}}_resfiles_{{p}}/res{{r}}/designable_aa_resfile",
        temps = " ".join(str(t) for t in PROSS_TEMPS),
    resources:
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/filterscan_{{sample}}_{{p}}_res{{r}}.log"
    shell:
        """
        mkdir -p $(dirname {params.path})
        # The image ships an MPI build of Rosetta. Snakemake runs rules inside an `srun`
        # job step, and a second MPI_Init in the same step makes PMIx abort and take the
        # step down with it. One Rosetta call per job avoids that, but Snakemake may put
        # several jobs in one step when they are grouped, so start Rosetta standalone.
        UNSET=$(env | grep -oE '^(PMIX|PMI|SLURM)_[A-Za-z0-9_]*' | sed 's/^/-u /')
        env $UNSET singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
            -parser:protocol {params.protocol} \
            -s {input.pdb} \
            -parser:script_vars sym={input.symm} \
            -parser:script_vars pdb_reference={input.pdb} \
            -parser:script_vars cst_full_path={input.cst} \
            -parser:script_vars cst_value=0.4 \
            -parser:script_vars pssm_full_path={input.pssm} \
            -parser:script_vars resfiles_path={params.path} \
            -parser:script_vars current_res={wildcards.r} \
            -out:path:all {WORKDIR}/pross \
            -beta \
            -overwrite > {log} 2>&1
        # A position whose profile admits nothing but the native residue has no mutation
        # to scan, so FilterScan writes no resfile at all. That is a real result, not a
        # failure -- the sharper ESM/ProteinMPNN profiles hit it where an alignment does
        # not -- so leave the empty resfiles behind for the merge to skip over.
        for t in {params.temps}; do
            [ -f "{params.path}.$t" ] || {{ echo nataa; echo start; }} > "{params.path}.$t"
        done
        """

# Stitch the per-residue resfiles back into the single file per threshold that
# pross_design expects. `input` arrives in residue order, so the body keeps the
# ascending numbering the serial scan produced.
rule merge_filterscan_resfiles:
    input:
        filterscan_residue_resfiles
    output:
        resfile = f"{WORKDIR}/pross/{{sample}}_resfiles_{{p}}/designable_aa_resfile.{{t}}"
    shell:
        """
        set -- {input}
        sed -n '1,/^start$/p' "$1" > {output.resfile}
        for f in {input}; do
            sed -n '/^start$/,$p' "$f" | tail -n +2 >> {output.resfile}
        done
        """

rule pross_design:
    input:
        resfile = f"{WORKDIR}/pross/{{sample}}_resfiles_{{p}}/designable_aa_resfile.{{t}}",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        cst = f"{WORKDIR}/pross/{{sample}}.cst",
        pssm = pssm_for_source
    output:
        sc = f"{WORKDIR}/pross/{{sample}}_pross_design_{{p}}_{{t}}.sc"
    params:
        protocol = f"{INPUTDIR}/pross/design/design.xml",
    resources:
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/pross_design_{{sample}}_{{p}}_{{t}}.log"
    shell:
        """
        singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
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
            -out:prefix pross_design_{wildcards.p}_{wildcards.t} \
            -beta > {log} 2>&1
        """

rule pross_design_wt:
    input:
        resfile = f"{WORKDIR}/pross/{{sample}}_resfiles_{{p}}/designable_aa_resfile.{{t}}",
        symm = f"{WORKDIR}/{{sample}}_2.symm",
        pdb = f"{WORKDIR}/relax_{{sample}}_new_0001_INPUT.pdb",
        cst = f"{WORKDIR}/pross/{{sample}}.cst",
        pssm = pssm_for_source
    output:
        sc = f"{WORKDIR}/pross/{{sample}}_pross_wt_{{p}}_{{t}}.sc"
    params:
        protocol = f"{INPUTDIR}/pross/design/design_WT.xml",
    resources:
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/pross_design_wt_{{sample}}_{{p}}_{{t}}.log"
    shell:
        """
        singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
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
            -out:prefix pross_wt_{wildcards.p}_{wildcards.t} \
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
        mem_mb=4000
    log:
        f"{WORKDIR}/logs/run_design_or_control_{{variant}}_{{sample}}_{{mode}}_{{m}}.log"
    shell:
        """
        MUTATION_LINE=$(cat {input.txt})

        mutpos=$(echo $MUTATION_LINE | cut -d'_' -f1)
        mutaa=$(echo $MUTATION_LINE | cut -d'_' -f2)

        mkdir -p {params.outdir} {params.pdbdir}
        singularity run -B {WORKDIR} -B {INPUTDIR} {ROSETTA_DIR}/rosetta_ml.sif rosetta_scripts \
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
    shell:
        """
        python {params.script} \
            -i1 $(dirname {input.design[0]}) \
            -i2 $(dirname {input.control[0]}) \
            -o {output.image}
        """
