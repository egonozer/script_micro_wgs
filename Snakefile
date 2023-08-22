import pandas as pd

configfile: "config.yaml"

## Read sample IDs and read files from a TSV file named "samples.tsv".
## This file must contain the headers "sample", "read1", and "read2".
## The paths to the read files can be absolute or relative.
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

def get_fq1(wildcards):
    return samples.loc[(wildcards.sample), ["read1"]].dropna()
def get_fq2(wildcards):
    return samples.loc[(wildcards.sample), ["read2"]].dropna()

## Set the output directory from the config file
outdir = config["outdir"]

## These are the steps that do not need to be run as cluster jobs
localrules: all, clustage_input, mqc

rule all:
    input:
        fastp_done=expand("{outdir}/results/checkpoints/{sample}.fastp", outdir=config["outdir"], sample=samples['sample']),
        fastqc_done=expand("{outdir}/results/checkpoints/{sample}.fastqc", outdir=config["outdir"], sample=samples['sample']),
        kraken_done=expand("{outdir}/results/checkpoints/{sample}.kraken", outdir=config["outdir"], sample=samples['sample']),
        snippy_done=expand("{outdir}/results/checkpoints/{sample}.snippy", outdir=config["outdir"], sample=samples['sample']),
        spades_done=expand("{outdir}/results/checkpoints/{sample}.spades", outdir=config["outdir"], sample=samples['sample']),
        fcsadapter_done=expand("{outdir}/results/checkpoints/{sample}.fcsadaptor", outdir=config["outdir"], sample=samples['sample']),
        mlst_done=expand("{outdir}/results/checkpoints/{sample}.mlst", outdir=config["outdir"], sample=samples['sample']),
        prokka_done=expand("{outdir}/results/checkpoints/{sample}.prokka", outdir=config["outdir"], sample=samples['sample']),
        clustage_done=expand("{outdir}/results/checkpoints/all.clustage", outdir=config["outdir"]),
        pagroup_done=expand("{outdir}/results/checkpoints/all.pagroup", outdir=config["outdir"]),
        snippycompile_done=expand("{outdir}/results/checkpoints/all.snippycompile", outdir=config["outdir"]),
        quast_done=expand("{outdir}/results/checkpoints/all.quast", outdir=config["outdir"]),
        multiqc_done=expand("{outdir}/results/checkpoints/all.multiqc", outdir=config["outdir"])

rule fastqc:
    input:
        read1=get_fq1,
        read2=get_fq2
    output:
        done="{outdir}/results/checkpoints/{sample}.fastqc"
    params:
        fqc_out="{outdir}/results/fastqc/{sample}"
    log: "{outdir}/results/fastqc/{sample}.log.txt"
    threads: 2
    resources: 
        jobname=lambda wildcards: f"smk-fastqc-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/fastqc/{wildcards.sample}-%j.out"
    conda: "envs/fastqc.yml"
    shell:
        "mkdir -p {outdir}/results/fastqc; "
        "mkdir -p {params.fqc_out}; "
        "fastqc --threads {threads} --outdir {params.fqc_out} --extract {input.read1} {input.read2} 2>&1 > {log}; "
        "touch {output.done}"

rule kraken:
    input:
        read1=get_fq1,
        read2=get_fq2
    output:
        report="{outdir}/results/kraken/{sample}.kraken_report.txt",
        topspecies="{outdir}/results/kraken/{sample}.kraken_topspecies.txt",
        done="{outdir}/results/checkpoints/{sample}.kraken"
    params:
        name="{sample}",
        db=config["krakendb"]
    log: "{outdir}/results/kraken/{sample}.log.txt"
    threads: 12
    resources: 
        jobname=lambda wildcards: f"smk-kraken-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/kraken/{wildcards.sample}-%j.out"
    conda: "envs/kraken.yml"
    shell:
        "export KRAKEN_DEFAULT_DB='{params.db}'; "
        "kraken --preload --threads {threads} --fastq-input --gzip-compressed "
        "--paired {input.read1} {input.read2} 2> {log} | kraken-report - 2>>{log} > {output.report}; "
        "perl scripts/kraken-topspecies.pl {output.report} {params.name} > {output.topspecies} 2>>{log}; "
        "touch {output.done}"

rule fastp:
    input:
        read1=get_fq1,
        read2=get_fq2
    output:
        paired1="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_1.fastq.gz",
        paired2="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_2.fastq.gz",
        unpaired1="{outdir}/results/trimmed_reads/{sample}_trimmed_unpaired_1.fastq.gz",
        unpaired2="{outdir}/results/trimmed_reads/{sample}_trimmed_unpaired_2.fastq.gz",
        loghtml="{outdir}/results/trimmed_reads/{sample}_fastp.html",
        logjson="{outdir}/results/trimmed_reads/{sample}_fastp.json",
        done="{outdir}/results/checkpoints/{sample}.fastp"
    log: "{outdir}/results/trimmed_reads/{sample}.log.txt"
    threads: 12
    resources:
        jobname=lambda wildcards: f"smk-fastp-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/fastp/{wildcards.sample}-%j.out"
    conda: "envs/fastp.yml"
    shell:
        "fastp "
        "--in1 {input.read1} --in2 {input.read2} "
        "--out1 {output.paired1} --unpaired1 {output.unpaired1} "
        "--out2 {output.paired2} --unpaired2 {output.unpaired2} "
        "-h {output.loghtml} -j {output.logjson} -w {threads} "
        "2>&1 > {log}; "
        "touch {output.done}"

rule spades:
    input:
        fastp_done="{outdir}/results/checkpoints/{sample}.fastp",
        read1="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_1.fastq.gz",
        read2="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_2.fastq.gz"
    output:
        contigs="{outdir}/results/spades/{sample}/contigs.fasta",
        done="{outdir}/results/checkpoints/{sample}.spades"
    params:
        name="{outdir}/results/spades/{sample}"
    log: "{outdir}/results/spades/{sample}.log.txt"
    threads: 12
    resources: 
        jobname=lambda wildcards: f"smk-spades-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/spades/{wildcards.sample}-%j.out"
    conda: "envs/spades.yml"
    shell:
        "spades.py -o {params.name} --careful -1 {input.read1} -2 {input.read2} "
        "-t {threads} --cov-cutoff auto 2>&1 > {log}; "
        "touch {output.done}"

# Remove contamination
rule fcsadaptor:
    input:
        spades_done="{outdir}/results/checkpoints/{sample}.spades",
        assembly="{outdir}/results/spades/{sample}/contigs.fasta"
    output:
        report="{outdir}/results/fcsadaptor/{sample}_fcsadaptor/fcs_adaptor_report.txt",
        done="{outdir}/results/checkpoints/{sample}.fcsadaptor"
    params:
        folder="{outdir}/results/fcsadaptor/{sample}_fcsadaptor",
        sif=config["fcssif"]
    log: "{outdir}/results/fcsadaptor/{sample}.fcsadaptor.log.txt"
    threads: 4
    resources: 
        jobname=lambda wildcards: f"smk-fcsadaptor-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/fcsadaptor/{wildcards.sample}-%j.out"
    envmodules: "singularity/3.8.1"
    shell:
        "scripts/run_fcsadaptor.sh "
        "--fasta-input {input.assembly} --output-dir {params.folder} "
        "--prok --container-engine singularity "
        "--image {params.sif} 2>&1 > {log} ; "
        "touch {output.done} "

rule fcsadaptor_filter:
    input:
        fcsadaptor_done="{outdir}/results/checkpoints/{sample}.fcsadaptor",
        report="{outdir}/results/fcsadaptor/{sample}_fcsadaptor/fcs_adaptor_report.txt",
        assembly="{outdir}/results/spades/{sample}/contigs.fasta"
    output:
        fixed="{outdir}/results/fcsadaptor/{sample}.fcsadaptor_fixed.fasta",
        done="{outdir}/results/checkpoints/{sample}.fcsadaptor_filter"
    log: "{outdir}/results/fcsadaptor/{sample}.filter.log.txt"
    threads: 4
    resources: 
        jobname=lambda wildcards: f"smk-fcs_filter-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/fcsadaptor_filter/{wildcards.sample}-%j.out"
    conda: "envs/perl.yml"
    shell:
        "perl scripts/ncbi_fcsadaptor_contamination_filter.pl "
        "-f {input.assembly} -c {input.report} "
        "> {output.fixed} 2> {log} ; "
        "module purge ;"
        "touch {output.done} "

rule qc_filter:
    input:
        fcsadaptor_done="{outdir}/results/checkpoints/{sample}.fcsadaptor_filter",
        assembly="{outdir}/results/fcsadaptor/{sample}.fcsadaptor_fixed.fasta",
        read1="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_1.fastq.gz",
        read2="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_2.fastq.gz",
    output:
        filtered="{outdir}/results/final_assemblies/{sample}.fasta",
        stats="{outdir}/results/spades_filtered/{sample}.qc_filter_stats.txt",
        done="{outdir}/results/checkpoints/{sample}.qc_filter"
    params:
        name="{outdir}/results/spades_filtered/{sample}",
        intermediate="{outdir}/results/spades_filtered/{sample}.filtered_sequences.fasta",
        phix=config["ref"]["phix"]
    log: "{outdir}/results/spades_filtered/{sample}.log.txt"
    threads: 4
    resources: 
        jobname=lambda wildcards: f"smk-qc_filter-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/qc_filter/{wildcards.sample}-%j.out"
    conda: "envs/qc_filter.yml"
    shell:
        "perl scripts/assembly_qc_filter.pl -c {input.assembly} "
        "-1 {input.read1} -2 {input.read2} "
        "-m 2 -l 200 -t {threads} -p {params.phix} -x 98 -o {params.name} "
        "> {output.stats} 2> {log}; "
        "mkdir -p {outdir}/results/final_assemblies; "
        "mv {params.intermediate} {output.filtered}; "
        "touch {output.done}"

rule quast:
    input:
        qcfilt_done=expand("{outdir}/results/checkpoints/{sample}.qc_filter", outdir=config['outdir'], sample=samples['sample'])
    output:
        done="{outdir}/results/checkpoints/all.quast",
        quast_report="{outdir}/results/final_assemblies/quast_output/report.txt"
    params:
        quast_in_dir="{outdir}/results/final_assemblies",
        quast_out_dir="{outdir}/results/final_assemblies/quast_output",
        ref=config["ref"]["fa"]
    log: "{outdir}/results/final_assemblies/quast.log.txt"
    threads: 12
    resources: 
        jobname="smk-quast",
        outlog="logs/quast/all-%j.out"
    conda: "envs/quast.yml"
    shell:
        "quast -o {params.quast_out_dir} -r {params.ref} -t {threads} {params.quast_in_dir}/*.fasta ; "
        "touch {output.done}"
    
rule pa_group:
    input: 
        qcfilt_done=expand("{outdir}/results/checkpoints/{sample}.qc_filter", outdir=config['outdir'], sample=samples['sample'])
    output:
        done="{outdir}/results/checkpoints/all.pagroup",
        group="{outdir}/results/summary/PA_groups.txt"
    params:
        list=expand("{outdir}/results/final_assemblies/{sample}.fasta", outdir=config['outdir'], sample=samples['sample'])
    threads: 2
    resources:
        jobname="smk-pa_group",
        outlog="logs/pa_group/all-%j.out"
    conda: "envs/perl.yml"
    shell:
        'perl scripts/PA_group_profiler.pl {params.list} > {output.group}; '
        'touch {output.done}'

rule mlst:
    input:
        #qc_done="{outdir}/results/checkpoints/{sample}.qc_filter",
        assembly="{outdir}/results/final_assemblies/{sample}.fasta"
    output:
        results="{outdir}/results/mlst/{sample}_mlst.txt",
        done="{outdir}/results/checkpoints/{sample}.mlst"
    params:
        profiles=config["mlst"]["profiles"],
        sequences=config["mlst"]["sequences"]
    log: "{outdir}/results/mlst/{sample}.log.txt"
    threads: 2
    resources:
        jobname=lambda wildcards: f"smk-mlst-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/mlst/{wildcards.sample}-%j.out"
    conda: "envs/mlst.yml"
    shell:
        "perl scripts/mlst_profiler/mlst_profiler.pl "
        "-b $CONDA_PREFIX/bin "
        "-p {params.profiles} "
        "-a {params.sequences} "
        "-f {input.assembly} > {output.results} 2> {log}; "
        "touch {output.done}"

rule prokka:
    input:
        #qc_done="{outdir}/results/checkpoints/{sample}.qc_filter",
        assembly="{outdir}/results/final_assemblies/{sample}.fasta"
    output:
        gbk="{outdir}/results/prokka/{sample}/{sample}.gbk",
        done="{outdir}/results/checkpoints/{sample}.prokka"
    params:
        name="{sample}",
        ref=config["ref"]["gbk"]
    log: "{outdir}/results/prokka/{sample}.log.txt"
    threads: 12
    resources: 
        jobname=lambda wildcards: f"smk-prokka-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/prokka/{wildcards.sample}-%j.out"
    conda: "envs/prokka.yml"
    shell:
        "prokka --outdir {outdir}/results/prokka/{params.name} "
        "--prefix {params.name} --locustag {params.name} "
        "--genus Pseudomonas --species aeruginosa "
        "--strain {params.name} --kingdom Bacteria "
        "--proteins {params.ref} --force "
        "--addgenes --mincontiglen 200 --cpus {threads} "
        "{input.assembly} 2>&1 > {log}; "
        "touch {output.done}"

rule agent:
    input:
        #prokka_done="{outdir}/results/checkpoints/{sample}.prokka",
        assembly="{outdir}/results/prokka/{sample}/{sample}.gbk",
    output:
        fasta="{outdir}/results/agent/{sample}.agent.accessory.fasta",
        annot="{outdir}/results/agent/{sample}.agent.accessory_loci.txt",
        done="{outdir}/results/checkpoints/{sample}.agent"
    params:
        name="{sample}",
        core=config["core"]["seq"]
    log: "{outdir}/results/agent/{sample}.log.txt"
    resources:
        jobname=lambda wildcards: f"smk-agent-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/agent/{wildcards.sample}-%j.out"
    conda: "envs/agent.yml"
    shell:
        "perl AGEnt/AGEnt.pl -q {input.assembly} -r {params.core} "
        "-n $CONDA_PREFIX/bin "
        "-o {outdir}/results/agent/{params.name} -p agent 2> {log}; "
        "touch {output.done}"

rule clustage_input:
    input:
        agent_done=expand("{outdir}/results/checkpoints/{sample}.agent", outdir=config['outdir'], sample=samples['sample'])
    output:
        files="{outdir}/results/clustage/clustage_seqs_in.txt",
        annot="{outdir}/results/clustage/clustage_annot_in.txt"
    params:
        list=expand("{sample}", sample=samples['sample'])
    resources:
        jobname="smk-clst_input",
        outlog="logs/clustage_input/all-%j.out"
    conda: "envs/perl.yml"
    shell:
        'perl scripts/clustage_input.pl -d {outdir} {params.list} > {output.files} 2> {output.annot}'
        
        # 'touch {output.files}; touch {output.annot}; '
        # 'for i in results/agent/*.agent.accessory.fasta; '
        # 'do j=`echo $i | perl -pe "s/^(.*)\.agent.accessory.fasta/\1/"`; echo $j'
        # 'echo -e "$i\t$j" >> {output.files}; '
        # 'echo -e "results/agent/$j.agent.accessory_loci.txt\t$j" >> {output.annot}; '
        # 'done '

rule clustage:
    input:
        files="{outdir}/results/clustage/clustage_seqs_in.txt",
        annot="{outdir}/results/clustage/clustage_annot_in.txt",
    output:
        done="{outdir}/results/checkpoints/all.clustage"
    params:
        ages=config["clustage"]["ages"],
        subelems=config["clustage"]["subelems"]
    log: "{outdir}/results/clustage/log.txt"
    threads: 12
    resources:
        time="24:00:00",
        jobname="smk-clustage",
        outlog="logs/clustage/all-%j.out"
    conda: "envs/clustage.yml"
    shell:
        "perl ClustAGE/ClustAGE.pl -f {input.files} --annot {input.annot} "
        "--blastpath $CONDA_PREFIX/bin -o {outdir}/results/clustage/clustage "
        "--age {params.ages} --subelem {params.subelems} --add_ages "
        "-p $CONDA_PREFIX/bin/gnuplot 2> {log}; "
        "touch {output.done}"

rule snippy:
    input:
        #fastqc_done="{outdir}/results/checkpoints/{sample}.fastqc",
        read1=get_fq1,
        read2=get_fq2
    output:
        done="{outdir}/results/checkpoints/{sample}.snippy",
        snps="{outdir}/results/snippy/{sample}/{sample}.tab",
        core_aa="{outdir}/results/snippy/{sample}/snps.processed_aa_variants.tsv",
        core_nt="{outdir}/results/snippy/{sample}/snps.processed_nt_variants.tsv"
    params:
        prefix="{sample}",
        dir="{outdir}/results/snippy/{sample}",
        ref=config["ref"]["gbk"],
        bed=config["core"]["bed"]
    log: "{outdir}/results/snippy/{sample}.log.txt"
    threads: 12
    resources: 
        jobname=lambda wildcards: f"smk-snippy-{wildcards.sample}",
        outlog=lambda wildcards: f"logs/snippy/{wildcards.sample}-%j.out"
    conda: "envs/snippy.yml"
    shell:
        "snippy --cpus {threads} --outdir {params.dir} --ref {params.ref} "
        "--R1 {input.read1} --R2 {input.read2} --force "
        "--prefix {params.prefix} "
        "2>&1 > {log}; "
        "perl scripts/snippy_process.pl -s {output.snps} -b {params.bed} 2>>{log}; "
        "touch {output.done}"

rule snippy_compile:
    input:
        snippy_done=expand("{outdir}/results/checkpoints/{sample}.snippy", outdir=config['outdir'], sample=samples['sample'])
    output:
        done="{outdir}/results/checkpoints/all.snippycompile",
        nt="{outdir}/results/summary/snippy.processed_nt_variants.tsv",
        aa="{outdir}/results/summary/snippy.processed_aa_variants.tsv"
    threads: 12
    resources:
        jobname="smk-snp_compile",
        outlog="logs/snippy_compile/all-%j.out"
    conda: "envs/perl.yml"
    shell:
        'perl scripts/snippy_compile.pl -o {outdir}/results/summary/snippy {outdir}/results/snippy; '
        'touch {output.done}'

rule mqc:
    input:
        fastp_done=expand("{outdir}/results/checkpoints/{sample}.fastp", outdir=config["outdir"], sample=samples['sample']),
        fastqc_done=expand("{outdir}/results/checkpoints/{sample}.fastqc", outdir=config["outdir"], sample=samples['sample']),
        kraken_done=expand("{outdir}/results/checkpoints/{sample}.kraken", outdir=config["outdir"], sample=samples['sample']),
        snippy_done=expand("{outdir}/results/checkpoints/{sample}.snippy", outdir=config["outdir"], sample=samples['sample']),
        spades_done=expand("{outdir}/results/checkpoints/{sample}.spades", outdir=config["outdir"], sample=samples['sample']),
        mlst_done=expand("{outdir}/results/checkpoints/{sample}.mlst", outdir=config["outdir"], sample=samples['sample']),
        prokka_done=expand("{outdir}/results/checkpoints/{sample}.prokka", outdir=config["outdir"], sample=samples['sample']),
        clustage_done=expand("{outdir}/results/checkpoints/all.clustage", outdir=config["outdir"]),
        pagroup_done=expand("{outdir}/results/checkpoints/all.pagroup", outdir=config["outdir"]),
        snippycompile_done=expand("{outdir}/results/checkpoints/all.snippycompile", outdir=config["outdir"]),
        quast_done=expand("{outdir}/results/checkpoints/all.quast", outdir=config["outdir"])
    output:
        done="{outdir}/results/checkpoints/all.multiqc",
        a_mqc="{outdir}/results/agent/agent_mqc.txt",
        report="{outdir}/results/multiqc_report.html"
    params:
        a_dir="{outdir}/results/agent",
        results_dir="{outdir}/results"
    conda: "envs/multiqc.yml"
    shell:
        "perl scripts/mqc_agent.pl {params.a_dir} > {output.a_mqc}; "
        "multiqc -f -c multiqc_config.yaml -o {params.results_dir} {params.results_dir}; "
        "touch {output.done}"
