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

## Set the output directory from teh config file
outdir = config["outdir"]

## These are the steps that do not need to be run as cluster jobs
localrules: all, clustage_input

rule all:
    input:
        expand("{outdir}/results/checkpoints/{sample}.fastqc", outdir=config["outdir"], sample=samples['sample']),
        expand("{outdir}/results/checkpoints/{sample}.kraken", outdir=config["outdir"], sample=samples['sample']),
        expand("{outdir}/results/checkpoints/{sample}.snippy", outdir=config["outdir"], sample=samples['sample']),
        expand("{outdir}/results/checkpoints/{sample}.mlst", outdir=config["outdir"], sample=samples['sample']),
        expand("{outdir}/results/checkpoints/all.clustage", outdir=config["outdir"]),
        expand("{outdir}/results/summary/snippy.processed_nt_variants.tsv", outdir=config["outdir"]),
        expand("{outdir}/results/summary/snippy.processed_aa_variants.tsv", outdir=config["outdir"]),
        expand("{outdir}/results/summary/PA_groups.txt", outdir=config["outdir"])

rule fastqc:
    input:
        read1=get_fq1,
        read2=get_fq2
    output:
        done="{outdir}/results/checkpoints/{sample}.fastqc"
    log: "{outdir}/results/fastqc/{sample}.log.txt"
    threads: 2
    conda: "envs/fastqc.yml"
    shell:
        "mkdir -p {outdir}/results/fastqc; "
        "fastqc --threads {threads} --outdir {outdir}/results/fastqc --extract {input.read1} {input.read2} 2>&1 > {log}; "
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
    conda: "envs/kraken.yml"
    shell:
        "export KRAKEN_DEFAULT_DB='{params.db}'; "
        "kraken --preload --threads {threads} --fastq-input --gzip-compressed "
        "--paired {input.read1} {input.read2} 2> {log} | kraken-report - 2>>{log} > {output.report}; "
        "perl scripts/kraken-topspecies.pl {output.report} {params.name} > {output.topspecies} 2>>{log}; "
        "touch {output.done}"

rule trimmomatic:
    input:
        read1=get_fq1,
        read2=get_fq2
    output:
        paired1="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_1.fastq.gz",
        paired2="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_2.fastq.gz",
        unpaired1="{outdir}/results/trimmed_reads/{sample}_trimmed_unpaired_1.fastq.gz",
        unpaired2="{outdir}/results/trimmed_reads/{sample}_trimmed_unpaired_2.fastq.gz",
        done="{outdir}/results/checkpoints/{sample}.trimmomatic"
    log: "{outdir}/results/trimmed_reads/{sample}.log.txt"
    threads: 12
    conda: "envs/trimmomatic.yml"
    shell:
        "trimmomatic PE -phred33 -threads {threads} {input.read1} {input.read2} "
        "{output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} "
        "ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 "
        "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 "
        "2>&1 > {log}; "
        "touch {output.done}"

rule spades:
    input:
        #trimmomatic_done="{outdir}/results/checkpoints/{sample}.trimmomatic",
        read1="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_1.fastq.gz",
        read2="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_2.fastq.gz"
    output:
        contigs="{outdir}/results/spades/{sample}/contigs.fasta",
        done="{outdir}/results/checkpoints/{sample}.spades"
    params:
        name="{outdir}/results/spades/{sample}"
    log: "{outdir}/results/spades/{sample}.log.txt"
    threads: 12
    conda: "envs/spades.yml"
    shell:
        "spades.py -o {params.name} --careful -1 {input.read1} -2 {input.read2} "
        "-t {threads} --cov-cutoff auto 2>&1 > {log}; "
        "touch {output.done}"

rule qc_filter:
    input:
        #spades_done="{outdir}/results/checkpoints/{sample}.spades",
        assembly="{outdir}/results/spades/{sample}/contigs.fasta",
        read1="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_1.fastq.gz",
        read2="{outdir}/results/trimmed_reads/{sample}_trimmed_paired_2.fastq.gz",
    output:
        filtered="{outdir}/results/spades_filtered/{sample}.filtered_sequences.fasta",
        stats="{outdir}/results/spades_filtered/{sample}.qc_filter_stats.txt",
        done="{outdir}/results/checkpoints/{sample}.qc_filter"
    params:
        name="{outdir}/results/spades_filtered/{sample}",
        phix=config["ref"]["phix"]
    log: "{outdir}/results/spades_filtered/{sample}.log.txt"
    threads: 12
    conda: "envs/qc_filter.yml"
    shell:
        "perl scripts/assembly_qc_filter.pl -c {input.assembly} "
        "-1 {input.read1} -2 {input.read2} "
        "-m 2 -l 200 -t {threads} -p {params.phix} -x 98 -o {params.name} "
        "> {output.stats} 2> {log}; "
        "touch {output.done}"

# rule vecscreen:
#     input:
#         qc_done="results/checkpoints/{sample}.qc_filter",
#         assembly="results/spades_filtered/{sample}.filtered_sequences.fasta"
#     output:
#         decon="results/vecscreen/{sample}.vecscreen.fasta",
#         log="results/vecscreen/{sample}.vecscreen_log.txt",
#         done="results/checkpoints/{sample}.qc_filter"
#     shell:
#         "perl scripts/vecscreen.pl ...; "
#         "touch {output.done}"

rule pa_group:
    input: 
        qcfilt_done=expand("{outdir}/results/checkpoints/{sample}.qc_filter", outdir=config['outdir'], sample=samples['sample'])
    output:
        group="{outdir}/results/summary/PA_groups.txt"
    params:
        list=expand("{outdir}/results/spades_filtered/{sample}.filtered_sequences.fasta", outdir=config['outdir'], sample=samples['sample'])
    threads: 2
    shell:
        'perl scripts/PA_group_profiler.pl {params.list} | '
        'perl -pe "s/.filtered_sequences//" > {output.group} '

rule mlst:
    input:
        #qc_done="{outdir}/results/checkpoints/{sample}.qc_filter",
        assembly="{outdir}/results/spades_filtered/{sample}.filtered_sequences.fasta"
    output:
        results="{outdir}/results/mlst/{sample}_mlst.txt",
        done="{outdir}/results/checkpoints/{sample}.mlst"
    params:
        profiles=config["mlst"]["profiles"],
        sequences=config["mlst"]["sequences"]
    log: "{outdir}/results/mlst/{sample}.log.txt"
    threads: 2
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
        assembly="{outdir}/results/spades_filtered/{sample}.filtered_sequences.fasta"
    output:
        gbk="{outdir}/results/prokka/{sample}/{sample}.gbk",
        done="{outdir}/results/checkpoints/{sample}.prokka"
    params:
        name="{sample}",
        ref=config["ref"]["gbk"]
    log: "{outdir}/results/prokka/{sample}.log.txt"
    threads: 12
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
        time="24:00:00"
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
        snps="{outdir}/results/snippy/{sample}/snps.tab",
        core_aa="{outdir}/results/snippy/{sample}/snps.processed_aa_variants.tsv",
        core_nt="{outdir}/results/snippy/{sample}/snps.processed_nt_variants.tsv"
    params:
        dir="{outdir}/results/snippy/{sample}",
        ref=config["ref"]["gbk"],
        bed=config["core"]["bed"]
    log: "{outdir}/results/snippy/{sample}.log.txt"
    threads: 12
    conda: "envs/snippy.yml"
    shell:
        "snippy --cpus {threads} --outdir {params.dir} --ref {params.ref} "
        "--R1 {input.read1} --R2 {input.read2} --force "
        "2>&1 > {log}; "
        "perl scripts/snippy_process.pl -s {output.snps} -b {params.bed} 2>>{log}; "
        "touch {output.done}"

rule snippy_compile:
    input:
        snippy_done=expand("{outdir}/results/checkpoints/{sample}.snippy", outdir=config['outdir'], sample=samples['sample'])
    output:
        nt="{outdir}/results/summary/snippy.processed_nt_variants.tsv",
        aa="{outdir}/results/summary/snippy.processed_aa_variants.tsv"
    threads: 12
    shell:
        'perl scripts/snippy_compile.pl -o {outdir}/results/summary/snippy {outdir}/results/snippy'


