# SCRIPT Microbiology Whole-Genome Sequencing Pipeline

## Installation  

If you don't already have Miniconda or Anaconda installed, do so:
[https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

Or if you're on Quest, load the anaconda module:  
```
module load anaconda
```

Set up the conda environment containing the required software and activate it:

```
conda env create --file environment_simple.yml   
conda activate script_micro_wgs
```

Download the directory and files from github:

```
git clone https://github.com/egonozer/script_micro_wgs
cd script_micro_wgs
gunzip ref/PA_subelements.txt
```

Download the support software  
  
  * Be aware that the minikraken database is 8GB in size and could take a while to download. If you already downloaded it previously you can comment those lines in the *install_dependencies.sh* 
file.   
  
```
bash install_dependencies.sh
```

Set up the config.yaml configuration file:

1. Change the `## Kraken database` path to point to the minikraken database directory

If running on Quest, also change the quest/config.yaml file:

1. Change the `default-resources` to the account, partition, and email address you'd like to use (please don't use my email address!)


## Analysis

**1.** Activate the conda environment

```
conda activate script_micro_wgs 
```

**2.** Create the sample sheet:

* Either modify the example *samples.tsv* file or create your own tab-separated-file
* First line of the file must have the headers `sample read1 read2` separated by tabs
* Fill in the rest of the file, one line per sample, with the sample id in the first column and full paths to first and second reads in the second and third columns 

**3.** Modify the *config.yaml* file in the *script_micro_wgs* directory

* If your *samples.tsv* is not in the *script_micro_wgs* directory, then give the full path to the *samples.tsv* file on the second line. For example:  
   `samples: "/path/to/samples.tsv"`  
* Set the path to the output directory. **WARNING:** This directory must already exist. Example:  
    `outdir: "/path/to/wgs_output"`

**4.** Run the analysis:

* If you are not using Quest (set `--cores` to the maximum number of cores on your system you would like to use):  
`snakemake --cores 12`
* If you are using Quest:  
`snakemake --profile quest`  
If running many analyses, you may want to run snakemakke using `screen` or `nohup` in case you accidentally disconnect or you don't want to stay connected to wait for the results:  
`nohup snakemake --profile quest &`

## Results 
Results will be locaed in a directory titled **results** within your designated output directory. This directory will contain several other directories:  

| Directory | Description |  
|:----------:|------------|
| *agent* | Accessory genome sequences |
|*checkpoints*| Checkpoint files generated by the pipeline. Can be ignored or deleted.|
|*clustage*| Clustered accessory genome sequences. Clustering is against reference accessory element set given in config.yaml|
|*fastqc*| Sequencing read quality control|
|*kraken*| Genus and species assignment |
|*mlst* | MLST allele and sequence type assignment |
|*prokka* | Annotated genome sequence | 
|*snippy* | Reference alignmenet and variant calls | 
|*spades* | Genome assembly |
|*spades_filtered* | Quality filtered genome assembly |
|*trimmed_reads* | Sequencing reads after quality trimming | 

