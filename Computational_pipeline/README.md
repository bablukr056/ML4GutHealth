# Metagenomic Data Download and Preprocessing

This workflow has two steps: first, raw sequencing reads are downloaded from SRA; second, each downloaded sample is processed through the preprocessing pipeline (quality trimming, host-read removal, taxonomic classification, and abundance estimation).

## Step 1 — Download raw reads (`sra_download.py`)

`sra_download.py` downloads a single SRA run and converts it to gzipped FASTQ files. Before running, activate the conda environment containing `prefetch` and `fastq-dump` (SRA Toolkit).

```bash
python sra_download.py -i <SAMPLE_ID> -o <OUTPUT_FOLDER>
```

**Example:**
```bash
python sra_download.py -i ERR11004600 -o /lustre/home/babluuniba2022/sra_downloads # run this on terminal
```

This creates a subfolder named `<SAMPLE_ID>` inside `<OUTPUT_FOLDER>`, containing the downloaded `.sra` file and the resulting FASTQ files (`<SAMPLE_ID>_pass_1.fastq.gz`, and `<SAMPLE_ID>_pass_2.fastq.gz` for paired-end data). A log file (`output.log`) records the download status for each sample processed.

## Step 2 — Run the preprocessing pipeline (`computational_pipeline.sh`)

Once the FASTQ files for a sample have been downloaded, `computational_pipeline.sh` processes that sample through four steps — quality trimming (Fastp), host-read removal (Bowtie2, GRCh38), taxonomic classification (Kraken2), and species-level abundance estimation (Bracken).

Before running, activate the conda environment containing these four tools, and confirm that `OUTPUT_FOLDER`, `KRAKEN2_DB`, and `BOWTIE2_DB` inside the script point to your own paths. The script expects raw FASTQ files named `<SAMPLE_ID>_pass_1.fastq.gz` (and `_pass_2.fastq.gz` for paired-end data) inside `<BIOPROJECT_ID>/<SAMPLE_FOLDER>/` — the same folder produced by Step 1.

The script takes three arguments — the BioProject accession, the subfolder containing that sample's FASTQ files, and the sample (run) accession:

```bash
chmod +x computational_pipeline.sh
./computational_pipeline.sh <BIOPROJECT_ID> <SAMPLE_FOLDER> <SAMPLE_ID>
```

**Example:**
```bash
./computational_pipeline.sh PRJEB36300 raw_fastq ERR11004600
```

For the full dataset, this pipeline was not run manually per sample. It was submitted as a batch of independent HTCondor jobs (2 CPUs, 64 GB RAM per job), one job per sample, with `computational_pipeline.sh` set as the job executable and the three per-sample arguments (`BIOPROJECT_ID`, `SAMPLE_FOLDER`, `SAMPLE_ID`) supplied through a job submission file listing one set of arguments per line. All samples in the dataset were processed this way as parallel Condor jobs. Errors at any step are logged to `error_log.txt` in the working directory.
