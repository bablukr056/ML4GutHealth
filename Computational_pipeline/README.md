# Metagenomic Preprocessing Pipeline

`computational_pipeline.sh` processes one raw metagenomic sample at a time through four steps — quality trimming (Fastp), host-read removal (Bowtie2, GRCh38), taxonomic classification (Kraken2), and species-level abundance estimation (Bracken). Before running, activate the conda environment containing these four tools, and confirm that `OUTPUT_FOLDER`, `KRAKEN2_DB`, and `BOWTIE2_DB` inside the script point to your own paths. The script expects raw FASTQ files named `<SAMPLE>_pass_1.fastq.gz` (and `_pass_2.fastq.gz` for paired-end data) inside `<BIOPROJECT>/<SAMPLE_FOLDER>/`.

Run it with three arguments — BioProject accession, the subfolder containing that sample's FASTQ files, and the sample (run) accession:

```bash
chmod +x bioinfor_pipeline.sh
./computational_pipeline.sh <BIOPROJECT> <SAMPLE_FOLDER> <SAMPLE>
```

For example: `./bioinfor_pipeline.sh PRJEB36300 raw_fastq ERR11004600`. For a full dataset, this script is not run manually per sample — it is submitted once per sample as an independent HTCondor job (32 CPUs, 272 GB RAM per job), with the per-sample arguments supplied through a job submission file. Errors at any step are logged to `error_log.txt` in the working directory.
