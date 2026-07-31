#!/lustrehome/babluuniba2022/miniconda3/envs/sra/bin/python

# ==============================================================================
# SRA Data Download Wrapper (prefetch + fastq-dump)
#
# Downloads a single SRA run (prefetch) and converts it to gzipped FASTQ
# files (fastq-dump), logging the success/failure status of each step.
#
# Usage:
#   python sra_download.py -i <SRA_RUN_ID> -o <OUTPUT_FOLDER>
#
# Arguments:
#   -i, --sra_id       SRA run accession to download (e.g. SRR1234567). Required.
#   -o, --out_folder    Output directory where the .sra file and resulting
#                        FASTQ files will be saved. Optional; defaults to the
#                        current working directory.
#
# Example:
#   python sra_download.py -i SRR1234567 -o /lustre/home/babluuniba2022/sra_downloads
#
# A log file (output.log) is created inside the output folder, recording the
# combined success/failure status of the prefetch and fastq-dump steps for
# each SRA run processed.
#
# For batch downloading many samples, this script is called once per sample,
# typically as part of a job-list-driven HTCondor submission.
# ==============================================================================

import subprocess
import os
import argparse
import argcomplete


def be_parser():
    import sys
    parser = argparse.ArgumentParser(description="Wrapper to download SRA project files",
                                     prefix_chars="-")
    parser.add_argument("-o", "--out_folder", type=str,
                        help="Output directory for the downloaded .sra and FASTQ files",
                        action="store",
                        required=False, default=os.getcwd())
    parser.add_argument("-i", "--sra_id", type=str, help="SRA run ID to download (e.g. SRR1234567)",
                        action="store", required=True)
    argcomplete.autocomplete(parser)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def prefetch_fastq_dump_dwn(data_list):
    _output_dir, _output_file, sra_id = data_list
    print(f"Currently downloading: {sra_id}")
    prefetch = f"/lustrehome/babluuniba2022/miniconda3/envs/sra/bin/prefetch --output-directory {_output_dir} {sra_id}"
    print(f"The command used was: {prefetch}")
    return_code = subprocess.call(prefetch, shell=True)
    prefetch_status = "success" if return_code == 0 else "fail"

    print(f"Generating fastq for: {sra_id}")
    sra_path = os.path.join(_output_dir, sra_id, f"{sra_id}.sra")
    fastq_dir = os.path.join(_output_dir, sra_id)
    os.makedirs(fastq_dir, exist_ok=True)
    fastq_dump = f"/lustrehome/babluuniba2022/miniconda3/envs/sra/bin/fastq-dump --outdir {fastq_dir} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-files --clip {sra_path}"
    print(f"The command used was: {fastq_dump}")
    return_code = subprocess.call(fastq_dump, shell=True)
    fastq_status = "success" if return_code == 0 else "fail"

    overall_status = "success" if prefetch_status == "success" and fastq_status == "success" else "fail"

    with open(_output_file, "a") as _f:
        _f.write(f"{sra_id}\t{overall_status}\n")


if __name__ == "__main__":
    args = be_parser()
    sra_id, output_dir = args.sra_id, args.out_folder

    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, "output.log")

    prefetch_fastq_dump_dwn([output_dir, output_file, sra_id])
