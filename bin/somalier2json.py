#!/usr/bin/env python3

#!/usr/bin/env python3
import csv
import json
import argparse
from datetime import datetime, timezone

sex_codes = {   '1' : "male", 
                '2' : "female",
                '0' : "unknown"
            }

def parse_samples(sample_args):
    """
    provides a list of sample:run_id
    Returns:
    {
        "sample1": "seq1",
        "sample2": "seq2"
    }
    """
    samples = {}
    for part in sample_args:
        if ":" in part:
            sample, seqrun = part.split(":", 1)
        else:
            exit("Need format s1:run_idX")
        samples[sample] = seqrun
    return samples

def load_somalier_sample(file_path):
    """"
    From the aggregated results from the tumor and normal samples
    Get the sex information from the both somalier and the orignal
    pedigree sex
    """
    ped_rows = []
    with open(file_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ped_rows.append(row)
    return ped_rows

def check_sex(row):
    sex_ok = False
    pedigree_sex = row.get("original_pedigree_sex")
    if sex_codes[row.get("sex")] == row.get("original_pedigree_sex"):
        sex_ok = True
    else:
        sex_ok = False
    return sex_ok, pedigree_sex

def output_json(results, output):
    filename = output
    with open(filename, "w") as f:
        json.dump(results, f, indent=2)

def write_cdm_load(sample, cdmassay, output_file, samples_dict, results_dir):
    seqrun = samples_dict.get(sample)
    filename = f"{sample}.peddy2cdm"
    with open(filename, "w") as f:
        f.write(
            f"--sequencing-run {seqrun} --assay {cdmassay} --sample-id {sample} --peddy {results_dir}/{output_file}"
        )

def main(somalier_sample_file, sample_arg, cdmassay, results_dir):
    somalier_samples = load_somalier_sample(somalier_sample_file)
    samples_dict = parse_samples(sample_arg)
    for item in somalier_samples:
        stautus, sex_original = check_sex(item)
        current_time = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        sample = item.get("sample_id")
        results = {
            "sex_check": {
                "is_correct_sex": stautus,
                "pedigree_sex": sex_original
                },
                "analysis_date": current_time
            }
        
        output_file = f"{sample}.somalier.json"
        output_json(results, output_file)
        write_cdm_load(sample, cdmassay, output_file, samples_dict, results_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PED and SEX check files from somalier with WGS")
    parser.add_argument("--somalier", required=True, help="id.samples.tsv")
    parser.add_argument("--samples", required=True, action="append", help="sample:run_id (can be used multiple times)")
    parser.add_argument("--cdmassay", required=True, help="cdm assay of sample")
    parser.add_argument("--results_dir", required=True, help="cdm assay of sample")

    args = parser.parse_args()

    main(args.somalier, args.samples, args.cdmassay, args.results_dir)
