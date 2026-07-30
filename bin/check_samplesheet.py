#!/usr/bin/env python3

# Released under the MIT license.
# See git repository (https://github.com/nf-core/raredisease) for full license text.


# if no header -> stop the pipeline
# if no assay -> stop the pipeline
#
# should  have valid sampleid
# should have valid fastq files and link
# if paired check for tumor and normal
'''
./check_samplesheet.py -c headerless.csv -o samplecheck.txt
'''

import csv
import sys
import argparse

CMDAssay = [
                "gmslymphomav3-0",
                "GMSMyeloidv1-0",
                "gmssolidtumorv3-0",
                "PARPinhibv1-0",
                "tumwgs-hema",
                "tumwgs-solid"
            ]

TUMOR_TYPES  = ('tumor', 'T')
NORMAL_TYPES = ('normal', 'N')

def process_linescsv_file(test):
    with open(test, mode='r') as file:
        csvFile = csv.DictReader(file)
        nrows = len(list(csvFile))
        return (nrows)

def check_group_pairing(rows):
    groups = {}
    for row in rows:
        groups.setdefault(row["group"], []).append(row)

    for group, group_rows in groups.items():
        types = [r["type"] for r in group_rows]
        n_tumor  = sum(1 for t in types if t in TUMOR_TYPES)
        n_normal = sum(1 for t in types if t in NORMAL_TYPES)

        if len(group_rows) == 1:
            if n_tumor != 1:
                print(f"group '{group}': a single-row group must be a tumor sample "
                      f"(type in {TUMOR_TYPES}), got type(s) {types}", file=sys.stderr)
                return False
        elif len(group_rows) == 2:
            if n_tumor != 1 or n_normal != 1:
                print(f"group '{group}': a two-row group must contain exactly one tumor "
                      f"and one normal sample, got type(s) {types}", file=sys.stderr)
                return False
        else:
            print(f"group '{group}': expected 1 (tumor-only) or 2 (tumor+normal) rows, "
                  f"got {len(group_rows)}", file=sys.stderr)
            return False

    return True

def process_csv_file(test,nrows):
    testAssay = []
    testId = []
    testType = []
    rows = []

    with open(test, mode='r') as file:
        csvFile = csv.DictReader(file)

        if nrows == 0:
            return None
        else:
            for row in csvFile:
                if len(row["assay"]) == 0 and row["assay"][0] not in CMDAssay:
                    print("invalid or missing 'assay' value", file=sys.stderr)
                    return None
                else:
                    testAssay.append(row["assay"])

                if len(row["id"]) == 0:
                    print("missing 'id' value", file=sys.stderr)
                    return None
                else:
                    testId.append(row["id"])

                if len(row["type"]) == 0:
                    print("missing 'type' value", file=sys.stderr)
                    return None
                else:
                    testType.append(row["type"])

                rows.append(row)

        if not check_group_pairing(rows):
            return None

        return nrows, testAssay, testId, testType

def writeFile (result,output):
    if result is not None:
        with open(output, mode= 'w') as outFile:
            lineCount, testAssay, testId, testType = result
            outFile.write(str(lineCount))
            outFile.write(str(testAssay))
            outFile.write(str(testId))
            outFile.write(str(testType))
            return outFile

def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--csv', dest = 'csv',  default = "test.csv", help = "Sample Sheet csv for the nextflow")
    parser.add_argument('-o', '--out', dest = 'output', default = "result", help = "Input csv structure and content signal")

    args = parser.parse_args()
    inputCsv = args.csv
    outputCsv = args.output
    count = process_linescsv_file(inputCsv)

    result = process_csv_file(inputCsv,count)

    writeFile(result, outputCsv)

if __name__ == "__main__":
    Main()
