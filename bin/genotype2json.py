#!/usr/bin/env python3

import json
import argparse


def parse_genotype_file(genotype_file_path):
    genotype_dict = {}
    with open(genotype_file_path, 'r') as genotype_file:
        for line in genotype_file:
            parts = line.split()
            if len(parts) == 2:
                genotype_dict[parts[0]] = parts[1]
            else:
                genotype_dict[parts[0]] = None  
    return genotype_dict


def main():
    parser = argparse.ArgumentParser(description="Convert genotype data file into a JSON blob.")
    parser.add_argument('genotype_file', type=str, help="Path to the genotype file")
    parser.add_argument('output_file', type=str, help="Path to save the output JSON file")
    
    args = parser.parse_args()

    genotype_data = parse_genotype_file(args.genotype_file)
    

    with open(args.output_file, 'w') as output_file:
        json.dump(genotype_data, output_file, indent=4)
    

if __name__ == '__main__':
    main()
