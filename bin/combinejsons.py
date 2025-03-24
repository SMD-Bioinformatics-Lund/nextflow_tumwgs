#!/usr/bin/env python
import json
import argparse

# Function to load JSON data from a file
def load_json_file(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Function to save JSON data to a file
def save_json_file(data, file_path):
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)

# Main function to combine the Info JSON and the genotype JSON
def combine_json_files(info_json_file, genotype_json_file, partner_run_json_file, output_file):
    # Load the Info and genotype JSON files
    info_data = load_json_file(info_json_file)
    genotype_data = load_json_file(genotype_json_file)
    
    # Combine the data by appending the genotype JSON
    combined_data = {
        **info_data,
        "id_snp_genotypes": genotype_data
    }

    # If partner_run_json_file is provided, include it in the combined data
    if partner_run_json_file:
        partner_run = load_json_file(partner_run_json_file)
        combined_data["partner_info"] = partner_run

    # Save the combined data to the output file
    save_json_file(combined_data, output_file)
    print(f"Combined JSON saved to {output_file}")

# Argument parser setup
def parse_args():
    parser = argparse.ArgumentParser(description='Append genotype JSON into Info JSON and save the result.')
    
    # Add arguments for the input JSON files and the output file
    parser.add_argument('info_json_file', help='Path to the Info JSON file')
    parser.add_argument('genotype_json_file', help='Path to the genotype JSON file')
    parser.add_argument('--partner_run_json_file', help='Path to the partner sample run JSON file (optional)', default=None)
    parser.add_argument('output_file', help='Path to the output combined JSON file')

    return parser.parse_args()

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()

    # Call the function to combine the JSON files
    combine_json_files(args.info_json_file, args.genotype_json_file, args.partner_run_json_file, args.output_file)
