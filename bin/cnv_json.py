import json
import argparse

def parse_tabular_to_json(lines):
    """Parse tabular data from a string and return a JSON-formatted string."""

    results = {}
    
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 8:
            continue
        try:
            chr_num, start, end, ratio = parts[0], int(parts[1]), int(parts[2]), float(parts[4])
            genes_info_list = parts[6].split(',') if parts[6] != '-' else []
            genes_annotation_list = parts[7].split(',') if parts[7] != '-' else []
            # Constructing the genes list
            genes = []
            for gene_info in genes_info_list:
                genes.append({"gene": gene_info})

            # Constructing the annotation list
            for gene_annotation in genes_annotation_list:
                annotation_parts = gene_annotation.split(':')
                genes.append({
                    "gene": annotation_parts[0],
                    "class": annotation_parts[1],
                    "cnv_type": annotation_parts[2]
                })

            results[f"{chr_num}:{start}-{end}"] = {
                "callers": "gatk",
                "ratio": ratio,
                "size": abs(end - start),
                "chr": chr_num,
                "start": start,
                "end": end,
                "genes": genes
            }
        
        except (ValueError, IndexError):
            print(f"Skipping malformed line: {line}")
    
    return json.dumps(results, indent=4)


def read_file(file_path):
    try:
        with open(file_path, "r") as file:
            return file.readlines()
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return []

def write_json(output_path, json_data):
    try:
        with open(output_path, "w") as file:
            file.write(json_data)
        print(f"JSON file successfully created: {output_path}")
    except IOError:
        print(f"Error writing to {output_path}.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", type=str, help="Set this flag to True")
    parser.add_argument("--json", type=str, help="Set this flag to True")
    args = parser.parse_args()

    input_file: str = args.bed
    output_file: str = args.json

    lines = read_file(input_file)

    if lines:
        json_output = parse_tabular_to_json(lines)
        write_json(output_file, json_output)

if __name__ == "__main__":
    main()


    # File paths
    #input_file = "/data/bnf/dev/saile/other/prj/results/tumwgs_validation/5218-21-validation.cnv.annotated.pane.bed"
    #output_file = "/data/bnf/dev/saile/other/prj/results/tumwgs_validation/5218-21-validation.cnv.annotated.pane.json"

# Read from a file and process the data







