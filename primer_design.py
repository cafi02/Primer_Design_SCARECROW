"""
Python Script to Download Fasta Files of cDNA and Design Primers for the Coding Region
Author: Carlotta Fiedler
Last updated: 17.04.2023
"""
# ==================================================IMPORT PACKAGES=====================================================

import os
from Bio import Entrez
import urllib.error
import sys
import datetime


# =============================================DEFINE HELPFUL FUNCTIONS=================================================

def separate_codons(sequence):
    """This function takes a sequence and inserts a whitespace after every third base to display separated codons:"""
    codon_seq = ""
    for i in range(0, len(sequence), 3):
        codon_seq += sequence[i:i + 3] + " "
    return codon_seq.strip()

def calculate_annealing_temp(seq):
    """This function calculates the annealing temperature between two complementary nucleotide strands."""
    temp = 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))
    return temp

def input_default(prompt, default):
    """This function takes a prompt and a default answer in order to return a default value when the input was empty:"""
    # Prompt user for input
    user_input = input(str(prompt))
    # Check if user input is empty
    if not user_input:
        # Set default value
        user_input = str(default)
    return user_input

def gc_content(seq):
    """Calculates the percentage of "C" and "G" in a given DNA or RNA sequence."""
    if len(seq) != 0:
        seq = seq.upper()
        gc_count = seq.count('C') + seq.count('G')
        percentage = gc_count / len(seq) * 100
        return round(percentage,2)
    else:
        return None


def print_append(holo, s):
    print(s)
    holo += "\n" + s
    return holo

# ==========================================DOWNLOAD FASTA FILES FROM NCBI==============================================

# Enter email for NCBI access
email = input_default("Please enter an email address: ", "example@email.de")

# Access NCBI with the email
try:
    Entrez.email = email
except ValueError:
    print("Invalid email, please try again.")


# Define the search terms (NCBI-Ids)
gene_list = []
input_genes = input_default(
    "Please enter one or multiple (comma-seperated) NCBI Reference Sequence (ex. NM_115282.4, NM_103925.6): ",
    "NM_115282.4, NM_103925.6")
try:
    gene_list = input_genes.split(", ")
    # check if each item in gene_list is a valid gene ID
    for ncbi_id in gene_list:
        if not ncbi_id.startswith("NM_"):
            raise ValueError(f"Invalid NCBI ID: {ncbi_id}")
        break
except ValueError as e:
    sys.exit(str(e))


# Fetch fasta files with the given references from NCBI and save at the current directory
for ncbi_id in gene_list:

    filename = f"{ncbi_id}.fasta"
    if not os.path.isfile(filename):

        # Fetch fasta file from NCBI
        try:
            net_handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")

            # Create new file in the current directory and parse the content from the fetched fasta file
            f = open(filename, "w")
            f.write(net_handle.read())
            f.close()
            net_handle.close()
            print(f"Saved {filename} in the current directory.")

        except urllib.error.HTTPError:
            sys.exit(f"Invalid NCBI ID: {ncbi_id} not found.")



# ==================================================DESIGN PRIMERS======================================================

# Create an export string to save the generated primers in the end
export = "Primer Design\n" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M") +\
         "\n______________________________________________________________________" + "\nDesign-Criteria:\n"

# Define the parameters for the primer design
min_len = int(input_default("\nEnter min. primer length (ex. 20): ", 20))
max_len = int(input_default("Enter max. primer length (ex. 23): ", 23))
export += f"Range of primer length:\tmin. {min_len}, max. {max_len}\n"

min_temp = float(input_default("Enter min. annealing the temperature in °C (ex. 55.0): ", 55))
max_temp = float(input_default("Enter max. annealing the temperature in °C (ex. 62.0): ", 62))
export += f"Range of annealing temperatures:\tmin. {min_temp} °C, max. {max_temp} °C\n"

tm_diff = float(input_default(
    "Enter max. difference in annealing temperatures between forward and reverse primers in °C (ex. 4.0): ", 4))
export += f"Maximal difference in annealing temperatures:\t{tm_diff} °C\n" + \
          "______________________________________________________________________"

# Read content from fasta files, extract coding regions and design primers
for ncbi_id in gene_list:

    # Read content from file
    f = open(f"{ncbi_id}.fasta")
    content = f.read()
    f.close()

    # Save header and strip sequence of new-line characters
    gene_info = content[:content.find("\n")]
    sequence = content[content.find("\n"):].replace("\n", "")

    # Find the coding region
    begin, start, stop = 0, 0, len(sequence)
    coding_region = ""

    for i in range(len(sequence)-3):
        # Find start codon
        if sequence[i:i+3] == 'ATG':
            start = i
            # Find stop codon in frame to the previous start codon
            for j in range(start, len(sequence)-3):
                if (j-start)%3 == 0 and sequence[j:j+3] in ['TAG', 'TGA', 'TAA']:
                    stop = j
        # Save the sequence between the start and stop as coding reqion if it is longer than the previously found region
        if len(sequence[start:stop+3]) > len(coding_region):
            coding_region = sequence[start:stop+3]


    # while begin <= len(sequence)-3:
    #
    #     # Find the start codon
    #     start = sequence[begin:].find('ATG') + begin
    #
    #     # Find the stop codon
    #     stop = -1
    #     for codon in ['TAG', 'TGA', 'TAA']:
    #         if sequence[start:].find(codon) != -1:
    #             stop_codon_pos = sequence[start:].find(codon)
    #             if stop_codon_pos % 3 == 0:
    #                 stop = start + stop_codon_pos
    #
    #     if len(sequence[start:stop]) > len(coding_region):
    #         coding_region = sequence[start:stop]
    #         begin = start

    # Save different primers and their annealing temperatures in dictionaries
    f = {}
    r = {}
    for i in range(min_len, max_len):

        for pos in range(start-(i-3), start+1):
            forward = coding_region[pos:pos + i]  # 5' --> 3' direction
            f[forward] = calculate_annealing_temp(forward)

        for pos in range(stop-(i-3), stop+1):
            reverse = coding_region[pos:pos+i] # 3' --> 5' direction
            reverse = reverse[::-1].translate(str.maketrans("ATCG", "TAGC"))  # 5' --> 3' direction
            r[reverse] = calculate_annealing_temp(reverse)


    # Find primer pair with the minimal difference in annealing temperatures
    min_diff = None
    min_seq_f = None
    min_seq_r = None

    for seq_f, temp_f in f.items():
        for seq_r, temp_r in r.items():
            diff = abs(temp_f - temp_r)
            if min_diff is None or diff < min_diff:
                min_diff = diff
                min_seq_f = seq_f
                min_seq_r = seq_r

    export = print_append(export, f'\n{gene_info}')
    export = print_append(export, f"Coding region: 5'-{separate_codons(coding_region[start:stop+3])}-3'\n")

    # Check predefined criteria
    criteria_checked = True
    if min_temp >= f[min_seq_f] or f[min_seq_f] >= max_temp and min_temp >= r[min_seq_r] or r[min_seq_r] >= max_temp:
        export = print_append(export,
            f"The annealing temperatures of both primers for {gene_info[:gene_info.find(' ')]} are not within the predefined range.")
        criteria_checked = False
    elif min_temp >= f[min_seq_f] or f[min_seq_f] >= max_temp:
        export = print_append(export,
            f"The annealing temperatures of the FORWARD primer for {gene_info[:gene_info.find(' ')]} is not within the predefined range.")
        criteria_checked = False
    elif min_temp >= r[min_seq_r] or r[min_seq_r] >= max_temp:
        export = print_append(export,
            f"The annealing temperatures of the REVERSE primer for {gene_info[:gene_info.find(' ')]} is not within the predefined range.")
        criteria_checked = False

    if min_diff > tm_diff:
        export = print_append(export,
            f"The difference in annealing temperatures of the primers for {gene_info[:gene_info.find(' ')]} exceeds the predefined maximum.")
        criteria_checked = False

    if not criteria_checked:
        export = print_append(export, "This is the best match for the primer pair:")
    else:
        export = print_append(export, f"This primer pair for {gene_info[:gene_info.find(' ')]} checks all predefined criteria:")

    # Print best primer pair matches
    export = print_append(export, f"FORWARD primer: 5'-{min_seq_f}-3' \t({f[min_seq_f]:.2f}°C, GC-content: {gc_content(min_seq_f)} %)")
    export = print_append(export, f"REVERSE primer: 5'-{min_seq_r}-3' \t({r[min_seq_r]:.2f}°C, GC-content: {gc_content(min_seq_r)} %)")
    export = print_append(export, f"Difference in annealing Temperature: \t{min_diff}°C\n")

# ==================================================EXPORT PRIMERS======================================================

save_output = input("Do you want to save the output? (y/n): ")

if save_output.lower() == "y":
    filename = "primers_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M") + ".txt"
    f = open(filename, "w", encoding="utf-8")
    f.write(export)
    f.close()
    print(f"The output has been saved to {filename}.")
else:
    print("The output has not been saved.")
