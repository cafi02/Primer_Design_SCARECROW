import os
from Bio import Entrez

# Enter email for NCBI access
# email = input("Please enter your email: ")
# Entrez.email = email
Entrez.email = "c.fiedler02@stud.uni-goettingen.de"

# Define the search term
# ncbi_id = input("Please enter a NCBI Reference Sequence (ex. NM_115282.4): ")
ncbi_id = "NM_115282.4"

filename = f"{ncbi_id}.fasta"
if not os.path.isfile(filename):

    # Fetch fasta file from NCBI
    net_handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")

    # Create new file in the current directory and parse the contend from the fetched fasta file
    f = open(filename, "w")
    f.write(net_handle.read())
    f.close()
    net_handle.close()
    print(f"Saved {filename} in the current directory.")




