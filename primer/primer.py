
def extract_seq_from_fasta(filepath):
    """This function opens a Fasta file, deletes the first line (header) and strips the content of all new line
     characters to yield a genomic sequence"""

    # read content from fasta file
    f = open(filepath)
    content = f.read()
    f.close()

    # read genomic sequence from content
    seq = content[content.find("\n"):]
    seq.strip("\n")

    # return genomic sequence
    return seq

# UNIT-TESTING



def extract_coding_region(seq):
    """This function finds the start codon ATG and then the next stop codon TAG, TGA, or TAA in the correct reading
    frame and extracts the coding region inbetween the start and the stop codon."""

    # find the start codon and delete the 5' overhang
    start = seq.index("ATG")

    # find the stop codon and delete the 3' overhang
    i = start
    stop = -1
    while i < len(seq):
        if seq[i:i+3] in ["TAG", "TGA", "TAA"]:
            stop = i
            break
        i += 3

    # return the coding region
    return seq[start:(stop+3)]

# UNIT-TESTING
gene = extract_coding_region("AGTACATGAAAGGGTTTCCCGAACCGTAAGTTTAAGCT")
if gene == "ATGAAAGGGTTTCCCGAACCGTAA":
    print("Correct")
else:
    print(gene)


def create_complement(seq):
    """This function creates the complementary strand in 3'-5' direction of a nucleotide sequence given in 5'-3'.
    The function is not case sensitive and will return the complementary sequence in all uppercase characters."""

    seq = seq.upper() # set all nucleotides to uppercase

    complement = "" # create an empty complementary strand in

    for nt in seq:
        if nt == "A":
            complement += "T"
        elif nt == "T":
            complement += "A"
        elif nt == "C":
            complement += "G"
        elif nt == "G":
            complement += "C"
        else:
            print("Check your input sequence")


def calculate_annealing_temp(seq):
    """This function calculates the annealing temperature between two complementary nucleotide strands."""
    temp = 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))
    return temp

# ERROR-HANDLING


def design_primer(seq):

    # creation of initial primers and temperatures
    forw = seq[:20]  # forward primer creates copies of the 5’-3’ strand
    rev = seq[-20:]  # reverse primer makes copies of the complementary 3’-5’ strand
    f_t = calculate_annealing_temp(forw) # initial annealing temperature for the forward primer
    r_t = calculate_annealing_temp(rev) # annealing temperature for the reverse primer

    # primer design
    while len(forw) < 23 and len(rev) < 23:

        if f_t < 55.0 or f_t > 62.0:
            forw = seq[:len(forw)+1]
            f_t = calculate_annealing_temp(forw)
        elif r_t < 55.0 or r_t > 62.0:
            rev = seq[-(len(rev)+1):]
            r_t = calculate_annealing_temp(rev)

        if abs(f_t - r_t) < 4.0:
            print(forw, f_t)
            print(rev, r_t)
        else:
            print(forw, f_t)
            print(rev, r_t)



# UNIT-TESTING



