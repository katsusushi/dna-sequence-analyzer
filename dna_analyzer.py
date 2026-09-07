"DNA Sequence Analyzer Module"

import os 
from Bio import Entrez
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 


                    #CODON TABLE
codon_table = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

valid_amino_acids = "ARNDCEQGHILKMFPSTWYV"

Entrez.email = "youremail@gmail.com"

def clean_sequence(seq):
    seq = seq.strip()
    seq = seq.upper()
    seq = seq.replace(" ", "")
    seq = seq.replace("\n", "")
    return seq

def is_valid_dna(seq):
    valid_bases = "ATGC"
    for base in seq:
        if base not in valid_bases:
            return False 
    return True

def count_nucleotides(seq):
    counts = {}
    counts["A"] = seq.count("A")
    counts["T"] = seq.count("T")
    counts["G"] = seq.count("G")
    counts["C"] = seq.count("C")
    return counts

def gc_content(seq):
    g_count = seq.count("G")
    c_count = seq.count("C")
    total = len(seq)
    gc = (g_count + c_count) / total * 100
    return gc

def reverse_complement(seq):
    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    complement = ""
    for base in seq: 
        complement = complement + complement_map[base]
    reversed_complement = complement[::-1]
    return reversed_complement

def transcribe(seq):
    seq = seq.replace("T", "U")
    return seq

def translate(seq, frame=0):
    protein = ""
    for i in range(frame, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            break
        amino_acid = codon_table[codon]
        if amino_acid == "*":
            break
        else:
            protein = protein + amino_acid 
    return protein

def read_fasta(filename):
    file = open(filename, "r")
    header = ""
    sequence = ""

    for line in file:
        line = line.strip()
        if line.startswith(">"):
            header = line
        else: 
            sequence = sequence + line 

    file.close()
    return header, sequence

def fetch_from_ncbi(accession_number):
    handle = Entrez.efetch(db = "nucleotide", id = accession_number, rettype = "fasta", retmode = "text")
    record = handle.read()
    handle.close()
    return record 

def plot_nucleotide_counts(counts):
    plt.bar(counts.keys(), counts.values())
    plt.xlabel("Nucleotide")
    plt.ylabel("Count")
    plt.title("Nucleotide Count: ")
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    plt.savefig(os.path.join(BASE_DIR, "static", "chart.png"))
    plt.close()

def plot_protein_count(counts):
    plt.bar(counts.keys(), counts.values())
    plt.xlabel("Amino Acid")
    plt.ylabel("Count")
    plt.title("Amino Acid Count: ")
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    plt.savefig(os.path.join(BASE_DIR, "static", "chart.png"))
    plt.close()

def six_frame_translation(seq):
    results = {}

    #forward the strands, 3 frames 
    rna = transcribe(seq)
    results["Frame +1"] = translate(rna, 0)
    results["Frame +2"] = translate(rna, 1)
    results["Frame +3"] = translate(rna, 2)

    #reverse complement, 3 frames
    rev_comp = reverse_complement(seq)
    rev_comp_rna = transcribe(rev_comp)
    results["Frame -1"] = translate(rev_comp_rna, 0)
    results["Frame -2"] = translate(rev_comp_rna, 1)
    results["Frame -3"] = translate(rev_comp_rna, 2)

    return results

def is_valid_protein(seq):
    valid_amino_acids = "ARNDCEQGHILKMFPSTWYV"
    for amino_acid in seq:
        if amino_acid not in valid_amino_acids:
            return False
    return True
        
def count_amino_acid(seq):
    counts = {}
    for amino_acid in seq:
        if amino_acid in counts:
            counts[amino_acid] = counts[amino_acid] + 1
        else: 
            counts[amino_acid] = 1
    return counts 

def analyze(seq):
    cleaned = clean_sequence(seq)

    if not is_valid_dna(cleaned):
        report = "ERROR: Sequence contains invalid characters (only A, T, G, C allowed)" 
        return report

    counts = count_nucleotides(cleaned)
    gc = gc_content(cleaned)
    rev_comp = reverse_complement(cleaned)
    rna = transcribe(cleaned)

    report = "\n" + "• DNA SEQUENCE ANALYSIS REPORT •\n"
    report += "\n" + "Sequence: " + cleaned + "\n"
    report += "Length: " + str(len(cleaned)) + " bp" + "\n"
    report += "Nucleotide Counts: " + str(counts) + "\n"
    report += "GC Content: {:.2f}%".format(gc) + "\n"
    report += "Reverse Complement: " + rev_comp + "\n"
    report += "RNA(transcribed): " + rna + "\n"
    report += "Protein(translated): " + translate(rna) + "\n"
    return report 
 
def analyze_protein(seq):
    cleaned = clean_sequence(seq)

    if not is_valid_protein(cleaned):
        report = "ERROR: Sequence contains invalid charecters (only the 20 valid letters allowed!)" 
        return report

    counts = count_amino_acid(cleaned)
    protein_length = len(cleaned)

    dna_like = set(cleaned.upper()) <= {"A", "T", "G", "C"}

    report = "" 
    if dna_like:
        report += "--- NOTE --- \n This sequence only contains A ,T, G ,C - it might actually be a DNA sequence.\n Double check you selected the right type.\n\n"

    report += "• PROTEIN SEQUENCE ANALYSIS REPORT •" + "\n"
    report += "\n" + "Protein Sequence: " + cleaned + "\n"
    report += "Length: " + str(protein_length) + "\n"
    report += "Amino Acid Count: " + str(counts) + "\n"
    return report 
    



# =========================================
#           TEST / DEMO SECTION
# =========================================

if __name__ == "__main__":

    # --- Fetch a real gene from NCBI and run the full pipeline ---
    result = fetch_from_ncbi("NM_000207")

    output_file = open("insulin.fasta", "w")
    output_file.write(result)
    output_file.close()

    header, sequence = read_fasta("insulin.fasta")
    analyze(sequence)



    # --- Multiple reading frame translation demo ---
    print(translate("AUGGCCUAG", 0))
    print(translate("AUGGCCUAG", 1))
    print(translate("AUGGCCUAG", 2))
    
    # --- Older stage tests, kept for reference ---
    #   analyze("GGCCTAAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC")
    #   header, sequence = read_fasta("test_sequence.fasta")
    #   print(header)
    # analyze(sequence)

    # --- 6-Frame Translation ---
    frames = six_frame_translation("ATGGCCTAG")
    print(frames)

    # --- Protein/Amino Acid Sequence Input --- 
    print(is_valid_protein("MAIVMGR"))
    print(is_valid_protein("MXYZ"))

    print(count_amino_acid("MAIVMGR"))
    analyze_protein("MAIVMGR")

    # --- Visualization --- 
    counts = count_nucleotides("ATGGCCATTGTAATGGGCCGCTGA")
    plot_nucleotide_counts(counts)

    # --- Flask test --- 
    print(analyze("ATGGCCATTGTAATGGGCCGCTGA"))
    print(analyze_protein("MAIVMGR"))