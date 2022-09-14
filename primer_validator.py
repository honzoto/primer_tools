#!/home/sbsuser/miniconda3/bin/python

import sys, os
from pathlib import Path
from Bio import SeqIO

version = "1.1"
program = Path(__file__).name
print(f"Starting Honzo's {program} v{version} for validating Norgen primers")

# Variable declarations
tf_partialhd = False
threshold = 'max' # min, max, or a number between 0 and 100

global dic_ambs, dic_ambs_rev
dic_ambs = {"A": ["A"], "C":["C"], "T":["T"], "G":["G"], "-":["-"],
    "M":["A","C"], "R":["A","G"], "W":["A","T"], "S":["C","G"],
    "Y":["C","T"], "K":["G","T"], "V":["A","C","G"], "H":["A","C","T"],
    "D":["A","G","T"], "B":["C","G","T"], "N":["A","C","G","T"]}
dic_ambs_rev = {"A":"A", "C":"C", "T":"T", "G":"G",
    "AC":"M", "AG":"R", "AT":"W", "CG":"S", "CT":"Y", "GT":"K",
    "ACG":"V", "ACT":"H", "AGT":"D", "CGT":"B", "ACGT":"N"}


help_menu = """
----------------------------------[ HELP MENU ]---------------------------------

    USAGE: 
    python {n} -i <filename> -r <reference> [options]

    ARGUMENTS:
    -i/--input      : <path> input primer file (.fasta)
    -r/--reference  : <path> input reference file (.fasta)
    -t/--threshold  : [default={t}] if 'max', then writing to file as pct,
                        if number (e.g. 80), output will write <80% instead
                        of the actual number

    OPTIONS:
    -p/--partial_hd : [default={p}] allow partial hamming distance additions
                        with ambiguous primers otherwise recorded as mismatch

--------------------------------------------------------------------------------
""".format(p=tf_partialhd, t=threshold, n=program)

for i, arg in enumerate(sys.argv):
    if arg == ("-i" or "--input"):
        pth_primers = Path(sys.argv[i+1])
    elif arg == ("-r" or "--reference"):
        pth_reference = Path(sys.argv[i+1])
    elif arg == ("-t" or "--threshold"):
        threshold = sys.argv[i+1]
    elif arg == ("-h" or "--help"):
        print(help_menu); quit()
    elif arg[0] == "-":
        print("[Warning] Unrecognized argument: {0}".format(arg))

# ============================================[ FUNCTIONS ]===========================================

class Query:
    def __init__(self, id, seq):
        self.id = id
        self.seq = str(seq).upper()
        self.revc = revc(self.seq)

    def input_results(self, index, hd, source, alignment, strand):
        self.index = index
        self.hd = hd
        self.source = source
        self.alignment = alignment
        self.strand = strand

def revc(seq):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
        'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K',
        'B':'V', 'V':'B', 'D':'H', 'H':'D', 'N':'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def ham_dist(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Undefined")
    hd = sum(chain1 != chain2 for chain1, chain2 in zip(seq1, seq2))
    return hd

def ham_dist_ambs(seq1, seq2):
    # arguments entered as (search, source)
    if len(seq1) != len(seq2):
        raise ValueError("Undefined")

    hd = 0
    for nuc1, nuc2 in zip(seq1, seq2):
        com1 = set(dic_ambs[nuc1])
        com2 = set(dic_ambs[nuc2])
        hd += len(com1^com2) / len(com1|com2)

def search_min_dist(reference, query, amb_hd=False):
    # accepts query object and reference, adds attributes to query
    reference = str(reference).upper()
    dists = []
    
    index = 0
    min_dist = 999
    min_substring = ""

    for search in [query.seq, query.revc]:
        l = len(search)
        for i in range(len(reference)-l+1):
            if amb_hd:
                d = ham_dist_ambs(search, reference[i:i+l])
            else:
                d = ham_dist(search, reference[i:i+l])

            if d < min_dist: 
                min_dist = d
                index = i
                min_substring = reference[i:i+l]

        dists.append(min_dist)

    if dists[0] <= dists[1]:
        aln = ''.join([' ' if n1 == n2 else '*' for n1, n2 in zip(query.seq, min_substring)])
        query.input_results(index=index, hd=min_dist, source=min_substring, alignment=aln, strand="sense")
    else:
        aln = ''.join([' ' if n1 == n2 else '*' for n1, n2 in zip(query.revc, min_substring)])
        query.input_results(index=index, hd=min_dist, source=min_substring, alignment=aln, strand="antisense")
    
    return query

# ============================================[ WORKFLOW ]============================================

print("Parameters:", sys.argv)

with open(pth_primers) as handle:
    print("Reading primer input: {0}".format(pth_primers))
    primer_fastas = list(SeqIO.parse(handle, "fasta"))

with open(pth_reference) as handle:
    print("Reading reference input: {0}".format(pth_reference))
    reference_fastas = list(SeqIO.parse(handle, "fasta"))

pth_output = pth_primers.parent / ("{i}-{r}_summary.csv".format(i=pth_primers.stem, r=pth_reference.stem))

with open(pth_output, "wt") as writer_out:
    #writer_out.write("reference,primer,{0}_homology\n".format(threshold))

    for primer in primer_fastas:
        print("\nProcessing primer: {0}".format(primer.id))
        query = Query(id=primer.id, seq=primer.seq)
        highest_homology = 0.00

        for reference in reference_fastas:
            search_min_dist(reference, query, amb_hd=tf_partialhd)
            current_homology = 1 - (query.hd / len(query.seq))
            if current_homology > highest_homology:
                highest_homology = current_homology
                record_reference = reference.description.split(",")[0]

        if threshold == 'max':
            writer_out.write("{0},{1},{2}%\n".format(record_reference, primer.id, round(highest_homology*100)))
        else: # threshold is an integer
            if highest_homology*100 < int(threshold):
                writer_out.write("{0},{1},<{2}%\n".format(record_reference, primer.id, threshold))
            else:
                writer_out.write("{0},{1},{2}%\n".format(record_reference, primer.id, round(highest_homology*100)))


print("Program complete.")