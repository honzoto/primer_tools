#!/home/sbsuser/miniconda3/bin/python

import sys, os
from pathlib import Path

from Bio import SeqIO

str_version = "4.3"
str_program = Path(__file__).name

"""
< PURPOSE >
The objective of this program is to determine where primers fit best in a genome sequence.

< VERSION UPDATES >
v1  - 2022-05-09 created program after suggestion from alex to perform primer alignments
    - if ID ends with "_p" or _r", check if reverse complement gets you lower mismatch
    - we are now checking both sense and antisense notations and keeping the one that matches better
v2  - created Query object to store information for alignments instead of list
    - we are no longer using a log file, instead writing straight to an output file
    - added support for primers with ambiguous bases
    - fixed bug where primers have to go in order of forward, probe, reverse for the frames to display correctly
    - removed file logging. program will now instead print to termmhinal
    - added option (-r) to isolate a subset of the entire gene by coordinates
    - added query name in addition to reference name in amplicon defline
v3  - support for ambiguous bases in reference, where partial scores are given to partial matches
    - option to not write to file where the primers are binding to with ()[] characters
    - prints alignment to file, with additional line to see where mismatches are (similar to clustalw webtool)
    - [bug-fix] was printing the forward strand of the reverse primer to the log file
    - update to support for ambiguous bases, reverse-complements are now 
    
v4  - attempts to get top hits, assemble them to generate most-likely amplicon
    - option to show description in printout instead of writing just the accession ID of reference
    - print primer alignments to terminal even when the amplicon size is too large

    - added +1 to primer naming scheme
"""

print("Starting honzo's {0} v{1} for determining primer positions.".format(str_program, str_version))

# ========================================[ GETTING ARGUMENTS ]=======================================

tf_modified = False
tf_partialhd = False
tf_showbinding = False
tf_writedesc = False
int_maxampsize = 1000
int_flank = 0
int_maxampnum = 5
region = "all"

help_menu = """
----------------------------------[ HELP MENU ]---------------------------------

    USAGE: 
    python get_locations.py -f <filename> -r <reference> [options]

    ARGUMENTS:
    -i/--input      : <path> input primer file (.fasta)
    -r/--reference  : <path> input reference file (.fasta)

    OPTIONS:
    -f/--flank      : [default={f}] include bases on the ends of the amplicon
                        in the output file
    -p/--partial_hd : [default={p}] allow partial hamming distance additions
                        with ambiguous primers otherwise recorded as mismatch
    -e/--region     : [default={e}] choose coordinates (0-based) to perform
                        alignment search (e.g. "-r 15000,15800")
    -a/--ampsize    : [default={a}] largest allowed amplicon size to write to 
                        file if the best amplicon position is larger than this
                        number, it will show as (no amplicon found)
    -b/--show_bind  : [default={b}] show where primer binds to the amplicon,
                        '(amp)' for sense, '[amp]' for antisense strand
    -d/--desc       : [default={d}] write description in addition to record ID

--------------------------------------------------------------------------------
""".format(a=int_maxampsize, m=tf_modified, f=int_flank, e=region, p=tf_partialhd, \
    b=tf_showbinding, d=tf_writedesc)

# User arguments
if len(sys.argv) > 1:
    print("[gl] Retrieving user arguments...")
    args = sys.argv

    for i in range(len(args)):
        if args[i] == ("-i" or "--input"):
            pth_primers = Path(args[i+1])
        elif args[i] == ("-r" or "--reference"):
            pth_reference = args[i+1]
        elif args[i] == ("-p" or "--partial_hd"):
            tf_partialhd = True
        elif args[i] == ("-e" or "--region"):
            region = [int(e) for e in sys.argv[i+1].split(",")]
        elif args[i] == ("-a" or "--ampsize"):
            int_maxampsize = int(args[i+1])
        elif args[i] == ("-f" or "--flank"):
            int_flank = int(args[i+1])
        elif args[i] == ("-b" or "--show_bind"):
            tf_showbinding = True
        elif args[i] == ("-d" or "--desc"):
            tf_writedesc = True

        elif args[i] == ("-h" or "--help"):
            print(help_menu); quit()
        elif args[i][0] == "-":
            print("[WARNING] unrecognized flag: {0}".format(args[i]))

else:
    print("[ERROR] No arguments were entered.")
    print(help_menu)

# Setting variables for absolute paths
dir_project = pth_primers.parent
pth_primers = dir_project / pth_primers
pth_reference = dir_project / pth_reference

# Global variable delcarations
global dic_ambs, dic_ambs_rev
dic_ambs = {"A": ["A"], "C":["C"], "T":["T"], "G":["G"], "-":["-"],
    "M":["A","C"], "R":["A","G"], "W":["A","T"], "S":["C","G"],
    "Y":["C","T"], "K":["G","T"], "V":["A","C","G"], "H":["A","C","T"],
    "D":["A","G","T"], "B":["C","G","T"], "N":["A","C","G","T"]}
dic_ambs_rev = {"A":"A", "C":"C", "T":"T", "G":"G",
    "AC":"M", "AG":"R", "AT":"W", "CG":"S", "CT":"Y", "GT":"K",
    "ACG":"V", "ACT":"H", "AGT":"D", "CGT":"B", "ACGT":"N"}

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

    return round(hd, 2)
    

# ============================================[ WORKFLOW ]============================================

print("Your parameters: {0}".format(" ".join(args)))

os.chdir(dir_project)
print("\nWorking directory: {0}".format(dir_project))
if tf_showbinding:
    pth_out = dir_project / (pth_primers.stem + "_amps.txt")
else:
    pth_out = dir_project / (pth_primers.stem + "_amps.fasta")
# reset the output file
open(pth_out, "wt").close()

# Header
print("Maximum amplicon size: {0}".format(int_maxampsize))
if not region == "all":
    print("Genome region: start={0}, end={1}".format(region[0], region[1]))

with open(pth_reference) as handle:
    lst_references = list(SeqIO.parse(handle, "fasta"))

with open(pth_primers) as handle:
    # currently, each primer file should only contain one set of primers
    lst_primers = list(SeqIO.parse(handle, "fasta"))

# --------------------------------------[ Generating alignments ]-------------------------------------

for r, rec_reference in enumerate(lst_references):
    if tf_writedesc: str_defline = rec_reference.description
    else: str_defline = rec_reference.id
    print("[gl] Reading reference-{0}: {1}".format(r+1, str_defline))
    lst_queries = []

    if region != "all":
        rec_reference.seq = rec_reference.seq[region[0]:region[1]]
        print("hz175", len(rec_reference.seq))

    for rec_primer in lst_primers:
        # since primers are all written in 5>3 notation, we don't know if it is sense or antisense yet
        query = Query(rec_primer.id, rec_primer.seq)
        search_min_dist(rec_reference.seq, query=query, amb_hd=tf_partialhd)
        lst_queries.append(query)

    int_ampstart = min([q.index for q in lst_queries]) - int_flank
    int_ampend = max([q.index + len(q.seq) for q in lst_queries]) + int_flank

    # to insert '[]' characters where queries are hitting, we must create a list cos strings are immutable
    lst_amp = list(rec_reference.seq)[int_ampstart:int_ampend]
    int_ampsize = len(lst_amp)

    str_entry = "Amplicon region: {0}-{1}, size={2}\n".format(int_ampstart, int_ampend, int_ampsize)

    bracket_positions = {}
    corrected_positions = {}

    for i, query in enumerate(lst_queries):
        start = query.index - int_ampstart
        end = start + len(query.seq) #+ 1

        # see if we can condense this chunk later
        if tf_showbinding:
            if query.strand == "sense":
                if start in bracket_positions: 
                    bracket_positions[start] += "("
                else: 
                    bracket_positions[start] = "("
                if end in bracket_positions: 
                    bracket_positions[end] = ")"+bracket_positions[end]
                else: 
                    bracket_positions[end] = ")"

            elif query.strand == "antisense":
                if start in bracket_positions: 
                    bracket_positions[start] += "["
                else: 
                    bracket_positions[start] = "["
                if end in bracket_positions: 
                    bracket_positions[end] = "]"+bracket_positions[end]
                else: 
                    bracket_positions[end] = "]"

        str_entry += "Primer-{a}: {id}\nStrand: {sd} / [{s}-{e}] / HD={hd}\n" \
            .format(a=i+1, id=query.id, sd=query.strand, hd=query.hd, s=query.index, e=query.index+len(query.seq))
        if query.strand == "sense":
            str_entry += "query     : {0}\n".format(query.seq)
        elif query.strand == "antisense":
            str_entry += "query     : {0}\n".format(query.revc)
        str_entry += "reference : {0}\n".format(query.source)
        str_entry += "alignment : {0}\n\n".format(query.alignment)

    if tf_showbinding:
        for i, pos in enumerate(sorted(bracket_positions.keys())):
            corrected_positions[pos+i] = bracket_positions[pos]
            
        for pos in sorted(corrected_positions.keys()):
            lst_amp.insert(pos, corrected_positions[pos])

    str_queries = "|".join([p.id for p in lst_primers])
    str_amp = "".join(lst_amp)

    # write sequence to fasta file
    with open(pth_out, "at") as wf_out:
        if int_ampsize > int_maxampsize:
            print("[{1}] Amplicon size too large: {0}".format(int_ampsize, rec_reference.id))
            wf_out.write(">{q}:{r}\n".format(q=pth_primers.stem, r=str_defline+"_amp"))
            wf_out.write("(best amplicon not within specifications)\n")
        else:
            print("[{1} Amplicon size: {0}".format(int_ampsize, rec_reference.id))
            wf_out.write(">{q}:{r}\n".format(q=pth_primers.stem, r=str_defline+"_amp"))
            wf_out.write(str_amp+"\n")
        print(str_entry)


print("Program complete.\n")
