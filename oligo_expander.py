#!/home/sbsuser/miniconda3/bin/python3

import sys
import itertools

from pathlib import Path
from Bio import SeqIO

str_program = Path(__file__).name

print("This is honzo's {0} for expanding IUPAC nucleotide sequences.\n".format(str_program))
help_menu = "USAGE: {0} <fasta file>".format(str_program)
if ("-h" or "--help") in sys.argv:
    print(help_menu); quit()

def get_combinations(str_fullseq):
    dic_ambs = {"M":["A","C"], "R":["A","G"], "W":["A","T"], "S":["C","G"],
        "Y":["C","T"], "K":["G","T"], "V":["A","C","G"], "H":["A","C","T"],
        "D":["A","G","T"], "B":["C","G","T"], "N":["A","C","G","T"]}

    # we create a dictionary to figure out where all these ambiguous nucleotides are
    dic_amblocs = {}
    for i, char in enumerate(str_fullseq):
        if char in dic_ambs.keys():
            dic_amblocs[i] = char

    j = 0
    lst_segments = []
    #print("fullseq", str_fullseq)
    for i in dic_amblocs.keys():
        str_subseq = str_fullseq[j:i+1]
        #print(i, str_subseq)
        lst_subseqs = [str_subseq[:-1]+nuc for nuc in dic_ambs[str_subseq[-1]]]
        lst_segments.append(lst_subseqs)
        j = i+1

    str_endseq = str_fullseq[j:]

    lst_seqsegments = [seq for seq in itertools.product(*lst_segments)]
    lst_seqs = ["".join(seqsegment)+str_endseq for seqsegment in lst_seqsegments]
    return lst_seqs

pth_input = Path(sys.argv[-1])
pth_output = pth_input.stem + "_expanded.fasta"
with open(pth_input) as handle:
    lst_records = list(SeqIO.parse(handle, "fasta"))

with open(pth_output, "wt") as wf_out:
    for record in lst_records:
        print("Record [{0}]:\n{1}".format(record.id, record.seq))
        for i, str_expandedseq in enumerate(get_combinations(str(record.seq))):
            wf_out.write(record.id+"~{0}\n".format(i+1))
            wf_out.write(str_expandedseq+"\n")
        print(i, "sequences made.\n")

print("Program complete.")