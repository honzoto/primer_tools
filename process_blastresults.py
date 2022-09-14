#!/home/sbsuser/miniconda3/bin/python


import sys, os, glob
from pathlib import Path
import pandas as pd

prog_name = Path(__file__).name
version = "1.1"

# if these are left blank, program will auto-detect hits/description files based on the current directory
pth_ht = ""
pth_dt = ""
pth_output = "hits_sn.txt"

help_menu = """
----------------------------------[ HELP MENU ]---------------------------------

    USAGE: 
    python {n} -ht <filename> -r <reference> [options]

    ARGUMENTS:
    -ht/--hits_table    : <path> input hits table (from BLAST search)
    -dt/--desc_table    : <path> input descriptions table (from BLAST search)
    -h/--help           : shows this menu.

--------------------------------------------------------------------------------
"""

print("Starting Honzo's {0} v{1} for processing BLAST outputs".format(prog_name, version))
print("Working in directory: {0}".format(os.getcwd()))
if len(sys.argv) < 2: 
    print(help_menu); quit()

for i, arg in enumerate(sys.argv):
    if arg == ("-h" or "--help"):
        print(help_menu); quit()
    elif arg == ("-ht" or "--hits_table"):
        pth_ht = Path(sys.argv[i+1])
    elif arg == ("-dt" or "--desc_table"):
        pth_dt = Path(sys.argv[i+1])
    elif arg[0] == "-":
        print("[Warning] unrecognized argument:", arg)

if "." in sys.argv:
    pth_ht = glob.glob("*HitTable.csv")[0]
    pth_dt = glob.glob("*Descriptions.csv")[0]

# Preliminary Processing
print("Hit table    : {0}".format(pth_ht))
# https://www.biostars.org/p/275241/
# changed 'subject ids' to Accession
str_cnames = "query id, Accession, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score"
lst_cnames = str_cnames.split(", ")
df_hits = pd.read_csv(pth_ht, header=None, names=lst_cnames)

print("Descriptions : {0}".format(pth_dt))
df_desc = pd.read_csv(pth_dt).dropna()
df_desc.columns = [v.strip() for v in df_desc.columns]

lst_pnames = sorted(list(set(df_hits["query id"])))
print("\nPrimers: {0}".format(lst_pnames))

# Starting workflow
# get accessions that are hit by forward, probe, reverse primers
dic_accessions = {}
for pname in lst_pnames:
    accessions = set(df_hits[df_hits["query id"] == pname]["Accession"])
    print("---> {0} unique hits found for primer: {1}".format(len(accessions), pname))
    dic_accessions[pname] = accessions

print("\nGetting accessions from all queries")
set_common = set.intersection(*dic_accessions.values())
print("---> {0} hits are common between all {1} primers".format(len(set_common), len(lst_pnames)))
df_commonhits = df_hits[df_hits["Accession"].isin(set_common)]
if len(set_common) > 5:
    print("First 5 hits: {0}, {1}, {2}, {3}, {4}".format(*list(set_common)[:5]))
else:
    print("Hits: {0}".format(list(set_common)))

df_combined = df_desc.merge(df_commonhits, on="Accession")
df_combined.to_csv("hits_combined.csv", index=False)
print("\nCombined dataframe:")
print(df_combined)

# Getting unique strains
print("\nGathering scientific names of hits")
set_unique = set(df_combined["Scientific Name"])

print("---> {0} unique names found".format(len(set_unique)))
print("Writing to file: {0}".format(pth_output))
with open(pth_output, "wt") as wf_out:
    for sn in set_unique:
        wf_out.write(sn+"\n")

print("\nProgram complete.")