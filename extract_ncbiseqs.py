#!/home/sbsuser/miniconda3/bin/python

import os, sys, glob, shutil

dir_data = os.getcwd()
if len(sys.argv) > 1:
    dir_data = sys.argv[-1]

print("Running honzo's sequence-extractor at: {0}".format(dir_data))
for dir in [d for d in os.listdir(".") if os.path.isdir(d)]:
    os.chdir(dir)
    print("Working in directory: {0}".format(dir))

    # get the name of the organism (to later use as the filename)
    lst_fastas = glob.glob("*.fna")
    with open(lst_fastas[0]) as rf_fasta:
        str_defline = rf_fasta.readline()
        str_filename = "_".join(str_defline.split()[0:3]).replace(">","")+".fasta"
    os.system("cat *.fna > {0}".format(str_filename))
    shutil.move(str_filename, "../"+str_filename)
    os.chdir("..")


print("Program complete.")