# primer_tools
Supporting scripts for primer design and validation.

## Generic scripts for data processing
### extract_ncbiseqs.py
Sometimes, when you download genomes in bulk from NCBI, they get saved into a bulk folder called "ncbi dataset". Within the folder, there is a subdirectory called "data", which holds all the FASTQ files organized into subdirectories by accession.

This script combines the chromosomes of each accession into one FASTA file named after the organism, and moves them into the data folder. (the script also needs to be run while you're in the data folder)

### oligo_expander.py
Takes a set of sequences containing non ACTG bases and creates a new file holding all the expanded versions of the sequence.

Sample input
\>sample_seq
ACGCATWCA

Sample output
\>sample_seq\~1
ACGCATACA
\>sample_seq\~2
ACGCATTCA

USAGE: oligo_expander.py <input.fasta>

### subset_fasta.py
Select "n" number of sequences from a large FASTA file, either randomly (-r) or evenly distributed across the file

USAGE: subset_fasta.py -f <fasta file> [options]

## Scripts for validation
### get_locations.py
This script extracts the most-likely amplicon from a reference genome given a set of forward, probe (optional), and reverse primers. The best location is determined based on minimum hamming-distance (HD) using a linear-complexity algorithm. 

The primers are saved in <primerfile>, reference in <referencefile> both FASTA-formatted. Amplicons are attempted for each record in the referencefile and written to output file. The exact alignment is written to stdout.

For generating internal controls for a primer kit, the (-f) flag is used to indicate how many flanking bases to add onto the ends of the amplicon.

BASIC USAGE: get_locations.py -i <primerfile> -r <referencefile>

### dimer_tester.py
Accepts a FASTA file <primer file> (resembling the primers/controls within a reaction kit) and does a complete biophysical/biochemical evaluation. Generates a summary file with homodimer formation (dG) and Tm values, as well as a heterodimer formation matrix in csv format.

Tm estimations are calculated based on Allawi-1997 publication, and dG estimations are calculated based on SantaLucia-1996 publications, identifying constants for thermodynamic properties through kmer profiling. Values are corrected based on IDT's OligoAnalyzer tool, using the following concentrations: 
- [oligo]=0.25uM
- [Na+]=20mM
- [Mg2+]=4.5mM
- [dNTP]=0.75mM

USAGE: dimer_tester.py <primer file>

### process_blastresults.py
If we want to quickly observe how sensitive a set of primers (fwd/prb/rev) are against its target, we can perform a BLAST search online. Note that this does not consider the positioning of forward/reverse primers to each other. It is a qualitative test to see whether or not the primer sequences align with the potential reference targets.

Generic BLASTn params:
- algorithm: megablast
- max e-value: 100 (specific) or 1000 (generic)
- word size: 16
- disable "automatically adjust parameters"

Instructions:
Download the hit table (csv) and descriptions table (csv) after performing the search for any primer
Navigate the the folder where the tables are found, then run the script. Alternatively, you can specify the absolute locations of the tables using the -ht and -dt flags
This script looks for accessions and names that all primers hit.
