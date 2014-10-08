#!/usr/bin/env python

# Author: Leah Roberts
# Affiliation: Scott Beatson Lab Group University of Queensland St Lucia
# Date: April 2014

##################### Script to count nucleotide length of fasta file ######################

### Description:

# This script takes in a fasta file as raw input and prints the fasta sequence ID and nucleotide length.
# This was designed for IS sequences, but can be used for any fasta file. 

### Usage:

# Type this command in the same directory as the fasta files you want to count:
# $ python Count_nt.py

# This will bring up the prompt "input .fa file: " 
# Type the name of the fasta file and hit enter. 
# This should print the sequence ID and length.


import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob


# Input the full name of the fasta file manually:
in_file = raw_input("input .fa file: ")

# Can also change the script so that it takes in a number of fasta file from a directory called 'Sequences':
# in_file = glob.glob('Sequences/*.fa')

for record in SeqIO.parse(in_file, "fasta"):
	cur = str(record.seq)
        ID = record.id
        total = len(cur)
                
        # print ID
        print ID
        print total
       


