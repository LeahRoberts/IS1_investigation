#!/usr/bin/env python

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob

# Script to count the number of "N"s in consensus mapping data generated for SeqFindR. 
# Total number of N's/total base pairs will give an indication of whether or not an IS is actually in that position.
# i.e. a large percentage of N's means that not much of that sequence aligned to the consensus, meaning that the IS is probably not there.

# Input the full name of the fasta file:
# If you want to try it on just one file (or try to make it loop):

#in_file = open("consensus_mapping_fa_files/B36EC_consensus.fa", "rU")

# Otherwise just input the full name of the fasta file:

#for in_file in glob.glob('IS_100_flank/*.fa'):
in_file = raw_input("input fa. file: ")
for record in SeqIO.parse(in_file, "fasta"):
	cur = str(record.seq)
        ID = record.id
        total = len(cur)
                
        # print ID
        print ID
        print total
       # print subset "/" total

# Prints a percentage of N in the total number of base pairs.
# Couldn't figure out how to get it to loop through the files, but this way turned out better for copying the output into excel.
# If you used the "in_file" command above, you need to uncomment the in_file.close() command for the script to run.

#in_file.close()
