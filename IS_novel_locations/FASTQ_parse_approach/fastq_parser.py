#!/usr/bin/env python

# Fastq parser
# Parses out reads that match those in a txt file

import glob
import os
from Bio import SeqIO
import sys

# Make a list of the desired read names from a txt file:

with open('list.txt', 'r') as f:
        read_ids = [line.strip() for line in f]
        print len(read_ids) + " read ids imported"
#       print len(read_ids)
#	print read_ids

# Takes the fastq read ID and compares it to the read names in the read_ids list
# If the read is found in the read_ids list, it is written out to a "new_reads.fastq" file

	in_file = sys.argv[1]
	print "Parsing " + in_file
        reads = SeqIO.parse(in_file, "fastq")
        output = open("new_reads.fastq", "w")
        for fastq_rec in reads:
#		print fastq_rec.id
		if fastq_rec.id in read_ids:
#			print fastq_rec.id + " exists!"
#			exit()
			SeqIO.write(fastq_rec, output, "fastq")
		
print "Finished parsing"	
