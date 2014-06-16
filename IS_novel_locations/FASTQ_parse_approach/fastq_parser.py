#!/usr/bin/env python

# Fastq parser
# Parses out reads that match those in a txt file

import glob
import os
from Bio import SeqIO

# Make a list of the read names from a txt file:

with open('list.txt', 'r') as f:
        read_ids = [line.strip() for line in f]
        print len(read_ids)
#	print read_ids

# Takes the fastq read ID and compares it to the read names in the read_ids list
# If the read is found in the read_ids list, it is written out to a "new_reads.fastq" file

	in_file = raw_input("fastq file: ")
        read = SeqIO.parse(in_file, "fastq")
        for fastq_rec in read:
#		print fastq_rec.id
		if fastq_rec.id in read_ids:
			print fastq_rec.id + " exists!"
			exit()
			output = open("new_reads.fastq", "a")
			SeqIO.write(fastq_rec, output, "fastq")
		
		#	if fastq_rec in read_ids:
                #	print "found one read"
                #	newreads.append(read)
		#	SeqIO.write(fastq_rec, output_file, "fastq")
