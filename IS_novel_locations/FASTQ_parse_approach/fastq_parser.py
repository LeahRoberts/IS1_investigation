#!/usr/bin/env python

# Fastq parser
# Parses out reads that match those in a txt file

# USAGE: 

# Simple overview:

# Execute in the directory with:
# -- list.txt
# -- list_2.txt
# -- strain_1.fastq
# -- strain_2.fastq

# Command: python fastq_parser.py strain_1.fastq strain_2.fastq

# Output: strain_reads_1.fastq (the parsed version of strain_1.fastq) and strain_reads_2.fastq (the parsed
# version of strain_2.fastq). 

# Script Requirements/Explanation:

# This script requires the generation of a list of reads names and the slight alteration of those list names. 
# The list of read names (i.e. the reads that you want parsed) can be generated from a bam file using samtools:
# $ samtools view <file.bam> | cut -f1 -d$'\t' | sort | uniq >> list.txt

# One minor complication is that the read names in the fastq files are followed by a '/1' or a '/2' depending
# on whether it is 'read 1' or 'read 2'. Depending on how the read list is parsed initially, the reads may or 
# may not have this ending.
# To account for this, two list files need to be generated and modified so that the read names have either a
# '/1' or '/2' at the end. This can very simply be done by using a sed command.

# Ultimately, you should have:
# -- Two list files, named list.txt and list_2.txt. These should contain the same read names, however, one should
# have '/1' appended to the end of all the reads, while the other should have '/2' appended to the end of the reads.
# -- the fastq files, which are entered as arguments in the command (see above). 

### NOTE: This script will delete the original fastq files when it finishes parsing them. For good measure, it is
# advisable to disable this feature until you're sure that the reads have parsed properly. To disable this, 
# comment out the two line that read "os.remove(in_file)".  

import glob
import os
from Bio import SeqIO
import sys

# Make a list of the desired read names from a txt file:

with open('list.txt', 'r') as f:
        read_ids = [line.strip() for line in f]
        print "%d reads imported" % (len(read_ids))  
#       print len(read_ids)
#	print read_ids

# Takes the fastq read ID and compares it to the read names in the read_ids list
# If the read is found in the read_ids list, it is written out to a "new_reads.fastq" file

	in_file = sys.argv[1]
	name = in_file.split("_")[0]
#       print name
	print "Parsing " + in_file
        reads = SeqIO.parse(in_file, "fastq")
        output = open(name + "_reads_1.fastq", "w")
        for fastq_rec in reads:
#		print fastq_rec.id
		if fastq_rec.id in read_ids:
#			print fastq_rec.id + " exists!"
#			exit()
			SeqIO.write(fastq_rec, output, "fastq")
		
print "Finished parsing " + in_file
os.remove(in_file)

# Next loop does the same thing except with the second list and the second fastq file

with open('list_2.txt', 'r') as f:
        read2_ids = [line.strip() for line in f]
	print "%d reads imported" % (len(read_ids))
#       print read2_ids

        in_file = sys.argv[2]
        name2 = in_file.split("_")[0]
        print "Parsing " + in_file
#	print name2
        reads = SeqIO.parse(in_file, "fastq")
        output = open( name2 + "_reads_2.fastq", "w")
        for fastq_rec in reads:
#               print fastq_rec.id
                if fastq_rec.id in read2_ids:
#                	print fastq_rec.id + " exists!"
#                       exit()
                        SeqIO.write(fastq_rec, output, "fastq")

print "Finished parsing " + in_file
os.remove(in_file)
