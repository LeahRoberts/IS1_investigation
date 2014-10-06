#!/bin/env python

# Author: Leah Roberts
# Affiliation: Scott Beatson Lab Group (University of Queensland St Lucia)
# Date: September 2014

################ Sliding Window Script ###########################
########### Used to created a graph for EasyFig ##################

# This script is designed to iteratively count the number of genes 
# within a predefined window size. The output is a csv file, which 
# gives the start and end location of the window, as well as the
# number of specified genes within that region.
# This csv file can also be easily parsed into a single text file 
# suitable as input into EasyFig. Input into Easyfig will result 
# in a graph above the genome figure corresponding to the number 
# of genes in each window across the genome. 

# To parse the results to a text file for EasyFig, type:
### cat results.txt | cut -f3 -d, > results_easyfig.txt

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation


# The infile need to be in a directory named "Genomes"
# The infile needs to be in embl format

in_files = glob.glob('Genomes/*.embl')	

# Defines the window size that the script will move through
# This can be changed based on preference, but will make
# the script run slower.
		
window_size = 100000

print_counter = 0

result = 0

for f in in_files:
	
	cur_genome = SeqIO.parse(f, "embl")			# Parses the input file
	for record in cur_genome:
		print "infile = %s" % (record.description)	# Prints the name of the infile (can also be "record.id")
		sequence_length = len(record)
		start_pos = 1					# Defines the starting position
		end_of_line = len(record) - window_size		# Defines the end of the "sliding window"	
	
		# Define the end position:
		# The end position is either:
			# The start position + the window size; or
			# The very end of the sequence (if the start position moves 
			# within the last sliding window region)

		if start_pos >= end_of_line:		
			end_pos = len(record)
		else:
			end_pos = start_pos + window_size
		print "starting analysis..."
		
		while start_pos <= sequence_length: 
			for feat in record.features:
				if feat.type == "operon":	# Defines which gene type to count - this can be changed
					if feat.location.end in range(start_pos, end_pos):
						result += 1
		
			newfile = open("results.txt", "a+")  	# Opens a new file to write results 
			
			newfile.write("%i, %i, %i\n" % (start_pos, end_pos, result))	# Writes out results
			
			IS_count = 0				# Resets the counter for the new region
			start_pos += 100			# Moves to the next region in predefined increments - these can be changed
			end_pos += 100				# Moves the end position - this needs to be the same as start_pos
			print_counter += 100			# This counts the number of positions that have been iterated through
			
			# The 'if' statement below is used to check that the script is still working
			# As described above, the print_counter variable tracks the number of positions
			# that have been iterated through.
			# Once the print counter reaches 500000 (or any predefined number), it prints 
			# a working statement.
			# The counter then resets to zero, and the counting starts over again. 
			
			if print_counter == 500000:
				print "500,000 positions have passed..."
				print_counter = 0
			
newfile.close()
print "finished!"				
