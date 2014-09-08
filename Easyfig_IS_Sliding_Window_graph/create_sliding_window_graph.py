import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

in_files = glob.glob('Genomes/*.embl')
window_size = 100000

print_counter = 0

IS_count = 0

for f in in_files:
	
	cur_genome = SeqIO.parse(f, "embl")
	for record in cur_genome:
		print "infile = %s" % (record.description)	
		sequence_length = len(record)
		start_pos = 1
		end_of_line = len(record) - window_size
		if start_pos >= end_of_line:
			end_pos = len(record)
		else:
			end_pos = start_pos + window_size
		print "starting analysis..."
		
		while start_pos <= sequence_length: 
			for feat in record.features:
				if feat.type == "operon":
					if feat.location.end in range(start_pos, end_pos):
						#print "found feature in range"
						IS_count += 1
						#print IS_count
			
			newfile = open("IS_counts.txt", "a+")
			
			newfile.write("%i, %i, %i\n" % (start_pos, end_pos, IS_count))
			
			IS_count = 0
			start_pos += 100
			end_pos += 100
			print_counter += 100
			
			if print_counter == 500000:
				print "500,000 positions have passed..."
				print_counter = 0
			
newfile.close()
print "finished!"				
