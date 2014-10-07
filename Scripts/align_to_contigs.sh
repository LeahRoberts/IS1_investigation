#!/bin/bash

# Author: Leah Roberts
# Affiliation: Scott Beatson Lab Group (University of Queensland St Lucia)
# Date: June 2014

#################### Align to Contigs Script ##########################

### This script is to be used in conjunction with the IS_Locator.sh script.
# It is designed to map the reads from 99 ST131 strains back to their velvet contigs using BWA

for f in *
do
	if [[ $f == *.fastq ]]
	then

		name=$(ls $f | cut -f1 -d_)
		echo "processing " $name
	
		REFERENCE=../../ST131_99_ord/$name\_[[:digit:]]*\_Contigs.fas
		echo "Reference is " $REFERENCE
	
		bwa index $REFERENCE
	
		if [[ $f == *_1.fastq.gz ]]
		then
			gunzip $f
			read1=$(echo $f | cut -f1-2 -d.)
			name1=$(echo $f | cut -f1 -d.)
			bwa aln $REFERENCE $read1 > $name1.sai
		
		elif [[ $f == *_2.fastq.gz ]]
		then
			gunzip $f
			read2=$(echo $f | cut -f1-2 -d.)
			name2=$(echo $f | cut -f1 -d.)	
			bwa aln $REFERENCE $read2 > $name2.sai
		
		elif [[ $f == *_1.fastq ]]
		then
			name1=$(echo $f | cut -f1 -d.)
                        bwa aln $REFERENCE $f > $name1.sai
		
		elif [[ $f == *_2.fastq ]]
		then	
			name2=$(echo $f | cut -f1 -d.)
                        bwa aln $REFERENCE $f > $name2.sai
	fi
done

# Write another loop that will take the .sai files from the previous section, and map them 
# to the reference strain generating .sam and .bam files.

for f in *
do
	if [[ $f == *fastq ]]
	then

# Parse out just the name of the strain (again, in the format $name_1.fastq). 
		name=$(echo $f | cut -f1 -d_)
		
# Want to only perform the alignment once - as there are two .sai files, this has the 
# potential to interate through twice. This if statement prevents the script from 
# iterating through more than once on the same strain:
		if [[ ! -e $name.bam ]]
		then
			echo "checking for pairs - " $name

            		if [[ $(echo $f | cut -f1 -d. | cut -f2 -d_) == "1" ]]
            		then
            			name1=$f
            			echo "first paired read = " $name1
            			
# Make sure that the names of paired files (.sai and .fastq) are the SAME and that the SAME strain files are being aligned 
# to the reference (again using BWA):
                		for g in *
                		do
                			if [[ $(echo $g | cut -f2 -d_) == "2.fastq" ]] && [[ $(echo $g | cut -f1 -d_) == $name ]]
                    			then
                    				echo $f "and" $g "are a pair - performing alignment"
                        			bwa sampe $REFERENCE $name\_1.sai $name\_2.sai $f $g > $name.sam
                        			samtools view -bS $name.sam > $name.bam
                        			samtools sort $name.bam $name.sorted
                        			samtools index $name.sorted.bam
						rm $name.bam $name.sam
						echo "finished read-mapping for " $name
                    			fi
                		done

# A second loop exactly the same as the first, except it takes in "read2" files:                		

			elif [[ $(echo $f | cut -f1 -d. | cut -f2 -d_) == "2" ]]
         		then
                    		name2=$f
                    		echo "second paired read = " $name2

                    		for g in *
                    		do
                    			if [[ $(echo $g | cut -f2 -d_) == "1.fastq" ]] && [[ $(echo $g | cut -f1 -d_) == $name ]]
                        		then
                        			echo $f "and" $g "are a pair"
                            			bwa sampe $REFERENCE $name\_1.sai $name\_2.sai $g $f > $name.sam
                            			samtools view -bS $name.sam > $name.bam
                            			samtools sort $name.bam $name.sorted
                            			samtools index $name.sorted.bam
 						rm $name.bam $name.sam
						echo "finished read-mapping for " $name
                        		fi
                		done
            		fi
		fi
	fi
	

done
