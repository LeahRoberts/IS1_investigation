#!/bin/bash

# Author: Leah Roberts
# Affiliation: Scott Beatson Lab Group (University of Queensland St Lucia)
# Date: June 2014

##################### IS1_Locator #####################

### Description:

# Script to map reads to a reference using bwa, filter for reads mapping to the reference and de novo assemble 
# those reads using velvet

### Installation Requirements:

# Need to have bwa, velvet and samtools installed 

### File Formats: 

# Reads are paired end but no interleaved and are in the format name_1.fastq, name_2.fastq.
# Reference is fasta format

### Run the script:

# The script needs to be executed in the directory which has the read files
# the reference needs to be one directory above

# Execute the script by typing:
# $ bash IS_locator.sh $REFERENCE '$REFERENCE_HEADER'

########################################################


REFERENCE=$1
IS=$2 # This is the header of the REFERENCE fasta file and is required to filter for reads that map to the IS

echo "indexing " $REFERENCE
bwa index $REFERENCE

# Write a loop that will generate .sai files by reading each of the paired strains to the 
# reference strain using the tool BWA. These are outputted as $name_1.sai or $name_2.sai.

for f in * 

do

# Used cut to parse out the strain name - THIS WILL ONLY WORK IF THE FASTQ FILE IS 
# FORMATTED CORRECTLY (i.e. $name_1.fastq).
# Prints the name of the strain currently being processed:

	echo "processing $(echo $f | cut -f1 -d.)"
	
# Formatted script so that it will take in both zipped and unzipped fastq files	

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
						echo "Filtering reads that map to the reference"
						samtools view $name.sorted.bam $IS > $name.sorted.mapped.sam
						
# Removes the original fastq files. These should be commented out if using the script for the first time or 
# if the user prefers to keep all the data:                        			
						rm $name.sam

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
 						echo "Filtering reads that map to the reference"
						samtools view $name.sorted.bam $IS > $name.sorted.mapped.sam
						rm $name.sam
                        		fi
                		done
            		fi
		fi
	fi
	

done


# Cleanup the files you don't need:
mkdir tmp
for f in *
do
        if [[ $f == *.fastq ]]
        then
                mv $f tmp/

        elif [[ $f != *.sorted.mapped.sam ]]
        then
                rm $f
        fi
        mv tmp/* ./
        rm -r tmp/
done


for f in *
do
	if [[ $f == *mapped.sam ]]
	then
		name=$(ls $f | cut -f1 -d.)
		echo "indexing fasta file to create .fai file"
		samtools faidx $REFERENCE
		echo "creating bam file from sam file"
		samtools view -bt $REFERENCE.fai $name.sorted.mapped.sam > $name.velvet.bam
		echo "finished creating bam file for velvet"
		rm $name.sorted.mapped.sam
	fi
done
	
# Can then use this $name.velvet.bam file in the velvet assemblies

# Run velvet assembly using the .bam file generated from just the reads mapping to the chromosome

contigs=../contigs.txt			# This is a list of the khmer size to be used in the velvet assembly
list=../list.txt			# This is a list of insert sizes and insert size standard deviation to be use in the velvet assembly
IS_reference=$3	# This is the IS reference file that contains the IS sequence plus 100 bp of flanking sequence either side. This reference will be compared against the newly generated contigs (from the IS reads assembly) to determine if the IS1 flanking regions are the same, or different. 

# The section below assembled the IS read files using velvet:

for f in *
do
	if [[ $f == *velvet.bam ]]
	then
		name=$(ls $f | cut -f1 -d.)
		echo $name

		kmer=$(cat $contigs | grep $name | cut -f2 -d_)
#		echo $kmer
		insert=$(cat $list | grep $name | cut -f2 -d,)
#		echo $insert
		sd=$(cat $list | grep $name | cut -f3 -d,)
#		echo $sd
		
		echo "running velveth on" $name
		velveth $name/ $kmer -bam -shortPaired $name.velvet.bam
		echo "creating contigs"
		velvetg $name/ -exp_cov auto -cov_cutoff auto -ins_length $insert -ins_length_sd $sd
		echo "finished running assembly"
		
# The next section runs a Nucmer comparison between the newly generated contigs and the $IS_reference, which contains
# the same IS from EC958 plus 100 bp of flanking region either side.
		
		cd $name/
		contig_num=$(cat contigs.fa | grep -c ">")
		nucmer $IS_reference contigs.fa
		nuc_results=$(cat out.delta | grep -c ">")
		echo $name "has" $contig_num "contigs" > log.txt
		echo $name "has" $nuc_results "matching contigs out of" $contig_num >> log.txt
		cd ../
		mv $name* $name/
		cd $name
		
# The next section requires the addition of two other script
# This section will parse the reads that map to the IS, and re-map them to the draft contig from the same strains 
# as the original read files:

		for f in *
		do
			if [[ $f == *.velvet.bam ]]
			then
				samtools view $f | cut -f1 -d$'\t' > list.txt
				cat list.txt > list_2.txt
				sed -i 's/$/\/1/g' list.txt
				sed -i 's/$/\/2/g' list_2.txt
				python ../../fastq_parser.py $name\_1.fastq $name\_2.fastq
				echo "renaming output fastq files"
				mv $name\_reads_1.fastq $name\_1.fastq
				mv $name\_reads_2.fastq $name\_2.fastq
				bash ../../align_to_contigs.sh
			fi
		done
		cd ../
	fi	
done


 
