#!/bin/bash

# Script to map reads to a reference using bwa, filter for reads mapping to the reference and de novo assemble those reads using velvet
# Reads are paired end but no interleaved
# reference is fasta format
# Results in bam file output of filtered reads that can be used in velvet; contigs of de novo assembly

# Need to have bwa, velvet and samtools installed 
# need to execute in the directory which has the reads
# the reference needs to be one directory above
# uncomment the bwa index step if the reference has not been indexed before
# The reads need to be in the format strainname_1.fastq for the script to work properly

REFERENCE=../IS1_yehA.fasta
IS='IS1_2384715..2385482'

bwa index $REFERENCE


for f in *
do

# Used cut to parse out the strain name - THIS WILL ONLY WORK IF THE FASTQ FILE IS FORMATTED CORRECTLY (i.e. $strainname_1.fastq) 
        echo "processing $(ls $f | cut -f1 -d.)"

                if [[ $f == *_1.fastq ]]
                then
			read1=$f
                        name1=$(ls $f | cut -f1 -d.)
                        bwa aln $REFERENCE $read1 > $name1.sai
#                       echo $name1
#                       echo $read1
                elif [[ $f == *_2.fastq ]]
		then
                        read2=$f
                        name2=$(ls $f | cut -f1 -d.)
#                        echo "paired-end read " $name2
                        bwa aln $REFERENCE $read2 > $name2.sai
#                       echo $name2
#                       echo $read2
        fi
done

# Write another loop that will take the .sai files from the previous section, and map them to the reference strain generating .sam and .bam files.

for f in *
do
	if [[ $f == *fastq ]]
	then
#  Parse out just the name of the strain (again, in the format $strainname_1.fastq) 
        name=$(ls $f | cut -f1 -d_)
        echo $name       
	        
		if [[ ! -e $name.bam ]]
        	then
         	echo "performing alignment on " $name

# Make sure that the names of paired files (.sai and .fastq) are the SAME and that the SAME strain files are being aligned to the reference (again using BWA):

                	if [[ $(ls $f | cut -f1 -d.) == *_1* ]]
                	then
                        	name1=$(ls $f | cut -f1 -d.)
#                       	echo $name1

                	elif [[ $(ls $f | cut -f1 -d.) == *_2* ]]
                	then
                        	name2=$(ls $f | cut -f1 -d.)
#                       	echo $name2             

			 	if [[ $(ls $name1 | cut -f1 -d_) == $(ls $name2 | cut -f1 -d_) ]]
                        	then
                                	echo "creating sam file"
                                	bwa sampe $REFERENCE $name1.sai $name2.sai $name1.fastq $name2.fastq > $name.sam
					echo "finished creating .sam file"
					samtools view -bS $name.sam > $name.bam
					samtools sort $name.bam $name.sorted
					samtools index $name.sorted.bam
					echo "Filtering reads that map to the reference"
					samtools view $name.sorted.bam $IS > $name.sorted.mapped.sam
				fi
			fi
		fi
	fi
done

# Cleanup the files you don't need (i.e. everything but the velvet.bam file):
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
	
# Can then use this $name.header.bam file in the velvet assemblies

# Run velvet assembly using the .bam file generated from just the reads mapping to the chromosome
contigs=../contigs.txt
list=../list.txt
IS_reference=../../EC958_IS1.fa

for f in *
do
	if [[ $f == *velvet.bam ]]
	then
		name=$(ls $f | cut -f1 -d.)
		echo $name

		kmer=$(cat $contigs | grep $name | cut -f2 -d_)
		echo $kmer
		insert=$(cat $list | grep $name | cut -f2 -d,)
		echo $insert
		sd=$(cat $list | grep $name | cut -f3 -d,)
		echo $sd
		
		echo "running velveth on" $name
		velveth $name/ $kmer -bam -shortPaired $name.velvet.bam
		echo "creating contigs"
		velvetg $name/ -exp_cov auto -cov_cutoff auto -ins_length $insert -ins_length_sd $sd
		echo "finished running assembly"
		
		cd $name/
		contig_num=$(cat contigs.fa | grep -c ">")
		nucmer $IS_reference contigs.fa
		nuc_results=$(cat out.delta | grep -c ">")
		echo $name "has" $contig_num "contigs" > log.txt
		echo $name "has" $nuc_results "matching contigs out of" $contig_num >> log.txt
		cd ../
		mv $name* $name/
		cd $name
		
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


 
