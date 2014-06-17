#!/bin/bash

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
        fi
done



