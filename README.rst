IS1 Investigation of the ST131 99 lineage
==========================================

Overview:
---------

This script is designed to map Illumina paired-end reads (using bwa) to an IS reference sequence (without flanking regions), parse out reads that map to the reference and de novo assemble them. The theory is that reads that overlap the edges of the IS will de novo assemble into the flanking regions which will give the contextual location of the IS.  

Nucmer was used to compare the contigs to firstly the EC958 IS (with 100 bp flanking region) and then to the EC958 chromosome (to find novel insertion sites - as described below). 

Method:
-------

Several files need to be arranged accordingly for the script to work.

1. The fastq read files (not interleaved) need to be in the one directory and need to be named: strain_1.fastq
2. A list of the kmers to be used for each strain needs to be given in the directory above the reads
3. A list of the insert size (as well as standard deviation) also needs to be given in the directory above the reads
4. The IS reference should be .fasta and in the directory above
5. The comparison file (containing IS1 and flanking regions-100 bp- from EC958) should be two directories above (or change the script)

The output will give:

1. The usual bwa output
2. A strain.velvet.bam file which contains only the reads that mapped to the IS reference.
3. A folder (same name as the strain) which has the velvet assembly data
4. A contigs.fa file that contains the assembled contigs
5. An out.delta file that contains the nucmer results of the comparison of the assembled contigs to the EC958 IS1 plus flanking region
6. A log.txt file that contains a summary of how many contigs were made, and how many match the IS1 and flanking regions in EC958 (NOTE: due to repeats in the out.delta file, number of matching contigs may not be 100% correct)
7. **NEW** A txt file with the names of the reads from the velvet.bam file

You will likely also need to know the length of your query EC958 IS to determine whether there are matches at the flanking regions (python script included). 

Manually investigating each of these will tell you whether the flanking regions match EC958 or not.

For contigs that don't match EC958, they were manually parsed into a new file (rem_contigs.fa) and compared to the EC958 chromosome using nucmer (rem_contigs.delta), on the assumption that these contigs were potentially novel flanking regions, and that the ST131 strains were highly similar to EC958 (potentially valid for clade C, probably not as much for clade A and B). 
