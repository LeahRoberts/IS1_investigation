IS1 Investigation of the ST131 99 lineage
==========================================

Overview:
---------

This script has 2 outcomes:

1. Mapping of reads to IS in order to resolve flanking regions
2. Parsing reads that match IS and mapping them back to their original contigs


For the first part, this script is designed to map Illumina paired-end reads (using bwa) to an IS reference sequence (without flanking regions), parse out reads that map to the reference and de novo assemble them. The theory is that reads that overlap the edges of the IS will de novo assemble into the flanking regions which will give the contextual location of the IS.  

Nucmer is then used to compare the contigs to firstly the EC958 IS (with 100 bp flanking region) and then manually to the EC958 chromosome (to find novel insertion sites - as described below). 

For the second part, the reads that are found to map to the IS in part 1 are parsed into separate fastq files and mapped back to their original assembled contigs. This should show where the reads are mapping to and give an indication of where IS are in the original assembly.

Setup:
-------

Several files need to be arranged accordingly for the script to work.

1. The fastq read files (not interleaved) need to be in the one directory and need to be named: strain_1.fastq
2. A list of the kmers to be used for each strain needs to be given in the directory above the reads
3. A list of the insert size (as well as standard deviation) also needs to be given in the directory above the reads
4. The IS reference should be .fasta and in the directory above
5. The comparison file (containing IS1 and flanking regions-100 bp- from EC958) should be in the same directory as the kmer list etc.
6. The "fastq_parser.py" script should also be in the directory above the raw reads
7. The "align_to_contigs.sh" script should also be in the directory above the raw reads

Run the Script:
----------------

Once you've completed the above requirements and you're in the directory containing all the raw illumina reads, to run the script simply type::

  bash ~/bin/bwa_velvet_assembly.sh
  
**NOTE** - All of the files mentioned in the Setup section need to be named according to their description in the bash script

Output:
--------

**NOTE** - To save space, much of the generated data is deleted, including the original bam, sam, and sai files, as well as the original fastq files.

1. A strain.velvet.bam file which contains only the reads that mapped to the IS reference.
2. A folder (same name as the strain) which has the velvet assembly data
3. A contigs.fa file that contains the assembled contigs
4. An out.delta file that contains the nucmer results of the comparison of the assembled contigs to the EC958 IS1 plus flanking region*
5. A log.txt file that contains a summary of how many contigs were made, and how many match the IS1 and flanking regions in EC958 (NOTE: due to repeats in the out.delta file, number of matching contigs may not be 100% correct)
6. Two text file with the names of the reads from the velvet.bam file - these are modified to either match the strain_1.fastq file or the strain_2.fastq file. 
7. Two new strain_1.fastq files - these are not the original fastq files - these only contain the reads that mapped to the reference IS
8. Bam, sam etc files - the usual output from bwa - which relate to the IS reads that have been mapped back to the strain's original contigs.

*You will likely also need to know the length of your query EC958 IS to determine whether there are matches at the flanking regions (python script included). 

Manually investigating each of these will tell you whether the flanking regions match EC958 or not.

For contigs that don't match EC958, they were manually parsed into a new file (rem_contigs.fa) and compared to the EC958 chromosome using nucmer (rem_contigs.delta), on the assumption that these contigs were potentially novel flanking regions, and that the ST131 strains were highly similar to EC958 (potentially valid for clade C, probably not as much for clade A and B). 

As for the .bam output, visualising this in Artemis against the ST131 contigs will show peaks where the reads are mapping. In the case of IS1, this won't show all IS1 due to the bias in using different IS1 references. 
