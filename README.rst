IS Investigation in the ST131 Lineage
=======================================

Several scripts were generated to fulfill the purpose of finding IS similar to EC958, and to also find novel IS1 locations in the ST131 lineage.

Overview of the IS_Locator script:
-----------------------------------

This script has 2 outcomes:

1. Mapping of reads to an IS reference in order to resolve flanking regions
2. Parsing reads that map to an IS reference and mapping them back to the draft genome contigs for that particular strain


For the first part, this script is designed to map Illumina paired-end reads (using bwa) to an IS reference sequence (without flanking regions), parse out reads that map to the reference and de novo assemble them. The theory is that reads that overlap the edges of the IS will de novo assemble into the flanking regions which will give the contextual location of the IS.  

Nucmer is then used to compare the contigs to firstly the EC958 IS (with 100 bp flanking region) and then manually to the EC958 chromosome (to find novel insertion sites - as described below). 

For the second part, the reads that are found to map to the IS in part 1 are parsed into separate fastq files and mapped back to their original assembled contigs. This should show where the reads are mapping to and give an indication of where IS are in the original assembly.

Installation Requirements
---------------------------

1. Burrows-Wheeler Aligner (bwa) (http://bio-bwa.sourceforge.net/)
2. SAMtools (http://samtools.sourceforge.net/)
3. Python (2.7 plus BioPython modules)
4. Velvet (https://www.ebi.ac.uk/~zerbino/velvet/)
5. Nucmer (http://mummer.sourceforge.net/)

Setup - Brief Overview:
------------------------

Several files need to be arranged accordingly for the script to work.

1. The fastq read files (not interleaved) need to be in the one directory and need to be named: strain_1.fastq
2. A list of the kmers to be used for each strain needs to be given in the directory above the reads, named **contigs.txt**
3. A list of the insert sizes (as well as standard deviation) also needs to be given in the directory above the reads, named **list.txt**
4. The IS reference should be .fasta and in the directory above
5. The comparison file (containing IS1 and flanking regions-100 bp- from EC958) should be in the same directory as the kmer list etc.
6. The "fastq_parser.py" script should also be in the directory above the raw reads
7. The "align_to_contigs.sh" script should also be in the directory above the raw reads

Run the Script:
----------------

Once you've completed the above requirements and you're in the directory containing all the raw illumina reads, to run the script simply type::

  bash ~/bin/bwa_velvet_assembly.sh ../$IS_REFERENCE_FOR_MAPPING.fa '$IS_reference_header' ../../IS_REFERENCE_FOR_COMPARISON.fa
  

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

As for the .bam output, visualising this in Artemis against the ST131 contigs will show peaks where the reads are mapping. 

Example File formats:
---------------------

1. contigs.txt file::

  B36EC_81_Contigs.fasta
  HVM1147_73_Contigs.fasta

The middle number represents the khmer length.

2. list.txt file::

  B36EC,248.06,93.39
  HVM1147,239.34,89.03
  
The first number is the insert size, and the second number is the SD of the insert size.

3. The IS reference for the first mapping step::

  >IS1_EC958
  GGTGATGCTGCCAACTTACTGATTTAGTGTATGATGGTGTTTTTGAGGTGCTCCAGTGGC
  TTCTGTTTCTATCAGCTGTCCCTCCTGTTCAGCTACTGACGGGGTGGTGCGTAACGGCAA
  
The header should also be included in the command as 'IS1_EC958' (including the '').

4. IS1 from EC958 reference file::

  >IS1,IS1..3288120..3288887(1),100bp flanked,[EC958 IS]
  CGGAAGAATCAGAGGCTGTGGTTTCAGACTGTCTGCCAGTACATTCCTCTCTCCGTTAAAAACCATAACGGGTTCATTATCTTCGTCTGTCAGCAGATTGGGTGATGCTGCCAACTTACTGATTT   AGTGTATGATGGTGTTTTTGAGGTGCTCCAGTGGCTTCTGTTTCTATCAGCTGTCCCTCCTGTTCAGCTACTGACGGGGTGGTGCGTAACGGCAAAAGCACCGCCGGACATCAGCGCTATCTCTG   CTCT
  >IS1,IS1..3290147..3290914(-1),100bp flanked,[EC958 IS]
  GAAAGATGGTGATAATGTGCTGCATTATACTGCGATTGTTAAGAAGTCGTCAGCCAATAATGCCCAAGTCACTGAGGGTGCTTTTTCTGCAGTCGCAACCGGTGATGCTGCCAACTTACTGATTT   AGTGTATGATGGTGTTTTTGAGGTGCTCCAGTGGCTTCTGTTTCTATCAGCTGTCCCTCCTGTTCAGCTACTGACGGGGTGGTGCGTAACGGCAAAAGCACCGCCGGACATCAGCGCTATCTCTG   CTCT
