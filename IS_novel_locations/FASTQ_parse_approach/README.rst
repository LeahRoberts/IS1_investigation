Designing a FASTQ parser using Biopython
==========================================
FASTQ_parse_approach:

Step 1: Design a python script to parse fastq records
------------------------------------------------------

Designed a fastq_parser (python script) which parses out reads that match those in a specified list.

Step 2: Generate a list of read names that you want to parse
---------------------------------------------------------------

The lists is generated from the velvet.bam file which has filtered reads that map to the IS reference::

	cat $name.velvet.bam | cut -f3 -d$'\t' > list.txt

This list then needs to be modified to match the fastq id::

	sed -i 's/^/@/g'; sed -i 's/$/\/1/g'
	
This list can be used for the $name_1.fastq file.
For the $name_2.fastq files, change the previous list to have "/2" at the end::

	sed -i 's/\/1/\/2/g' 
	
Use these lists when parsing the fastq files.


Step 3: Parse fastq file
--------------------------

Not working atm


Step 4: Map these reads back to the contigs
---------------------------------------------

I already have a script to do this: "align_to_contigs.sh"

Just need to change the fastq file being used to align.
