Fastq Parser
=============

This script parses predefined reads from a fastq file and puts them in a new fastq file. 

Step 1: Generate a list of read names that you want to parse
---------------------------------------------------------------
**Example lists can be found in the Example_files folder**

The lists is generated from the velvet.bam file which has filtered reads that map to the IS reference::

	cat $name.velvet.bam | cut -f3 -d$'\t' > list.txt

This list then needs to be modified to match the fastq id::

	sed -i 's/$/\/1/g'
	
**NOTE** - The read id does not need to have an '@' symbol at the beginning.
This list can be used to parse reads from the $name_1.fastq file.

For the $name_2.fastq files, create another list idential to the previous list, except with a "/2" at the end::

	sed -i 's/\/1/\/2/g' 
	
Use these lists when parsing the fastq files.

These lists should be named **list.txt** and **list_2.txt**.

Step 2: How-to: Run Fastq Parser
---------------------------------

Run the fastq_parser.py script as follows in the directory with the fastq reads::

	python fastq_parser.py strain_1.fastq strain_2.fastq

The two lists containing the names of the reads you want to parse need to be in the same directory as the reads ([] refers to which argument they relate to)::

	list.txt [1]
	list_2.txt [2]

This will give the output::

	strain_reads_1.fastq
	strain_reads_2.fastq

These new fastq files contain only the reads stipulated in the pre-made lists. 

**Example output can be found in the Example_files folder**

*Note: The original fastq files will be deleted - this can be prevented by commenting out that line in the script.*

Step 3: Map these reads back to the contigs
---------------------------------------------

See the "align_to_contigs.sh" script. 
