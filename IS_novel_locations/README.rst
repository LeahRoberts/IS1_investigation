Finding location of IS reads among contigs of ST131 strains
=============================================================

How to map ST131 reads back to their contigs, parse the reads that are known to be related to IS1, and count the number of occurrences:
-----------------------------------------------------------------------------------------

**NOTE** - There is a problem with this approach in the analysis - It is very difficult see where the reads are actually mapping. It would be easier to look at a bam pileup in Artemis to see where IS reads are and where their paired read is mapping. Please see "FASTQ_parse_approach" folder for alternate method.

Step 1:
--------

Firstly, you need to know what reads map to the IS element.
We have already mapped reads back to IS and filtered this to give a bam file of only reads mapped to IS.
By parsing this bam file we can get a list of the read names that match IS::

	samtools view B36EC.velvet.bam | cut -f1 -d$'\t' > B36EC.readnames.txt



Step 2:
--------

Mapping the reads back to their contigs means that you need to have a script that will use multiple references.

The "align_to_contigs.sh" script will do this provided that the necessary information is in the correct location.

This script is in ~/bin on binf-training.


Step 3:
---------

We now need to parse the reads from the .bam file generated in step 2 that match the list of read names we generated in step 1.
This can take a while - B36EC (the test strain) had 2.5 million reads, of which roughly 5000 were parsed.::

	cat B36EC.readnames.txt | while read name; do samtools view B36EC.sorted.bam | grep $name >> IS1_related_reads.sam; done



Step 4:
--------

Once the .sam file has been created with all the appropriate reads, using the command::

	cat $IS1_related_reads.sam | cut -f3 -d$'\t' > locations.txt

Creates a file called locations.txt that has all of the NODE_* information regarding where the reads map to. 

Next, using the command::

	cat locations.txt | sort | uniq -c

Sorts all of the NODE_* locations, and counts the number of occurrences.
This can be parsed into a separate file::

	cat locations.txt | sort | uniq -c > location_occur.txt
