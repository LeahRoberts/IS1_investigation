Finding location of IS reads among contigs of ST131 strains
=============================================================

How to map ST131 reads back to their contigs, parse the reads that are known to be related to IS1, and count the number of occurrences:
-----------------------------------------------------------------------------------------------------------------------------------------

*This was a preliminary method to the now finalised "IS_Locator.sh" method. As such, this README should be used for reference only.*

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

Creates a file called locations.txt that has all of the NODE_* information regarding where the reads map to. This looks like::

	NODE_601_length_177_cov_56.146893
	NODE_601_length_177_cov_56.146893
	NODE_601_length_177_cov_56.146893 

Next, using the command::

	cat locations.txt | sort | uniq -c

Sorts all of the NODE_* locations, and counts the number of occurrences.
This looks like::

	12 *
	2 NODE_108_length_19307_cov_9.113430
	2 NODE_112_length_20804_cov_8.447606
	2 NODE_122_length_3609_cov_9.052369
	70 NODE_128_length_6478_cov_8.512504
	2 NODE_18_length_26204_cov_8.977714
	2 NODE_1_length_30927_cov_9.063634
	34 NODE_225_length_1306_cov_8.650077
	6 NODE_239_length_8961_cov_9.351747
	2 NODE_241_length_7970_cov_9.011417

Where the number of the left is the number of occurrences of reads mapping to the contig on the right.

This can be parsed into a separate file::

	cat locations.txt | sort | uniq -c > location_occur.txt

Problem with this approach:
----------------------------

When it came to the analysis, I realised a problem with this approach - whilst I knew which contigs were involved, I had no way of knowing where on the contig the reads had mapped to. I could look at the bam file again to find these locations, but that would be very labour intensive and take too long. 

Instead, I decided to go back to the fastq parser approach. By parsing only the reads I want and mapping those back to the contigs, I can open the bam file in artemis to visualise where the reads are mapping. 
