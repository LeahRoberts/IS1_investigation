Example Files for use in IS1 investigation
==========================================

Used when mapping reads to a reference IS and de novo assembling those reads.

List of Files:
---------------

1. contigs.txt
2. EC958_complete.fasta
3. EC958_IS1.fasta
4. IS1.fasta
5. list.txt

Description:
-------------

1. This file has the strain name and the kmer size which is used when using velvet to assemble the reads. It looks like this::

	B36EC_81_Contigs.fasta
	HVM1147_73_Contigs.fasta

Where the number indicates the kmer size.

2. This is the complete EC958 chromosome in fasta format used to test the "unknown" contigs which have assembled but have no similarity to the EC958 IS nor their flanking regions.

3. This is a multi-fasta file that contains the IS1 from EC958 with 100 bp flanking region. To determine whether the assembled contigs match these flanking regions, you firstly need to know the length of these sequences. To do this you can use the Count_IS_nt.py script in the Scripts folder.

4. This is the IS1 reference used to map all the reads to. 

5. This is a list of the read insert size and standard deviation, which is used in the velvet assembly. It looks like this::

	B36EC,248.06,93.39
	HVM1147,239.34,89.03

With the strain name at the beginning, the read insert size in the middle, and the standard deviation at the end. 
