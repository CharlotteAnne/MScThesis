# MSc Thesis
A collection of scripts for my Bioinformatics MSc thesis
  * peptide_database_MUSCLE.R - a script to take multiple STAR SJ.out.tab junction files, extend the coordinates and get the fasta sequences. These DNA files can then be translated using EMBOSS Transeq and cleaned up with fasta_parsing/trim_fasta.py
  * fasta_parsing/trim_fasta.py - for changing fasta headers on translated database files, trimming to the first stop codon, and trimming based on user defined minimum length (ie. 7aa for exons and the extension length for junctions). The other short python programs in this folder are for counting reads in a fasta file and changing header names
