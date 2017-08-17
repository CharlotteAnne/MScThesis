# MSc Thesis
A collection of scripts for my Bioinformatics MSc thesis
  * peptide_database_MUSCLE.R - a script to take multiple STAR SJ.out.tab junction files, extend the coordinates and get the fasta sequences. These DNA files can then be translated using EMBOSS Transeq and cleaned up with cleanup.py
  * cleanup.py - for changing fasta headers on translated database files, trimming to the first stop codon, and trimming based on user defined minimum length (ie. 7aa for exons and the extension length for junctions)
