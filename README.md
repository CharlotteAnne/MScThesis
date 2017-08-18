# MSc Thesis
A collection of scripts for my Bioinformatics MSc thesis

Note! To view R markdown HTML files you can use RawGit https://rawgit.com/

  * peptide_database_MUSCLE.R - a script to take multiple STAR SJ.out.tab junction files, extend the coordinates and get the fasta sequences. These DNA files can then be translated using EMBOSS Transeq and cleaned up with fasta_parsing/trim_fasta.py
  * fasta_parsing/trim_fasta.py - for changing fasta headers on translated database files, trimming to the first stop codon, and trimming based on user defined minimum length (ie. 7aa for exons and the extension length for junctions). The other short python programs in this folder are for counting reads in a fasta file and changing header names
  * tissue_specific_MAJIQ_analysis.Rmd - code for parsing MAJIQ tsv output file of multiple tissue comparisons, reformatting to junctions and scoring tissue specificity, graphs using the tissue-scoring system are found here
  * SPLICING_SUMMARY.nb.html - a summary of all the MAJIQ analyses with a few graphs
  * NET_categories_upset.nb.html - using upset plots to categorise net genes
  * deseq_pca.nb.html - PCA plots for all samples and code for generating lists of expressed genes
  * MAJIQ_paper_NET_data_analysis.nb.html - heatmaps generated from the MAJIQ paper data, with proportion of nets spliced
  * all_tissue_splice_Analysis.nb.html - analysis for rat tissue comparisons, with heatmap
  * Supporting_Reads_Dist.nb.html - exploring the coverage of novel junctions to decide on filters for peptide database
