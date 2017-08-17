# Spliced Peptide Database from STAR junction file, using information about implied novel exons (from MAJIQ output)#
# Charlotte Capitanchik #
# Summer 2017 #


#library('stringi')
library('data.table')
library('dplyr')
library('stringr')
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))


# First read in all STAR junction files for the specific tissue. They are read into a single dataframe.
#--------------------------------------------------------------------------------------------------------

directory <- "/homes/s1207699/data/rna/rat_muscle/STAR_junction_files/"
file_names <- paste0(directory,dir(directory))
file_names <- file_names[grep('*SJ.out.tab',file_names)]
raw_star <- do.call(rbind,lapply(file_names,fread,header=FALSE, sep="\t",data.table=FALSE,stringsAsFactors = FALSE))

#--------------------------------------------------------------------------------------------------------
# The format for star junction files is tsv, with each column as follows:
# 1 = chromosome, 2 = junction start (first base of intron, 1 based), 
# 3 = junction end (last base of intron, 1 based), 4 = strand (0: undefined, 1: +, 2: -)
# 5 = intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT
# 6 = 0: unannotated, 1: annotated (but this can't be used for our analyis because in the second pass 
# novel junctions are considered annotated)
# 7 = number of uniquely mapping reads crossing the junction
# 8 = number of multi-mapping reads crossing the junction
# 9 = maximum spliced alignment overhang
# We want to first filter for > x unique read mapping , y multimaps & greater than z nt intron length.
# Then junctions will be filtered, then duplicates will be removed - duplicates aren't removed first because a 
# junction might be well supported in one sample but not another, & this info is lost if we condense the 
# list first.
#--------------------------------------------------------------------------------------------------------

num_reads <- 6
num_multi <- 0 
intron_length <- 60

log_file <- paste0("Minimum unique reads: ",num_reads,
                   "\nMaximum multimapping reads: ", num_multi,
                   "\nMinimum intron length: ",intron_length,"\n")

# filter junctions based on above criteria

temp_filt_star = raw_star[(raw_star$V7>num_reads)&(raw_star$V8<=num_multi)&((raw_star$V3+1)-raw_star$V2>intron_length),]

log_file <- paste0(log_file,"Number of junctions in raw concatenated dataframe: ", nrow(raw_star),"\n")
log_file <- paste0(log_file, "Number of junctions in filtered dataframe: ", nrow(temp_filt_star),"\n")

# Now we just take out the columns we need: chr, start, stop, strand & remove duplicates

junctions <- data.frame(temp_filt_star$V1, temp_filt_star$V2, temp_filt_star$V3, temp_filt_star$V4, stringsAsFactors = FALSE)
junctions <- unique(junctions[order(junctions[[2]]),])
colnames(junctions) <- c("chr","start","stop","strand")

# Change the strand info to +, - or .

junctions$strand <- gsub('0', '.',junctions$strand)
junctions$strand <- gsub('1', '+',junctions$strand)
junctions$strand <- gsub('2', '-',junctions$strand)


# Find out which junctions are novel and add a column annotating them as such

rat_introns <- fread("/homes/s1207699/data/rna/rat_genome_prim/sjdbList.out.tab",stringsAsFactors = FALSE)
rat_introns_srt <- rat_introns[order(rat_introns[[2]]),]
colnames(rat_introns_srt) <- c("chr","start","stop","strand")

# Get the annotated juncs

rat_introns_srt <- data.frame(rat_introns_srt)
annotated_juncs <- merge(junctions,rat_introns_srt)
annotated_juncs$novel <- rep(FALSE,nrow(annotated_juncs))

# Get the novel juncs

novel_juncs <- anti_join(junctions,annotated_juncs)
novel_juncs$novel <- rep(TRUE,nrow(novel_juncs))

# Check that the numbers add up

nrow(novel_juncs)+nrow(annotated_juncs) == nrow(junctions)

# Assemble the new dataframe
# At this stage can just take novel, or take both frames.
#junctions <- rbind(novel_juncs, annotated_juncs)
junctions <- novel_juncs
junctions <- junctions[order(junctions[[2]]),]


# We also want to get gene annotations for these junctions
# This rat_genes file was downloaded from biomart and then sort -u for duplicates
# Needed to sort out the strand information to be in the same format 
# awk  -F, 'BEGIN {FS=OFS="\t"} {gsub("-1","-",$4); gsub("1", "+",$4); print}' all_uniq_rat_genes.tsv > all_uniq_rat_genes_str.tsv

rat_genes_file <- "/homes/s1207699/data/rna/rat_annotation/all_uniq_rat_genes_str.tsv"
rat_genes <- fread(rat_genes_file, header=FALSE,sep="\t",data.table=FALSE)
rat_genes <- rat_genes[order(rat_genes[[2]]),]

# write the junctions to file so I can make a system call to bedtools
junctionsTEMP <- paste0(directory,"junctionsTEMP")
junctionsNEW <- paste0(directory,"junctionsNEW")
fwrite(junctions, file=junctionsTEMP,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)


# Explanation of the bedtools options:
# -wa -wb, want the matching rows of both files to be concatenated together
# I could have put -f 1.0, means the whole junction has to be contained within the annotated gene, 
# so chimeric junctions won't be assigned to multiple genes (duplicate rows), but this is at the cost 
# of not assigning a gene name potentially to alternative start and stop sites outwith the bounds of
# known genes, which I think is more important. So with this command there will be duplicated
# rows 
# -loj means that junctions that don't overlap with a gene will still printed  

system(paste0("bedtools intersect -wa -wb -loj -a ",junctionsTEMP," -b ",rat_genes_file," > ",junctionsNEW))
junctions_temp <- fread(junctionsNEW, header=FALSE, sep="\t",data.table=FALSE)


# Following the above comment we can add a column to let us know if a given junction has been assigned
# to more than one gene

junctions_temp$duplicated <- duplicated(junctions_temp[,1:4]) | duplicated(junctions_temp[,1:4], fromLast=TRUE)


# Downstream when we make the bed12 file having duplicates of junction coordinates might be a 
# problem -bed12 spec says that intervals shouldn't be repeated, also don't want duplicate sequences
# in the splice peptide database > so just keep the first time the junc appears with that annotation
# it will still be marked duplicate so we'll know later it stretches across two genes

junctions_temp <- junctions_temp %>% group_by(V1, V2, V3, V4) %>% filter(row_number() == 1)


junctions <- data.frame(junctions_temp$V1,junctions_temp$V2,junctions_temp$V3,junctions_temp$V4,junctions_temp$V5,junctions_temp$V9,junctions_temp$V10, junctions_temp$V11,junctions_temp$duplicated)
colnames(junctions) <- c("chr","start","stop","junc strand","novel","gene strand","gene name","gene id","duplicated")


log_file <- paste0(log_file,"Number of non-redundant filtered junctions: ", nrow(junctions),"\n")


#--------------------------------------------------------------------------------------------------------
# Make new dataframes for 'upstream' translations and 'downstream' translations
# I want to find out where the downstream end and upstream end (After extending by extension are still
# in annotated exons)
# These are made by extending the junctions x nucleotides, taking into consideration the bed format
# where start of interval is 0-based and end of interval is 1-based
#--------------------------------------------------------------------------------------------------------

extension <- 66

upstream <- data.frame(junctions$chr,junctions$start-extension-1,junctions$start-extension,junctions$'junc strand')
upstream$junctions.chr <- as.character(upstream$junctions.chr)
downstream <- data.frame(junctions$chr,junctions$stop-1+extension,junctions$stop+extension,junctions$'junc strand')
downstream$junctions.chr <- as.character(downstream$junctions.chr)



# To get a bed file of exons with gene IDs from a gtf file I used:
# grep exon Rattus_norvegicus.Rnor_6.0.80.gtf | cut -f1,4,5,7,9 | 
# sed 's/;.*//' | sed 's/gene_id//' | sed 's/"//' | sed 's/"//' | sort -u > rat_exons.bed
# Read in this exon bed file

exon_file <- "/homes/s1207699/data/rna/rat_annotation/rat_unique_exons.update.bed"
exons <- fread(exon_file, header=FALSE, sep="\t",data.table=FALSE)

# Sort the exons by start location

exons <- exons[order(exons$V2),]

----------------------------------------------------------------------------------------------------
# Find out whether downstream/upstream extensions end within annotated exons
----------------------------------------------------------------------------------------------------

# Temp bed files path

upstreamTEMP <- paste0(directory,"upstreamTEMP")
downstreamTEMP <- paste0(directory,"downstreamTEMP")

# Temp overlap files path

upstreamOVERLAP <- paste0(directory,"upstreamOVERLAP")
downstreamOVERLAP <- paste0(directory,"downstreamOVERLAP")


# Write the temp file and use bedtools intersect to get the exons that overlap

fwrite(upstream, file=upstreamTEMP,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)
system(paste0("bedtools intersect -wa -wb -a ",upstreamTEMP," -b ",exon_file," > ",upstreamOVERLAP))

# Read in the overlap temp file
upstreamOVERLAP_temp <- fread(upstreamOVERLAP, header=FALSE, sep="\t",data.table=FALSE)

# Rearrange the dataframe to get rid of redundant columns
upstreamOVERLAP_df <- data.frame(upstreamOVERLAP_temp$V1, 
                                 upstreamOVERLAP_temp$V2,
                                 upstreamOVERLAP_temp$V3,
                                 upstreamOVERLAP_temp$V4,
                                 upstreamOVERLAP_temp$V6,
                                 upstreamOVERLAP_temp$V7,
                                 upstreamOVERLAP_temp$V9)
colnames(upstreamOVERLAP_df) <- c("chr","start ext","end ext","strand","exon start","exon stop","gene id")


# Repeat the steps for the downstream coordinates
fwrite(downstream, file=downstreamTEMP,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)
system(paste0("bedtools intersect -wa -wb -a ",downstreamTEMP," -b ",exon_file," > ",downstreamOVERLAP))
downstreamOVERLAP_temp <- fread(downstreamOVERLAP, header=FALSE, sep="\t",data.table=FALSE)
downstreamOVERLAP_df <- data.frame(downstreamOVERLAP_temp$V1, 
                                   downstreamOVERLAP_temp$V2,
                                   downstreamOVERLAP_temp$V3,
                                   downstreamOVERLAP_temp$V4,
                                   downstreamOVERLAP_temp$V6,
                                   downstreamOVERLAP_temp$V7,
                                   downstreamOVERLAP_temp$V9)
colnames(downstreamOVERLAP_df) <- c("chr","start ext","end ext","strand","exon start","exon stop","gene id")

# Now in our original junction dataframe 'junctions' we can seperate into four seperate
# catagories:
# 1. Both downstream and upstream coordinates overlap with annotated exons
# 2. Upstream coordinate overlaps but downstream does not
# 3. Downstream coordinate overlaps but upstream does not
# 4. Neither coordinate overlaps
# We can seperate out these possibilities

# 1. Both downstream and upstream coordinates overlap with annotated exons

both_overlap_temp <- junctions[which((junctions$stop+extension) %in% downstreamOVERLAP_df$`end ext`),]
both_overlap <- both_overlap_temp[which((both_overlap_temp$start-extension) %in% upstreamOVERLAP_df$`end ext`),]

# 2. Upstream coordinate overlaps but downstream does not

up_overlap_temp <- junctions[which((junctions$start-extension) %in% upstreamOVERLAP_df$`end ext`),]
up_overlap <- up_overlap_temp[which((up_overlap_temp$stop+extension) %not in% downstreamOVERLAP_df$`end ext`),]

# 3. Downstream coordinate overlaps but upstream does not

down_overlap_temp <- junctions[which((junctions$stop+extension) %in% downstreamOVERLAP_df$`end ext`),]
down_overlap <- down_overlap_temp[which((down_overlap_temp$start-extension) %not in% upstreamOVERLAP_df$`end ext`),]

# 4. Neither coordinate overlaps

none_overlap_up_temp <- junctions[which((junctions$start-extension) %not in% upstreamOVERLAP_df$`end ext`),]
none_overlap_up <- none_overlap_up_temp[which((none_overlap_up_temp$stop+extension) %not in% downstreamOVERLAP_df$`end ext`),]

none_overlap_down_temp <- junctions[which((junctions$stop+extension) %not in% downstreamOVERLAP_df$`end ext`),]
none_overlap_down <-none_overlap_down_temp[which((none_overlap_down_temp$start-extension) %not in% upstreamOVERLAP_df$`end ext`),]
  
none_overlap <- rbind(none_overlap_up,none_overlap_down)  
none_overlap <- unique(none_overlap)
                       
# If this is correct then nrow(both+up+down+none) should equal nrow(junctions)

nrow(none_overlap)+nrow(down_overlap)+nrow(up_overlap)+nrow(both_overlap) == nrow(junctions)
                                                                                
# The percentage of junctions that would have intronic sequence translated if we blindly extend
# Note: This is likely an overestimate because a proportion of 'none_overlap' will be in novel exons,
# So the extension would be fine

(nrow(none_overlap)+nrow(down_overlap)+nrow(up_overlap))/nrow(junctions) * 100

# Percentage and number of each junction type that are novel

a <- nrow(both_overlap[which(both_overlap$novel==TRUE),])/nrow(both_overlap)
b <- nrow(both_overlap[which(both_overlap$novel==TRUE),])
c <- nrow(up_overlap[which(up_overlap$novel==TRUE),])/nrow(up_overlap)
d <- nrow(up_overlap[which(up_overlap$novel==TRUE),])
e <- nrow(down_overlap[which(down_overlap$novel==TRUE),])/nrow(down_overlap)
f <- nrow(down_overlap[which(down_overlap$novel==TRUE),])
g <- nrow(none_overlap[which(none_overlap$novel==TRUE),])/nrow(none_overlap)
h <- nrow(none_overlap[which(none_overlap$novel==TRUE),])

# Add to log file

log_file <- paste0(log_file,"% both-overlap that are novel junctions: ",a,"\n",
                   "Number of both-overlap novel junctions: ",b,"\n",
                   "% upstream-overlap that are novel junctions: ",c,"\n",
                   "Number of upstream-overlap novel junctions: ",d,"\n",
                   "% downstream-overlap novel junctions: ",e,"\n",
                   "Number of downstream-overlap novel junctions: ",f,"\n",
                   "% none-overlap novel junctions: ",g,"\n",
                   "Number of none-overlap novel junctions: ",h,"\n")

#------------------------------------------------------------------------------------------------
# Now that we have catagorised the junctions, we need to create complementary bed files that 
# can be used to extract the genomic sequence from the genome using bedtools 'get fasta' utility
#------------------------------------------------------------------------------------------------

# For now, I haven't written anything to process the different categories further but
# This could be future work, so in this case I'm just going to collapse all the categories

#------------------------------------------------------------------------------------------------
# To use the getfasta utility with 'split' to join the two chunks together I need to make a bed12
# file. This has the format:
# 1.chr 2.start 3.stop 4.name 5.score 6.strand 7.thickstart 8.thickstart 9.rgb 10.num blocks 
# 11.block sizes (comma seperated) 12.block starts, 0 first, relative to [[2]]
# The name column will be the fasta header

both_overlap<- rbind(both_overlap,up_overlap,down_overlap,none_overlap)
both_overlap <- both_overlap[order(both_overlap$start),]
bed12_both_overlap <- data.frame(both_overlap$chr,both_overlap$start-extension-1,both_overlap$stop+extension)

# I want the fasta header to look like:
# >ENSRNOG00000048039.19:1077-1299|19:1010-1076,19:1299-1365,novel:FALSE,gene:LOC679149,duplicated:FALSE 


bed12_both_overlap$name <- paste0(both_overlap$'gene id',".",both_overlap$chr,":",both_overlap$start,"-",both_overlap$stop,"|",both_overlap$chr,":",both_overlap$start-extension-1,"-",both_overlap$start-1,",",
       both_overlap$chr,":",both_overlap$stop,"-",both_overlap$stop+extension,",",
       "novel:",both_overlap$novel,",","gene:",both_overlap$'gene name',",",
       "duplicated:",both_overlap$duplicated)

# The next column is 'score' so we put all zeros

bed12_both_overlap$score <- rep(0,nrow(bed12_both_overlap))

# Next is strand, I'll use the junction strand - could use the gene strand here but there might
# be a case where the junction has been incorrectly annotated so safer to use the junction
# strand here

bed12_both_overlap$strand <- both_overlap$`junc strand`

# Next is three columns of zeros

bed12_both_overlap$thickstart <- rep(0,nrow(bed12_both_overlap))
bed12_both_overlap$thickstop <- rep(0,nrow(bed12_both_overlap))
bed12_both_overlap$itemRGB <- rep(0,nrow(bed12_both_overlap))

# Next is block count, will be 2 for all

bed12_both_overlap$blockcount <- rep(2,nrow(bed12_both_overlap))

# Next is block sizes, this was specified earlier as 'extension'

bed12_both_overlap$blocksizes <- rep(paste0(extension,",",extension,","),nrow(bed12_both_overlap))

# Next is block starts, with chr start as 0, so the first has to be 0
bed12_both_overlap$blockstarts <- paste0(0,",",(both_overlap$stop)-(both_overlap$start-extension-1))

# Seperate out the junctions where the strand is unknown (corresponding to non-canonical junction)
# These need to be translated in six frames
log_file <- paste0(log_file,"Number of non-canonical filtered junctions: ", nrow(bed12_both_overlap[bed12_both_overlap$strand==".",]),"\n")
bed12_both_overlap_SIX <- bed12_both_overlap[bed12_both_overlap$strand==".",]
bed12_both_overlap <- bed12_both_overlap[bed12_both_overlap$strand=="+" | bed12_both_overlap$strand=="-",]

# Write the bed file for three frames
bed12_both <- paste0(directory,"bed12_both_overlap")
fwrite(bed12_both_overlap, file=bed12_both,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)
# NOTE!!! -s means that getfasta will give you the reverse sequence if strand is minus!!!!!!!
# This is good because it saves faffing around later with the translation step.
system(paste0("bedtools sort -i ",bed12_both," > ",bed12_both,".sorted"))

# write bed file for six frames
bed12_both_SIX <- paste0(directory,"bed12_both_overlap_SIX")
fwrite(bed12_both_overlap_SIX, file=bed12_both_SIX,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)
# NOTE!!! -s means that getfasta will give you the reverse sequence if strand is minus!!!!!!!
# This is good because it saves faffing around later with the translation step.
system(paste0("bedtools sort -i ",bed12_both_SIX," > ",bed12_both_SIX,".sorted"))

#get fasta for three frames
system(paste0("bedtools getfasta -fi /homes/s1207699/data/rna/rat_genome_fasta/Rattus_norvegicus.Rnor_6.0.dna.primary.fa -bed ",
       bed12_both,".sorted -s -split -name -fo /homes/s1207699/data/rna/rat_muscle/STAR_junction_files/NOVEL_bed12_both_overlap.fasta"))

#get fasta for six frames
system(paste0("bedtools getfasta -fi /homes/s1207699/data/rna/rat_genome_fasta/Rattus_norvegicus.Rnor_6.0.dna.primary.fa -bed ",
              bed12_both_SIX,".sorted -split -name -fo /homes/s1207699/data/rna/rat_muscle/STAR_junction_files/NOVEL_bed12_both_overlap_SIX.fasta"))


# So bed12_both_overlap.fasta contains all the junction sequences

#--------------------------------------------------------------------------------------------------------------------
# Get the novel exon sequences from the MAJIQ PSI output & getFasta
#--------------------------------------------------------------------------------------------------------------------

# Read in the psi.tsv files, merge and split the exon coord column by ';'
muscle_psi_1 <- fread("/homes/s1207699/data/rna/MAJIQ_analysis/rat/GSE53960/muscle_liver_analysis/muscle_psi/muscle.psi_psi.tsv",stringsAsFactors = FALSE)
muscle_psi_2 <- fread("/homes/s1207699/data/rna/MAJIQ_analysis/rat/GSE41637/muscle_liver_analysis/muscle_psi/muscle.psi_psi.tsv",stringsAsFactors = FALSE)
muscle_psi_1 <- rbind(muscle_psi_1,muscle_psi_2)
muscle_exons <- data.frame(muscle_psi_1$chr, str_split_fixed(muscle_psi_1$`Exons coords`,";", Inf), muscle_psi_1$strand,muscle_psi_1$`Gene ID`,stringsAsFactors = FALSE)

# Find out how many exon columns we have now (number of columns minus chr, strand and id column)
exon_col_number <- ncol(muscle_exons)-3

# Now I want to seperate each exon column into a seperate 'bed' style dataframe, rbind and remove empty rows
Musc_exons <- data.frame(chr=character(),exon=character(),strand=character(),gene_id=character(), stringsAsFactors = FALSE)
for (i in 1:exon_col_number){
  hey <- 1+i
  musc_ex <- data.frame(muscle_exons[[1]],muscle_exons[[hey]],muscle_exons$muscle_psi_1.strand,muscle_exons$muscle_psi_1..Gene.ID.)
  colnames(musc_ex) <- c("chr","exon","strand","gene_id")
  Musc_exons <- rbind(musc_ex, Musc_exons)
}

# Get rid of empty & duplicated rows
Musc_exons <- unique(Musc_exons[!(Musc_exons$exon==""), ])

# Split the exon column into start and stop

Musc_exons <- data.frame(Musc_exons[[1]], str_split_fixed(Musc_exons[[2]],'-',2),Musc_exons[[3]],Musc_exons[[4]],stringsAsFactors = FALSE)

# Now subset for unannotated exons, earlier we read in a file of annotated exons as 'exons'
# Need to use bedtools overlap to find out which exons overlap with annotated exons
# As MAJIQ-defined exons will include alternative start and stop sites as entirely new exon boundaries
# As a first filter we can get rid of exons that are exactly the same...

annotated_Exons <- exons[,1:4]
colnames(annotated_Exons) <- c("chr","start","stop","strand")
colnames(Musc_exons) <- c("chr","start","stop","strand","gene_id")
Musc_exons$start <- as.numeric(Musc_exons$start)
Musc_exons$stop <- as.numeric(Musc_exons$stop)
novel_exons <- anti_join(Musc_exons, annotated_Exons)
novel_exons <- unique(novel_exons)

# To put the files in bed format have to convert to 0-based indexing
novel_exons$start <- novel_exons$start-1
annotated_Exons$start <- annotated_Exons$start-1

# Now use bedtools overlap to get rid of the overlapping ones
fwrite(novel_exons, file=paste0(directory,"novel_exon.bed"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)
fwrite(annotated_Exons, file=paste0(directory,"annotated_exon.bed"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)
system(paste0("bedtools intersect -v -a ",paste0(directory,"novel_exon.bed")," -b ",paste0(directory,"annotated_exon.bed")," > ",paste0(directory,"true_novel_exon.bed")))

novel_exons <- fread(paste0(directory,"true_novel_exon.bed"))
colnames(novel_exons) <- c("chr","start","stop","strand","gene_name")

# Seperate out exons with length 11
eleven_exons <- novel_exons[(novel_exons$stop-novel_exons$start==11),]
not_eleven_exons <- novel_exons[!(novel_exons$stop-novel_exons$start==11),]
nrow(not_eleven_exons)+nrow(eleven_exons) == nrow(novel_exons)

# Slight aside here: exons less than 21bp will be translated to aa sequence less than 7aa so will be filtered out
# (minimum match length), so to include them you'll need the flanking exon sequences too
# Below I take some of these NET sequences that are divisible by three and write them to a file so I can manually
# deal with them

less_than_seven <- novel_exons[(novel_exons$stop-novel_exons$start<21),]
greater <- novel_exons[(novel_exons$stop-novel_exons$start>132),]
NE_nets_names <- fread("/homes/s1207699/stuff_fromeric/RAT_NEMMnetlist_nets_3.bed")
NE_nets_names <- NE_nets_names$V4
less_than_seven_NETs <- less_than_seven[(less_than_seven$gene_name %in% NE_nets_names),]
less_than_seven_NETs <- less_than_seven_NETs[(less_than_seven_NETs$stop-less_than_seven_NETs$start)%%3==0,]
nrow(less_than_seven_NETs)

fwrite(less_than_seven_NETs, file=paste0(directory,"lessthansevennets.bed"),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)

# not_eleven_exons are good to go, but for the 11 exons need to expand both directions up to arbritary number (Say 100)
# get coverage and then take the interval with the biggest coverage as the interval for the sequence
# make a unique ID so we can match up the coverages 

not_eleven_exons$gene_id <- paste0(not_eleven_exons$chr,":", not_eleven_exons$start,"-",not_eleven_exons$stop)
eleven_exons$gene_id <- paste0(eleven_exons$chr,":", eleven_exons$start,"-",eleven_exons$stop)

# now for elevent exons make a downstream frame and an upstream frame
exon_extension <- 89

up_exon <- data.frame(eleven_exons$chr,as.numeric(eleven_exons$start-exon_extension), eleven_exons$stop, eleven_exons$strand,eleven_exons$gene_id)
down_exon <- data.frame(eleven_exons$chr,eleven_exons$start, as.numeric(eleven_exons$stop+exon_extension), eleven_exons$strand,eleven_exons$gene_id)
up_exon_file <- paste0(directory,"up_exon.bed")
down_exon_file <- paste0(directory,"down_exon.bed")
fwrite(up_exon,file=up_exon_file,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)
fwrite(down_exon,file=down_exon_file,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8) 

# need the bam files for this step

bam_files <- c("/homes/s1207699/data/rna/rat_muscle/GSE41637/rat_b_musc/rat_b_musc_prim_Aligned.sortedByCoord.out.filtered.bam",
               "/homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_M21_34_prim_Aligned.sortedByCoord.out.filtered.bam",
               "/homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_M21_12_prim_Aligned.sortedByCoord.out.filtered.bam",
               "/homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_F21_34_prim_Aligned.sortedByCoord.out.filtered.bam",
               "/homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_F21_12_prim_Aligned.sortedByCoord.out.filtered.bam"
               )
# this step can take some time 

system(paste0("multiBamSummary BED-file --BED ",up_exon_file," --bamfiles /homes/s1207699/data/rna/rat_muscle/GSE41637/rat_b_musc/rat_b_musc_prim_Aligned.sortedByCoord.out.filtered.bam /homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_M21_34_prim_Aligned.sortedByCoord.out.filtered.bam /homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_M21_12_prim_Aligned.sortedByCoord.out.filtered.bam /homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_F21_34_prim_Aligned.sortedByCoord.out.filtered.bam /homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_F21_12_prim_Aligned.sortedByCoord.out.filtered.bam -out ",directory,"up_exon_cov --outRawCounts ",directory,"up_exon_cov.tab"))
up_exon_cov <- fread(paste0(directory,"up_exon_cov.tab"))
# add up the coverage from every bam file
up_exon_cov <- data.frame(up_exon_cov[,1:3], rowSums(up_exon_cov[,4:ncol(up_exon_cov)]))

system(paste0("multiBamSummary BED-file --BED ",down_exon_file," --bamfiles /homes/s1207699/data/rna/rat_muscle/GSE41637/rat_b_musc/rat_b_musc_prim_Aligned.sortedByCoord.out.filtered.bam /homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_M21_34_prim_Aligned.sortedByCoord.out.filtered.bam /homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_M21_12_prim_Aligned.sortedByCoord.out.filtered.bam /homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_F21_34_prim_Aligned.sortedByCoord.out.filtered.bam /homes/s1207699/data/rna/rat_muscle/GSE53960/rat_musc_F21_12_prim_Aligned.sortedByCoord.out.filtered.bam -out ",directory,"down_exon_cov --outRawCounts ",directory,"down_exon_cov.tab"))
down_exon_cov <- fread(paste0(directory,"down_exon_cov.tab"))
# add up the coverage from every bam file
down_exon_cov <- data.frame(down_exon_cov[,1:3], rowSums(down_exon_cov[,4:ncol(down_exon_cov)]))

# give each file gene id columns
down_exon_cov$gene_id <- paste0(down_exon_cov[[1]],":",down_exon_cov[[2]],"-",down_exon_cov[[3]]-exon_extension)
up_exon_cov$gene_id <- paste0(up_exon_cov[[1]],":",up_exon_cov[[2]]+exon_extension,"-",up_exon_cov[[3]])

# merge
final_cov <- join_all(list(eleven_exons,up_exon_cov,down_exon_cov), by = 'gene_id')
colnames(final_cov) <- c("chr","start","stop","strand","gene_name","gene_id","up_chr","up_start","up_stop","up_cov","down_chr","down_start","down_stop","down_cov")
nrow(unique(final_cov))==nrow(unique(eleven_exons))

final_cov1 <- final_cov %>%
  filter(up_cov > down_cov) %>%
  select(chr,up_start,up_stop,strand,gene_name)
  
final_cov2 <- final_cov %>%
  filter(down_cov > up_cov) %>%
  select(chr,down_start,down_stop,strand,gene_name)

final_cov3 <- final_cov %>%
  filter(down_cov == up_cov) %>%
  select(chr,start,stop,strand,gene_name)

colnames(final_cov1) <- c("chr","start","stop","strand","gene_name")
colnames(final_cov2) <- c("chr","start","stop","strand","gene_name")
colnames(final_cov3) <- c("chr","start","stop","strand","gene_name")

nrow(final_cov1)+nrow(final_cov2)+nrow(final_cov3)==nrow(final_cov)

final_cov <- rbind(final_cov1,final_cov2,final_cov3)

# check we havent lost any exons along the way
nrow(unique(final_cov)) == nrow(unique(eleven_exons))

novel_exons_final <- rbind(final_cov,not_eleven_exons[,1:5])
nrow(novel_exons_final) == nrow(novel_exons)
novel_exon_filename <- paste0(directory,"novel_exons.bed")
# need to adjust novel_exons_final so that it's in bed12 format (Strand in the 6th column, name in the 4th)
novel_exons_final$name <- paste0(novel_exons_final$gene_name,".",novel_exons_final$chr,":",novel_exons_final$start,"-",novel_exons_final$stop,"|","novel_exon")
novel_exons_final <- data.frame(novel_exons_final[,1:3], novel_exons_final$name,rep(".",nrow(novel_exons_final)),novel_exons_final$strand)

# get rid of repeated interval, keep first gene annotation

novel_exons_final <- novel_exons_final %>% group_by(chr, start, stop) %>% filter(row_number() == 1)



fwrite(novel_exons_final,file=novel_exon_filename,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)

# Get fasta
system(paste0("bedtools getfasta -fi /homes/s1207699/data/rna/rat_genome_fasta/Rattus_norvegicus.Rnor_6.0.dna.primary.fa -bed ",
              novel_exon_filename," -s -name -fo /homes/s1207699/data/rna/rat_muscle/STAR_junction_files/novel_exons.fasta"))

# Print number of novel exons to log file
log_file <- paste0(log_file, "Number of novel exons: ",nrow(novel_exons_final),"\n")

#-------------------------------------------------------------------------------------------------------------------
# Get the retained intron sequences from the MAJIQ PSI output & getFasta
#-------------------------------------------------------------------------------------------------------------------

# Previously we have a dataframe containing the muscle psi output called 'muscle_psi_1'
# Subset for where IR Coords doesn't equal NA or nothing

IR <- muscle_psi_1[!is.na(muscle_psi_1$`IR coords`)&!(muscle_psi_1$`IR coords`==""),]

# Split coords into seperate columns 
IR_coords <- unique(data.frame(IR$chr,str_split_fixed(IR$`IR coords`,'-',2),IR$strand,IR$`Gene ID`, stringsAsFactors = FALSE))

# Filter intron size according to the intron_length variable set at the beginning to filter the
# SJ.out.file
IR_coords <- IR_coords[as.numeric(IR_coords$X2)-as.numeric(IR_coords$X1) > intron_length,]

# Make fasta header line
IR_name <- paste0(IR_coords$IR..Gene.ID.,".",IR_coords$IR.chr,":",IR_coords$X1,"-",IR_coords$X2,"|retained_intron")
IR_filename <- paste0(directory,"retained_introns.bed")

# Put into bed format
IR_file <- data.frame(IR_coords$IR.chr,as.numeric(IR_coords$X1)-1,IR_coords$X2,IR_name, rep('.',nrow(IR_coords)), IR_coords$IR.strand)
IR_file <- unique(IR_file)

# Write file and get fasta
fwrite(IR_file, file=IR_filename,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,nThread=8)
system(paste0("bedtools getfasta -fi /homes/s1207699/data/rna/rat_genome_fasta/Rattus_norvegicus.Rnor_6.0.dna.primary.fa -bed ",
              IR_filename," -s -name -fo /homes/s1207699/data/rna/rat_muscle/STAR_junction_files/retained_introns.fasta"))

# Print number of intron retentions to log file
log_file <- paste0(log_file, "Number of intron retentions: ",nrow(IR_file))

#-------------------------------------------------------------------------------------------------------------------
# Translate fasta files
#-------------------------------------------------------------------------------------------------------------------

# At this point we have generated temporary fasta files as follows:
# novel_exons.fasta
# retained_introns.fasta
# bed12_both_overlap.fasta
# Will translate each seperately as they need different processing rules
# using EMBOSS Transeq, with the 'trim' feature enabled
# transeq -sequence novel_exons.fasta -frame F -table 0 -sformat pearson -outseq novel_exons_trans.fasta
# Then - clean up! 



# Finish by writing log file
fileConn<-file(paste0(directory,"log_file.txt"))
writeLines(log_file, fileConn)
close(fileConn)

