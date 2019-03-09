# -----------------------------------------------------------------------------#
# Sargassum Fungi Processing Raw Reads - Depends on "itsx_fastq_extractor.r"
# Author: Geoffrey Zahn
# Requirements: ITSx v 1.1b1
# -----------------------------------------------------------------------------#

# Prepare demultiplexed reads in Bash... ####
# Convert Fastq to Fasta for ITSx
# Run ITSx on all fasta files
# subset original Fastq files to only matching reads from ITSx output
# for ITS1fasta in *.ITS1.ITS1.fasta
# do grep "^>" $ITS1fasta > tempfile #get read names that ITSx wrote out
# sed -i 's/^>/@/' tempfile #change > to @
# grep -A3 -Fwf tempfile $(basename $ITS1fasta .fastq.fna.ITS1.ITS1.fasta).fastq > fastq_temp #Look up each read name in original fastq and write all 4 lines of each read to new file
# mv fastq_temp $(basename $ITS1fasta .fastq.fna.ITS1.ITS1.fasta).condensed.fastq # rename the file
# done


# Extract fastqs from ITSx output ####
fastas <- list.files(path = "./Fastqs", pattern = ".fastq.fna.ITS1.ITS1.fasta", full.names = TRUE)
rawfqs <- list.files(path = "./Fastqs", pattern = "condensed.fastq$", full.names = TRUE)
NewNames <- paste0(strsplit(fastas,".fastq.fna.ITS1.ITS1.fasta"),".ITS1.fastq")



#same order and lengths???
identical(length(rawfqs),length(fastas))
fastas[11]
rawfqs[11]
NewNames[11]

# No ITS1 found in REV reads

# Jack Darcy's script ... takes -f ITSx_fasta -q Original_fastq -o New_name 
# (original options changed to allow easy looping: changed defaults to "i")

for(i in 1:length(fastas)){
  source("./itsx_fastq_extractor.r") 
}


# Load packages ####
# library("devtools")
# devtools::install_github("benjjneb/dada2")
# devtools::install_github("joey711/phyloseq")
# devtools::install_github("bryandmartin/corncob")
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
library(DECIPHER)
library(decontam)
library(phangorn)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(gridExtra)
library(ggpubr)
library(dada2); packageVersion("dada2")
library(stringr)
library(purrr)
library(corncob)

# Load dada2 and prep ####


# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "./Fastqs" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- list.files(path)
fastqs <- fns[grepl("ITS1.fastq$", fns)] # CHANGE if different file extensions or to target only certain sequences
rm(fns)

fnFs <- sort(list.files(path, pattern="ITS1.fastq", full.names = TRUE))
# fnRs <- sort(list.files(path, pattern="_R2_001", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

# visualize a couple of fwd read quality profiles to help you decide reasonable filtration parameters
plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])

# Filter and trim ####
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
# filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, # fnRs, filtRs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=16) # On Windows set multithread=FALSE

head(out)
no.pass = which(as.data.frame(out)$reads.out == 0)
out[no.pass,]





# learn error rates ####
# Since some samples had zero reads pass QC, reassign filtFs and filtRs
filtFs <- sort(list.files(filtpath, full.names = TRUE, pattern = "_F_"))
# filtRs <- sort(list.files(filtpath, full.names = TRUE, pattern = "_R_"))

errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST = 20, nbases = 1e10)
# errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST = 20, nbases = 1e10)

# sanity check
plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)


# Dereplication ####
derepFs <- derepFastq(filtFs, verbose=TRUE)
# derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
# Since some samples were removed (no reads passed QC), reassign sample.names
sample.names <- sapply(strsplit(basename(filtFs), "-"), `[`, 1)
# sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
# names(derepRs) <- sample.names


# Sample inference ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
# dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


# Make a sequence table ####
seqtab <- makeSequenceTable(dadaFs)


# Remove Chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# Track Reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$filter.loss = (track[,1]-track[,2])/track[,1]
write.csv(track, file = "./Output/read_counts_at_each_step.csv", row.names = TRUE)




# import metadata ####
meta = read.csv("./metadata.csv")[,1:6]
row.names(meta) <- meta$IlluminaName
row.names(seqtab.nochim)

# reorder metadata
meta = meta[order(row.names(meta)),]

# Find controlsamples (extraction negatives) and clean ####
meta$controls <- meta$Island == "Blank"

# remove missing samples excluded due to poor QC 

row.names(seqtab.nochim) <- map(strsplit(row.names(seqtab.nochim), split = "_"),1)
good.samples = (row.names(meta) %in%  row.names(seqtab.nochim))
meta <- (meta[good.samples,])
rm(good.samples)


# find contaminants
# seqtab.nochim = readRDS(file = "./Output/clean_dada2_seqtable.RDS")
contams = isContaminant(seqtab.nochim, neg = meta$controls, normalize = TRUE)
table(contams$contaminant)  # No control samples passed QC...only 1 sequence made it through


# remove them
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$controls == FALSE,]
meta = meta[meta$controls == FALSE,]

# Remove all seqs with fewer than 100 nucleotides
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# Assign Taxonomy ####

# Save intermediate seqtab and re-load before assigning taxonomy to reduce virtual memory usage
saveRDS(seqtab.nochim, file = "./Output/clean_dada2_seqtable_NCBI-taxonomy.RDS")

seqtab.nochim <- readRDS(file = "./Output/clean_dada2_seqtable_NCBI-taxonomy.RDS")

taxa <- assignTaxonomy(seqtab.nochim, "./Taxonomy/sh_general_release_dynamic_s_01.12.2017_w_all-outgroups.fasta", multithread=20)

# Save intermediate files
write.csv(as.data.frame(seqtab.nochim), file = "./Output/SeqTable_no-chimera_no-contams.csv", row.names = TRUE, quote = FALSE)
saveRDS(seqtab.nochim, file = "./Output/clean_dada2_seqtable_NCBI-taxonomy.RDS")
saveRDS(taxa, file = "./Output/UNITE_Taxonomy_from_dada2.RDS")
seqs <- getSequences(seqtab.nochim)

# Hand off to Phyloseq ####

seqtab.nochim = readRDS(file = "./Output/clean_dada2_seqtable_NCBI-taxonomy.RDS")

taxa = readRDS(file = "./Output/UNITE_Taxonomy_from_dada2.RDS")

unique(taxa[,4])

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))


# Save RDS object for Phyloseq
saveRDS(ps, file = "./Output/clean_phyloseq_object.RDS")
ps = readRDS(file = "./Output/clean_phyloseq_object.RDS")
