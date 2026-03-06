#Load required packages----
library(tidyverse)
library(rentrez)
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostrings)
#BiocManager::install("muscle")
library(muscle)
#BiocManager::install("DECIPHER")
library(DECIPHER)
library(ape)
library(utils)
library(tidyverse)
library(phytools)
library(ggplot2)
source("scripts/functions.R")

#Data Acquisition and Exploration----

#Define taxon and gene of interest
taxon <- "Phaethornis"
gene <- "ND2"

raw_fasta_path <- sprintf("data_raw/%s/%s/gene_fetch/", taxon, gene)
processed_path  <- sprintf("data_processed/%s", taxon)

if (!dir.exists(raw_fasta_path)) {
  dir.create(raw_fasta_path, recursive = TRUE)
}

if (!dir.exists(processed_path)) {
  dir.create(processed_path, recursive = TRUE)
}

search_term <- sprintf("%s[ORGN] AND %s[Gene]", taxon, gene)

#Determine number of hits from searching nuccore database for specified taxon and gene of interest
gene_search_explore <- entrez_search(db = "nuccore", term = search_term)

#Count number of hits
max_hits <- gene_search_explore$count

#Use web history to fetch sequences
gene_search <- entrez_search(db = "nuccore", term = search_term, retmax = max_hits, use_history = T)

#Download Entrez_Functions.R in current directory and load to extract FASTA files from web_history objects
gene_fetch <- FetchFastaFiles(searchTerm = search_term, seqsPerFile = 100, fastaFileName = paste0(raw_fasta_path, "gene_fetch"))

#Merge seqeunce files into one dataframe
gene_seqs <- MergeFastaFiles(filePath = raw_fasta_path, filePattern = "gene_fetch*")

#Determine range of sequences to narrow search and reduce variability of sequence length
summary(nchar(gene_seqs$Sequence))
bp_hist(gene_seqs, "Sequence", binwidth = 150)

#Filter the sequence dataframe for sequences of similar lengths.
min_length <- 1000
max_length <- 1050
gene_seqs <- gene_seqs %>%
  filter(nchar(Sequence) >= min_length, nchar(Sequence) <= max_length)

#Check if filtering by sequence length worked. The result should contain sequences of less variability in length.
bp_hist(gene_seqs, "Sequence")
summary(nchar(gene_seqs$Sequence))

#Determine how many unknown nucleotides and gaps are present in data
table(str_count(gene_seqs$Sequence, "N"))
table(str_count(gene_seqs$Sequence, "-"))

#View sequences on a text edit application to look out for improper data
write.table(gene_seqs, file = sprintf("%s/%s_view_sequences.txt", processed_path, gene), sep = "\t", col.names = TRUE, row.names = FALSE)

#Filter sequences to remove sequences with N's at the beginning and end, and  sequences containing 0.1% N's or greater. Put filtered sequences into a seperate column
missing_data <- 0.01
gene_seqs<- gene_seqs %>%
  mutate(Sequence_Filtered = str_remove_all(Sequence, "^N+|N+$")) %>%
  filter(str_count(Sequence_Filtered, "N") <= (missing_data * str_count(Sequence)))

#Check if sequences were modified by trimming N's
all(are_equal <- gene_seqs$Sequence == gene_seqs$Sequence_Filtered)

#No modification other than filtering out sequences. Removing this column since it contains the same information.
gene_seqs$Sequence_Filtered <- NULL

#Extract species name and accession numbers into separate columns
gene_seqs$Species_Name <- word(gene_seqs$Title, 2L, 3L)
gene_seqs$AccessionNumber <- word(gene_seqs$Title, 1L)

#Rearrange the columns
gene_seqs <- gene_seqs[, c("Title","AccessionNumber", "Species_Name", "Sequence")]

#Count number of unique species in this dataframe
length(unique(gene_seqs$Species_Name))

#Randomly select 1 sequence to represent each species in this genus for downstream analysis
set.seed(1234)
gene_seqs_subset <- gene_seqs %>% 
  group_by(Species_Name) %>% 
  sample_n(1)

#Looking at the sequence lengths of our subset
bp_hist(gene_seqs_subset, "Sequence")

#Sequence Alignment----

#Prepare data for alignment by converting data types
gene_seqs_subset <- as.data.frame(gene_seqs_subset)
gene_seqs_subset$Sequence <- DNAStringSet(gene_seqs_subset$Sequence)
names(gene_seqs_subset$Sequence) <- gene_seqs_subset$Species_Name

# save sequences
write_csv(gene_seqs_subset, sprintf("data_processed/%s/%s_seq.csv", taxon, gene))

#Align sequences with default muscle settings
gene_seqs_subset.alignment <- DNAStringSet(muscle::muscle(gene_seqs_subset$Sequence), use.names = TRUE)

#View alignment
BrowseSeqs(gene_seqs_subset.alignment)

#Export alignment into a FASTA file. 
writeXStringSet(gene_seqs_subset.alignment, file = paste0("data_processed/Phaethornis/",gene, "_alignment", format = ".fasta"))
#File was viewed in MEGA. Since this gene codes for a protein, the sequences were translated with the genetic code set to vertebrate mitochondrial. Gaps were minimal, sequences were similar, and stop codons were not found in the reading frames. This indicates accuracy in the alignment, no contamination of data, and lack of reverse compliments.

#Distance Matrix, Clustering, and Mapping to Phylogenetic Tree----

#Convert alignment to a DNAbin data class in order to create a distance matrix to build a tree
dist_matrix <- dist.dna(as.DNAbin(gene_seqs_subset.alignment), model = "TN93")

#Create a neighbour joining tree with the distance matrix and plot it to new a dendogram
nj_tree <- nj(dist_matrix)
plot(nj_tree)

# save tree

