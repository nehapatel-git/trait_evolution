#Load required packages----
library(rentrez)
library(Biostrings)
library(tidyverse)
library(muscle)
library(DECIPHER)
library(ape)
library(utils)
library(tidyverse)
library(phytools)
library(ggplot2)
source("scripts/Entrez_Functions.R")

#Functions----

#Function to generate a histogram for # bp of sequences
plotbpLengthHist <- function(df, sequenceColumn, binwidth = 5) {
  ggplot(df, aes(x = nchar(.data[[sequenceColumn]]))) +
    geom_histogram(binwidth = binwidth, fill = "blue", color = "black") +
    labs(x = "Sequence Length (bp)", y = "Frequency") +
    ggtitle(paste("Distribution of Sequence Length"))
}

#Data Acquisition and Exploration----

#Define taxa and gene of interest
taxa <- "Phaethornis"
gene <- "ND2"

search_term <- sprintf("%s[ORGN] AND %s[Gene]", taxa, gene)

#Determine number of hits from searching nuccore database for specified taxa and gene of interest
gene_search_explore <- entrez_search(db = "nuccore", term = search_term)

#Count number of hits
maxHits <- gene_search_explore$count

#Use web history to fetch sequences
gene_search <- entrez_search(db = "nuccore", term = search_term, retmax = maxHits, use_history = T)

#Download Entrez_Functions.R in current directory and load to extract FASTA files from web_history objects
gene_fetch <- FetchFastaFiles(searchTerm = search_term, seqsPerFile = 100, fastaFileName = "raw_data/gene_fetch/gene_fetch.fasta")

#Merge seqeunce files into one dataframe
gene_seqs <- MergeFastaFiles(filePath = "raw_data/gene_fetch", filePattern = "gene_fetch*")

#Determine range of sequences to narrow search and reduce variability of sequence length
summary(nchar(gene_seqs$Sequence))
plotbpLengthHistogram(gene_seqs, "Sequence", binwidth = 150)

#Filter the sequence dataframe for sequences of similar lengths.
min_length <- 1000
max_length <- 1050
gene_seqs <- gene_seqs %>%
  filter(nchar(Sequence) >= min_length, nchar(Sequence) <= max_length)

#Check if filtering by sequence length worked. The result should contain sequences of less variability in length.
plotbpLengthHist(gene_seqs, "Sequence")
summary(nchar(gene_seqs$Sequence))

#Determine how many unknown nucleotides and gaps are present in data
table(str_count(gene_seqs$Sequence, "N"))
table(str_count(gene_seqs$Sequence, "-"))

#View sequences on a text edit application to look out for improper data
write.table(dfGene, file = "view_sequences.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

#Filter sequences to remove sequences with N's at the beginning and end, and  sequences containing 0.1% N's or greater. Put filtered sequences into a seperate column
missing.data <- 0.01
dfGene <- dfGene %>%
  mutate(Sequence_Filtered = str_remove_all(Sequence, "^N+|N+$")) %>%
  filter(str_count(Sequence_Filtered, "N") <= (missing.data * str_count(Sequence)))

#Check if sequences were modified by trimming N's
all(are_equal <- dfGene$Sequence == dfGene$Sequence_Filtered)
#No modification other than filtering out sequences. Removing this column since it contains the same information.
dfGene$Sequence_Filtered <- NULL

#Extract species name and accession numbers into separate columns
dfGene$Species_Name <- word(dfGene$Title, 2L, 3L)
dfGene$AccessionNumber <- word(dfGene$Title, 1L)

#Rearrange the columns
dfGene <- dfGene[, c("Title","AccessionNumber", "Species_Name", "Sequence")]

#Count number of unique species in this dataframe
length(unique(dfGene$Species_Name))

#Randomly select 1 sequence to represent each species in this genus for downstream analysis
set.seed(1234)
dfGene_Subset <- dfGene %>% 
  group_by(Species_Name) %>% 
  sample_n(1)

#Looking at the sequence lengths of our subset
plotSequenceLengthHistogram(dfGene_Subset, "Sequence")

#Sequence Alignment----

#Prepare data for alignment by converting data types
dfGene_Subset <- as.data.frame(dfGene_Subset)
dfGene_Subset$Sequence <- DNAStringSet(dfGene_Subset$Sequence)
names(dfGene_Subset$Sequence) <- dfGene_Subset$Species_Name

#Align sequences with default muscle settings
dfGene_Subset.alignment <- DNAStringSet(muscle::muscle(dfGene_Subset$Sequence), use.names = TRUE)

#View alignment
BrowseSeqs(dfGene_Subset.alignment)

#Export alignment into a FASTA file. 
writeXStringSet(dfGene_Subset.alignment, file = "alignment.fasta", format = "fasta")
#File was viewed in MEGA. Since this gene codes for a protein, the sequences were translated with the genetic code set to vertebrate mitochondrial. Gaps were minimal, sequences were similar, and stop codons were not found in the reading frames. This indicates accuracy in the alignment, no contamination of data, and lack of reverse compliments.

#Distance Matrix, Clustering, and Mapping to Phylogenetic Tree----

#Convert alignment to a DNAbin data class in order to create a distance matrix to build a tree
dist_matrix <- dist.dna(as.DNAbin(dfGene_Subset.alignment), model = "TN93")

#Create a neighbour joining tree with the distance matrix and plot it to new a dendogram
nj_tree <- nj(dist_matrix)
plot(nj_tree)


#Proceeding to map traits onto the phylogenetic tree. Download trait data from AVONET dataset (titled AVONET Supplementary dataset 1.xlsx). To convert to a format that is more compatible with R, save worksheet of interest (AVONET_Raw_Data) as a CSV file with the name below. Import file into R using a relative path assuming the downloaded data is in your current working session. 
avonet_file_path <- "AVONET Supplementary dataset 1_import.csv"
AVONET <- read_csv(avonet_file_path, na = "NA")

#Setting name of genus to genus of interest to extract appropriate data from database
genus <- "Phaethornis"

#AVONET dataset has species names from multiple datasets. The following code was also tested with species names from BirdTree and resulted in a missing value for one of the species. Analysis was proceeded with BirdLife species names instead.

#Extract tarsus length values from BirdLife from genus of interest into a dataframe
genus_trait <- AVONET[grep(paste0("^",genus," "), AVONET$Species1_BirdLife), c("Species1_BirdLife","Tarsus.Length")]

#Create a dataframe with species names and trait data of interest, filter out NA values, and determine mean value for each species in the dataset
genus_mean_trait <- genus_trait %>%
  filter(!is.na(Tarsus.Length)) %>%
  group_by(Species1_BirdLife) %>%
  summarize(Mean_Tarsus_Length = mean(Tarsus.Length))

#Merge dataframes to extract only trait data for species included in the previously created phylogenetic tree
dfGene_mean_trait <- merge(dfGene_Subset, genus_mean_trait, by.x = "Species_Name", by.y = "Species1_BirdLife", all.x = TRUE)

#Confirming the sequence lengths has stayed consistent 
plotSequenceLengthHistogram(dfGene_mean_trait, "Sequence")

#Create a named vector from the dataframe for compatibility with contMap from the phytools package
trait_vector <- setNames(dfGene_mean_trait$Mean_Tarsus_Length, dfGene_mean_trait$Species_Name)

#IUCN red list data was downloaded from https://www.iucnredlist.org/ by searching for keyword "Phaethornis". Import data into R using a relative path assuming the downloaded data is in your current working session.
iucn_file_path <- "IUCN_data.csv"
IUCN <- read_csv(iucn_file_path)

#Extract species names and corresponding red list category from dataframe
IUCN_subset <- IUCN[,c("scientificName","redlistCategory")]

#Merge red list category data for species present in the dataframe used to make the phylogenetic tree
dfGene_mean_trait_IUCN <- merge(dfGene_mean_trait, IUCN_subset, by.x = "Species_Name", by.y = "scientificName", all.x = TRUE)

#Create a named vector from the dataframe for compatibility with dotTree from the phytools package
IUCN_vector <- setNames(dfGene_mean_trait_IUCN$redlistCategory, dfGene_mean_trait_IUCN$Species_Name)

#Use the contMap function from the phytools package to map traits onto phylogenetic tree. Example provided in the contMap documentation was used to set features for this tree
?contMap
trait.contMap<-contMap(nj_tree, trait_vector,plot=FALSE,res=200)
trait.contMap<-setMap(trait.contMap, c("white","#FFFFB2","#FECC5C","#FD8D3C","#E31A1C"))


#Prepare to add IUCN data by setting colours to each category displayed in IUCN dataframe
cols = setNames(c("navy", "royalblue","skyblue"), c("Near Threatened","Vulnerable","Least Concern"))

#Code in lines 177-185 from: http://www.phytools.org/Cordoba2017/ex/15/Plotting-methods.html and http://blog.phytools.org/2017/01/overlaying-contmap-style-continuous.html

#Plot red list category as traits to the phylogenetic tree with the dotTree function from phytools
dotTree(nj_tree,IUCN_vector,legend = FALSE,length=5,fsize=0.8,lwd=7,ftype="i", colors = cols)

#Plot the contMap tree onto the dotTree
plot(trait.contMap$tree,colors=trait.contMap$cols,add=TRUE,ftype="off",lwd=5,
     xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
     ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)

#Add a legend for the contMap tree to represent tarsus lengths
add.color.bar(0.3*max(nodeHeights(nj_tree)),trait.contMap$cols,title="Tarsus Length (mm)",lims=trait.contMap$lims,digits=3,prompt=FALSE,x=0.08*max(nodeHeights(nj_tree)),y=0.4*(1+par()$usr[3]),lwd=4,fsize=0.8,subtitle="")

#Add a legend for the dotTree to represent red list categories
legend("bottomright", legend = c("Near Threatened", "Vulnerable", "Least Concern"), col = c("navy", "royalblue", "skyblue"), title = "Red List Category", pch = 16, pt.cex=1.8, cex=0.8, bty="n")


#Determine lambda parameter using phylosig from the phytools package. Code from: http://www.phytools.org/static.help/phylosig.html
lambda <- phylosig(nj_tree, trait_vector, method="lambda",test=TRUE,nsim=1000)
lambda
#Phylogenetic signal lambba is 0.35
