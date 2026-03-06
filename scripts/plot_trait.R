library(tidyverse)
library(phytools)
library(ggplot2)
library(plotrix) # needed for the dotTree function from phytools
source("scripts/functions.R")

# genus of interest
taxon <- "Phaethornis"
gene <- "ND2"

# read in sequences resulting from get_tree.R
seq_subset <- read_csv(sprintf("data_processed/%s/%s_seq.csv", taxon, gene))

#Proceeding to map traits onto the phylogenetic tree. Download trait data from AVONET dataset (titled AVONET Supplementary dataset 1.xlsx). To convert to a format that is more compatible with R, save worksheet of interest (AVONET_Raw_Data) as a CSV file with the name below. Import file into R using a relative path assuming the downloaded data is in your current working session. 
avonet_file_path <- "data_raw/Phaethornis/AVONET_Supplementary_dataset_1_import.csv"
AVONET <- read_csv(avonet_file_path, na = "NA")


#AVONET dataset has species names from multiple datasets. Birdlife was used because data was available for all species in this analysis.

#Extract tarsus length values from BirdLife from genus of interest into a dataframe
trait <- AVONET[grep(paste0("^",taxon," "), AVONET$Species1_BirdLife), c("Species1_BirdLife","Tarsus.Length")]

#Create a dataframe with species names and trait data of interest, filter out NA values, and determine mean value for each species in the dataset
mean_trait <- trait %>%
  filter(!is.na(Tarsus.Length)) %>%
  group_by(Species1_BirdLife) %>%
  summarize(Mean_Tarsus_Length = mean(Tarsus.Length))

#Merge dataframes to extract only trait data for species included in the previously created phylogenetic tree
seq_subset_mean_trait <- merge(seq_subset, mean_trait, by.x = "Species_Name", by.y = "Species1_BirdLife", all.x = TRUE)

#Confirming the sequence lengths has stayed consistent 
bp_hist(seq_subset_mean_trait, "Sequence")

#Create a named vector from the dataframe for compatibility with contMap from the phytools package
trait_vector <- setNames(seq_subset_mean_trait$Mean_Tarsus_Length, seq_subset_mean_trait$Species_Name)

#IUCN red list data was downloaded from https://www.iucnredlist.org/ by searching for keyword "Phaethornis". Import data into R using a relative path assuming the downloaded data is in your current working session.
iucn_file_path <- sprintf("data_raw/%s/IUCN_data.csv", taxon)
IUCN <- read_csv(iucn_file_path)

#Extract species names and corresponding red list category from dataframe
IUCN_subset <- IUCN[,c("scientificName","redlistCategory")]

#Merge red list category data for species present in the dataframe used to make the phylogenetic tree
seq_subset_mean_trait_IUCN <- merge(seq_subset_mean_trait, IUCN_subset, by.x = "Species_Name", by.y = "scientificName", all.x = TRUE)

#Create a named vector from the dataframe for compatibility with dotTree from the phytools package
IUCN_vector <- setNames(seq_subset_mean_trait_IUCN$redlistCategory, seq_subset_mean_trait_IUCN$Species_Name)

#Use the contMap function from the phytools package to map traits onto phylogenetic tree. Example provided in the contMap documentation was used to set features for this tree

trait.contMap<-contMap(nj_tree, trait_vector,plot=FALSE,res=200)
trait.contMap<-setMap(trait.contMap, c("white","#FFFFB2","#FECC5C","#FD8D3C","#E31A1C"))


#Prepare to add IUCN data by setting colours to each category displayed in IUCN dataframe
cols = setNames(c("navy", "royalblue","skyblue"), c("Near Threatened","Vulnerable","Least Concern"))

# plot from: http://www.phytools.org/Cordoba2017/ex/15/Plotting-methods.html and http://blog.phytools.org/2017/01/overlaying-contmap-style-continuous.html

if (!dir.exists("figures")) {
  dir.create("figures")
}

png(sprintf("figures/%s_%s_plot.png", taxon, gene), width=8, height=6, units="in", res = 300)

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

dev.off()

#Determine lambda parameter using phylosig from the phytools package. Code from: http://www.phytools.org/static.help/phylosig.html
lambda <- phylosig(nj_tree, trait_vector, method="lambda",test=TRUE,nsim=1000)
lambda
#Phylogenetic signal lambba is 0.35
