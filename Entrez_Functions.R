# Entrez functions.

FetchFastaFiles <- function(searchTerm, seqsPerFile = 100, fastaFileName) {
  
  # This function will fetch FASTA files from NCBI nuccore based on a provided search term.
  
  # searchTerm = character vector containing Entrez search term
  # seqsPerFile = number of sequences to write to each FASTA file
  # fastaFileName = character vector containing name you want to give to the FASTA files you are fetching
  
  # Initial search for finding maximum number of hits
  search1 <- entrez_search(db = "nuccore", term = searchTerm)
  # Second search for obtaining max number of hits and their IDs
  search2 <- entrez_search(db = "nuccore", term = searchTerm, retmax = search1$count, use_history = T)
  
  # Fetch the sequences in FASTA format using the web_history object.
  for (start_rec in seq(0, search2$retmax, seqsPerFile)) {
    fname <- paste(fastaFileName, start_rec, ".fasta", sep = "")
    recs <- entrez_fetch(db = "nuccore", web_history = search2$web_history, rettype = "fasta", retstart = start_rec, retmax = seqsPerFile)
    write(recs, fname)
    print(paste("Wrote records to ", fname, sep = ""))
  }
  
  return(search2)
  
}

MergeFastaFiles <- function(filePattern) {
  
  # This function merges multiple FASTA files into one dataframe.
  
  # filePattern = Character vector containing common pattern in FASTA file names
  
  # Read the FASTA files in.
  fastaFiles <- list.files(pattern = filePattern)
  l_fastaFiles <- lapply(fastaFiles, readDNAStringSet)
  
  # Convert them into dataframes.
  l_dfFastaFiles <- lapply(l_fastaFiles, function(x) data.frame(Title = names(x), Sequence = paste(x) ))
  
  # Combine the list of dataframes into one dataframe.
  dfSeqs <- do.call("rbind", l_dfFastaFiles)
  
  return(dfSeqs)
  
}

