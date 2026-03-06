# Credit: Jacqueline May (FetchFastaFiles, MergeFastaFiles)
# Modified by Neha Patel

# Fetch FASTA files from NCBI nuccore based on a provided search term.
fetch_fasta_files <- function(search_term, seqs_per_file = 100, fasta_file) {
  
  # search_term = character vector containing Entrez search term
  # seqs_per_file = number of sequences to write to each FASTA file
  # fasta_file = character vector containing name you want to give to the FASTA files you are fetching
  
  # Initial search for finding maximum number of hits
  search1 <- entrez_search(db = "nuccore", term = search_term)
  # Second search for obtaining max number of hits and their IDs
  search2 <- entrez_search(db = "nuccore", term = search_term, retmax = search1$count, use_history = T)
  
  # Fetch the sequences in FASTA format using the web_history object.
  for (start_rec in seq(0, search2$retmax, seqsPerFile)) {
    fname <- paste(fasta_file, start_rec, ".fasta", sep = "")
    recs <- entrez_fetch(db = "nuccore", web_history = search2$web_history, rettype = "fasta", retstart = start_rec, retmax = seqs_per_file)
    write(recs, fname)
    print(paste("Wrote records to ", fname, sep = ""))
  }
  
  return(search2)
  
}

# This function merges multiple FASTA files into one dataframe.
merge_fasta_files <- function(file_path, file_pattern) { 
  
  # filePath = Character vector containing path to FASTA files
  # filePattern = Character vector containing common pattern in FASTA file names
  
  # Read the FASTA files in.
  fasta_files <- list.files(path = file_path, pattern = file_pattern, full.names = TRUE)
  l_fastaFiles <- lapply(fasta_files, readDNAStringSet)
  
  # Convert them into dataframes.
  l_df_fasta_files <- lapply(l_fasta_files, function(x) data.frame(Title = names(x), Sequence = paste(x) ))
  
  # Combine the list of dataframes into one dataframe.
  df_seqs <- do.call("rbind", l_df_fasta_files)
  
  return(df_seqs)
  
}


bp_hist <- function(df, seq_col_name, binwidth = 5) {
  ggplot(df, aes(x = nchar(.data[[seq_col_name]]))) +
    geom_histogram(binwidth = binwidth, fill = "sky blue", color = "black") +
    labs(x = "Sequence Length (bp)", y = "Frequency") +
    ggtitle(paste("Distribution of Sequence Length")) +
    theme_bw()
}
