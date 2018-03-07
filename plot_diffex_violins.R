#!/usr/bin/Rscript

suppressMessages(library(optparse))


default_expr_mat = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/sunhye_census_matrix_20170407.csv"
default_annotation = "~/single_cell_pipeline/scde_input/shl_0407_w_centroids_cell_info.csv"
default_clusters = "~/single_cell_pipeline/scde_input/diffex_by_trs_clusters_1_4/"
default_out = "./"

#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-e", "--expr_mat"), type="character", default=default_expr_mat,
              help="gene expression input filename [default= %default]", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=default_annotation,
              help="metadata about cells in input file [default= %default]", metavar="character"),
  make_option(c("-cl", "--clusters"), type="character", default=default_clusters,
              help="file that defines comparison groups [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NA,
              help="output file name [default= %default]", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  stop("Please provide all necessary arguments.", call.=FALSE)
}

suppressMessages(library(tidyverse))
suppressMessages(library(biomaRt))
suppressMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

genes_names <- c('733_234_transcripts_filtered.csv', '733_345_transcripts_filtered.csv', '737_transcripts_filtered.csv', 'ctrl_transcripts_filtered.csv')

new_genes_names <- paste0(opt$clusters, genes_names)

census_matrix = read.table(opt$expr_mat, sep = "\t", header = TRUE)
annotation = read.table(opt$annotation, sep = "\t", header = TRUE)


lookup_genes <- function(txname){
  
  txs <- transcripts(edb, filter = TxIdFilter(txname),
                     columns = c("symbol"))
  return(txs$symbol)
  
}

lookup_transcripts <- function(genename){
  
  txs <- transcripts(edb, filter = GenenameFilter(genename),
                     columns = c("tx_id"))
  return(txs$tx_id)
  
}

plot_gene_by_treatment_and_cluster <- function(census_matrix, transcripts){

  # input = readline(prompt="Enter gene symbol to plot (ex. MYCN): ")

  plot_transcript <- function(transcript){
    filter_census_matrix <- census_matrix[rownames(census_matrix) %in% transcript,]
    filter_census_matrix0 <- tidyr::gather(filter_census_matrix, "cell_id", "counts") %>% 
      inner_join(annotation) %>% 
      mutate(cluster_733_234 = rowSums(.[, c(8,9,10)], na.rm = TRUE)) %>% 
      mutate(cluster = as.character(cluster_733_234)) %>% 
      mutate(treatment_group = ifelse(treatment_group == "sh733", "sh733_234", as.character(treatment_group)))
    
    subset_733_345 <- filter_census_matrix0[filter_census_matrix0$treatment_group=="sh733_234",] %>%
      mutate(treatment_group = "sh733_345") %>% 
      mutate(cluster = centroid_733_345)
    
    filter_census_matrix0 <- rbind(filter_census_matrix0, subset_733_345)
    
    filter_census_matrix0 <- filter_census_matrix0[filter_census_matrix0$cluster != 0,] %>% 
      mutate(cluster = as.character(cluster)) %>% 
      dplyr::filter(!is.na(cluster))

    vplot <- ggplot(data = filter_census_matrix0, aes(x=cluster, y=counts)) + 
      geom_violin() +
      geom_jitter(height = 0, width = 0.1) +
      facet_grid(. ~ treatment_group) + 
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      labs(title = lookup_genes(transcript), subtitle = transcript)
    return(vplot)
  }
  
  # prep_transcripts <- lapply(transcripts, prep_transcript)
  vplots <- lapply(transcripts, plot_transcript)
  return(vplots)
}

question = "Choose from following:\n[P] Plot gene\n[L] Plot list of genes that are differentially expressed\n[X] Exit "

prompt_genes_names <- paste(basename(new_genes_names), collapse = " ")

while (TRUE) {
  # if (ctr < 3){
  #   print(question)
  # }
  
  cat(question)
  action <- toupper(readLines("stdin",n=1))
  # action <- readline(prompt=question)
  
  if(action == "X"){
    break
  } else if (action == "P"){
    
    # plot gene expression across all 16 sets ---------------------------------
    
    gene_question <- "Enter gene symbol to plot (ex. MYCN): "
    cat(gene_question)
    input <- readLines("stdin", n=1)
    transcripts <- lookup_transcripts(input)
    
    vplots <- plot_gene_by_treatment_and_cluster(census_matrix, transcripts) 
    pdf_out <- paste0(opt$out, input, ".pdf")
    pdf(pdf_out)
    invisible(lapply(vplots, print))
    dev.off()
    print(paste0("saving ", pdf_out))
  } else if (action == "L"){
    
    treatment_question <- paste("Enter treatment_group to analyze (", prompt_genes_names, "): ")
    
    cat(treatment_question)
    diffex_group <- readLines("stdin",n=1)
    # diffex_group <- readline(prompt=treatment_question)
    diffex_genes <- read.csv(new_genes_names[grep(diffex_group, new_genes_names)])
    
    top_n_question <- "how many genes do you want to plot? "
    cat(top_n_question)
    top_n <- readLines("stdin", n=1)
    
    transcripts <- rownames(diffex_genes[c(1:top_n),])
    vplots <- plot_gene_by_treatment_and_cluster(census_matrix, transcripts) 
    pdf_out <- paste0(opt$out, diffex_group, ".pdf")
    pdf(pdf_out)
    invisible(lapply(vplots, print))
    dev.off()
    print(paste0("saving ", pdf_out))
  } else if (action == "I"){
    browser()
    print(a)
  } 
}

