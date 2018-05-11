#!/usr/bin/Rscript

suppressMessages(library(optparse))


default_expr_mat = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/sunhye_census_matrix_20170407.csv"
default_annotation = "~/single_cell_pipeline/scde_input/shl_0407_w_centroids_cell_info.csv"

# default_expr_mat = "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/transcripts.tpm_census_matrix.csv"
# default_annotation = "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/FACS_20171031_sunlee_sample_sheet.csv"

# default_clusters = "~/single_cell_pipeline/scde_input/diffex_by_trs_clusters_1_4/"
default_out = "./"

#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-e", "--expr_mat"), type="character", default=default_expr_mat,
              help="gene expression input filename [default= %default]", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=default_annotation,
              help="metadata about cells in input file [default= %default]", metavar="character"),
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
suppressMessages(library(gtools))
suppressMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

# make_option(c("-cl", "--clusters"), type="character", default=default_clusters,
#             help="file that defines comparison groups [default= %default]", metavar="character"),

genes_names <- c('733_234_transcripts_filtered.csv', '733_345_transcripts_filtered.csv', '737_transcripts_filtered.csv', 'ctrl_transcripts_filtered.csv')

new_genes_names <- paste0(opt$clusters, genes_names)

print("loading census matrix")
census_matrix = cataract::safe_read(opt$expr_mat)
print("loading cell metadata")
annotation = cataract::safe_read(opt$annotation)

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

plot_gene_by_treatment_and_day <- function(census_matrix, transcripts){
  # input = readline(prompt="Enter gene symbol to plot (ex. MYCN): ")
  
  plot_transcript <- function(transcript){
    filt_cm <- census_matrix[rownames(census_matrix) %in% transcript,]
    filt_cm <- tidyr::gather(filt_cm, "sample_id", "counts") %>% 
      inner_join(annotation)
    new_levels <- mixedsort(levels(filt_cm$day))
    filt_cm$day <- factor(filt_cm$day, levels = new_levels)
    bplot <- ggplot(data = filt_cm, aes(x=day, y=counts)) + 
      geom_boxplot() +
      geom_jitter(height = 0, width = 0.1) +
      # scale_x_discrete(mixedsort(levels(filt_cm$day))) +
      facet_grid(. ~ treatment_group) + 
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      labs(title = lookup_genes(transcript), subtitle = transcript)
    return(bplot)
  }
  
  # prep_transcripts <- lapply(transcripts, prep_transcript)
  bplots <- lapply(transcripts, plot_transcript)
  return(bplots)
}

question = "Choose from following:\n[P] Plot gene\n[L] Plot list of genes that are differentially expressed\n[X] Exit \n"

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
    transcripts <- transcripts[which(transcripts %in% rownames(census_matrix))]
    bplots <- plot_gene_by_treatment_and_day(census_matrix, transcripts) 
    pdf_out <- paste0(opt$out, input, ".pdf")
    pdf(pdf_out)
    invisible(lapply(bplots, print))
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
    transcripts <- transcripts[which(transcripts %in% rownames(census_matrix))]
    bplots <- plot_gene_by_treatment_and_day(census_matrix, transcripts) 
    pdf_out <- paste0(opt$out, diffex_group, ".pdf")
    pdf(pdf_out)
    invisible(lapply(bplots, print))
    dev.off()
    print(paste0("saving ", pdf_out))
  } else if (action == "I"){
    browser()
    print("a")
  } 
}

