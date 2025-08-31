library("optparse")

#Read command-line options
option_list <- list(
  make_option(c("--annot"), type="character", default=NULL,
              help="Path to the annotations file made by prepareAnnotations.R.", metavar = "path"),
  make_option(c("--batch"), type="character", default = "NULL",
              help = "Specify analysis batch. Value '1 200' would split the gene ids into 200 batches and run the first batch.", metavar = "path"),
  make_option(c("--out"), type="character", default = "NULL",
              help = "Path to the output folder.", metavar = "path"),
  make_option(c("--fill"), type="logical", default = TRUE,
              help = "Fill alternative internal exons for promoter and 3'end events..", metavar = "path"),
  make_option(c("--cage"), type="character", default = NULL,
              help = "Path to the CAGE annotations file.", metavar = "path"),
  make_option(c("--start_end_diff"), type="integer", default = 20,
              help = "Minimal difference (in basepairs) between the alternative promoters or 3'ends.", metavar = "path")
)
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

annot_file = opt$annot
batch_string = opt$batch
out_dir = opt$out
fill_internal_exons = opt$fill
cage_file = opt$cage
start_end_diff = opt$start_end_diff
dir.create(out_dir)

#Import other dependencies
library("txrevise")
library("purrr")
suppressMessages(library("dplyr"))
suppressMessages(library("rtracklayer"))

constructEventMetadata <- function(transcript_ids){
  # Parse IDs without assuming no periods in gene_id or transcript_id
  # Pattern: (gene_id).(grp_id).(event_type).(transcript_id)
  # where grp_id starts with "grp_" and event_type has no periods
  
  regex <- "^(.*?)\\.(grp_\\d+)\\.(upstream|downstream|contained)\\.(.*)$"
  ensembl_gene_id <- sub(regex, "\\1", transcript_ids)
  grp_id <- sub(regex, "\\2", transcript_ids)
  event_type <- sub(regex, "\\3", transcript_ids)
  ensembl_transcript_id <- sub(regex, "\\4", transcript_ids)
  
  dplyr::tibble(
    transcript_id = as.character(transcript_ids),
    ensembl_gene_id = ensembl_gene_id,
    grp_id = grp_id,
    event_type = event_type,
    ensembl_transcript_id = ensembl_transcript_id
  ) |>
    dplyr::mutate(gene_id = paste(ensembl_gene_id, grp_id, event_type, sep = "."))
}

#Import prepared transcript annotations
txrevise_data = readRDS(annot_file)

#Import CAGE promoter annotations and merge them into Ensembl annotations
if(!is.null(cage_file)){
  cage_data = readRDS(cage_file)
  new_exons = c(txrevise_data$exons, cage_data$exons)
  new_metadata = dplyr::select(txrevise_data$transcript_metadata,
                               ensembl_gene_id, ensembl_transcript_id,
                               longest_start, longest_end, cds_start_NF,
                               cds_end_NF, cds_start_end_NF)
  new_metadata = dplyr::bind_rows(new_metadata, cage_data$transcript_metadata)
  txrevise_data = list(exons = new_exons, cdss = txrevise_data$cdss, transcript_metadata = new_metadata)
}

#### Split genes into batches ####
batch_vector = as.integer(unlist(strsplit(batch_string, split = " ")))
gene_ids = unique(txrevise_data$transcript_metadata$ensembl_gene_id)
batch_size = ceiling(length(gene_ids)/batch_vector[2])
batches = txrevise::splitIntoBatches(length(gene_ids), batch_size)
selection = batches == batch_vector[1]
selected_gene_ids = gene_ids[selection]

#Set up output file names
batch_id = paste(batch_vector, collapse = "_")
error_file = file.path(out_dir, paste0("failed_genes.batch_",batch_id, ".txt"))

#Only proceed with event construction if there are any genes in the list
if (length(selected_gene_ids) > 0){
  #Construct events
  gene_ids_list = setNames(as.list(selected_gene_ids), selected_gene_ids)

  #Construct alternative events and remove failed genes
  safe_construct = purrr::safely(constructAlternativeEventsWrapper)
  alt_events = purrr::map(gene_ids_list, ~safe_construct(., txrevise_data$transcript_metadata,
                                                         txrevise_data$exons,
                                                         txrevise_data$cdss,
                                                         max_internal_diff = 10,
                                                         max_start_end_diff = start_end_diff,
                                                         fill_internal = fill_internal_exons)$result)
  failed_genes = purrr::map_lgl(alt_events, is.null)
  alt_events = alt_events[!failed_genes] #Remove failed genes

  #Flatten
  alt_events = purrr::flatten(alt_events) %>% txrevise::flattenAlternativeEvents(min_alt_event_count = 1)

  #Construct event metadata
  event_metadata = constructEventMetadata(names(alt_events))

  #Iterate over groups and positions and export annotations to gff files
  for (group in c("grp_1", "grp_2")){
    for (event in c("upstream", "downstream", "contained")){
      selected_events = dplyr::filter(event_metadata, grp_id == group, event_type == event)
      out_file = file.path(out_dir, paste("txrevise",group, event, batch_id, "gff3", sep = "."))
      #Only write output if there are any events
      if(nrow(selected_events) > 0){
        gff_granges = txrevise::transcriptsToAnnotations(alt_events[selected_events$transcript_id], event_metadata)
        rtracklayer::export.gff3(gff_granges, out_file)
      } else { #Write empty output files
        file.create(out_file)
      }
    }
  }

  #Write failed genes
  failed_names = names(which(failed_genes))
  if(length(failed_names) > 0){
    write.table(failed_names, error_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
} else {
  #If there are now genes in the batch then write empty output files
  for (group in c("grp_1", "grp_2")){
    for (event in c("upstream", "downstream", "contained")){
      out_file = file.path(out_dir, paste("txrevise",group, event, batch_id, "gff3", sep = "."))
      file.create(out_file)
    }
  }
}
