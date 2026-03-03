######################## LIBRARIES

library(Biostrings) #readDNAStringSet, pairwiseAlignment, mismatchTable, writeXStringSet
library(pwalign) #mismatchTable
library(dplyr) # data wrangling (mutate, select, bind_rows, group_by, summarise, count)
library(tidyr) # separate

# Bioconductor packages -- install once with:
# install.packages("BiocManager")
# BiocManager::install("Biostrings", update = FALSE, ask = FALSE)
# BiocManager::install("pwalign",    update = FALSE, ask = FALSE)





######################## READ DNA FILES AND FORMAT SEQ LABELS
# save to new directory to preserve originals
output_dir = "fastas-dna-clean"
dir.create(output_dir, showWarnings = FALSE)

# read all fasta files at once
fasta_files = list.files("fastas-dna", pattern = "\\.fa$", full.names = TRUE)

# format fasta labels to gene_population_chr_start_end
# neandertal camelCase must come before tolower() to preserve distinction between 3 types
lapply(fasta_files, function(f) {
    dna = readDNAStringSet(f)
    labels = names(dna)

    # Easy fixes
    labels = gsub("[:\\-]", "_", labels)
    labels = tolower(labels)
    labels = gsub("hum\\.fa", "ref", labels)

    labels = gsub("neandertal_altai", "neandertalAltai", labels)
    labels = gsub("neandertal_chagyrskaya", "neandertalChagyrskaya", labels)
    labels = gsub("neandertal_vindija", "neandertalVindija", labels)

    # Extract gene name and coords from ref entry
    ref_label = labels[grepl("_ref_", labels)]
    gene = gsub("_ref_.*", "", ref_label)
    coords = gsub(".*_ref_", "", ref_label)

    # Prefix labels missing the gene name
    labels = ifelse(!startsWith(labels, gene), paste0(gene, "_", labels), labels)

    # Append coords to short labels still ending in a string
    labels = ifelse(grepl("[a-z]$", labels), paste0(labels, "_", coords), labels)

    names(dna) = labels
    writeXStringSet(dna, file.path(output_dir, basename(f)))
})


# combine all cleaned fastas into single dataframe and save
dna_clean = bind_rows(lapply(list.files(output_dir, pattern = "\\.fa$", full.names = TRUE), function(f) {
    dna = readDNAStringSet(f)
    data.frame(label = names(dna), sequence = as.character(dna))
}))

write.csv(dna_clean, "data-dna clean.csv")






######################## VARIANT CALLING
#look for variants
mismatches <- function(query, ref) {
    result <- pwalign::pairwiseAlignment(ref, query, substitutionMatrix = "BLOSUM100") |>
        pwalign::mismatchTable()
    if (nrow(result) == 0) {
        return(data.frame(ID=names(query), Reference_DNA=NA, Sample_DNA=NA, Pos=NA))
    }
    result |>
        mutate(ID=names(query),
               Pos=PatternStart,
               Reference_DNA=as.character(PatternSubstring),
               Sample_DNA=as.character(SubjectSubstring)) |>
        select(ID, Reference_DNA, Sample_DNA, Pos)
}

fasta_files_clean = list.files("fastas-dna-clean", pattern = "\\.fa$", full.names = TRUE)
dna_list = setNames(lapply(fasta_files_clean, readDNAStringSet),
                     gsub("-all with humref.fa", "", basename(fasta_files_clean)))

df = bind_rows(lapply(dna_list, function(dna) {
  bind_rows(lapply(seq_along(dna[-1]), function(i) mismatches(dna[i+1], dna[1])))
}))

#check that all labels have the expected five columns
table(lengths(strsplit(df$ID, "_")))

#split into five columns
df = df |> separate(ID, into = c("gene", "population", "chr", "start", "end"), sep = "_")

write.csv(df, "data-dna variants.csv")





######################## PERCENT VARIATION
# how many variants are present on each gene
# calculate gene length from start/end, inclusive count
percent_variation = df |>
    mutate(BasePairs = as.numeric(end) - as.numeric(start) + 1) |>
    dplyr::count(gene, population, chr, start, BasePairs, sort = FALSE)
percent_variation$pervar = ((percent_variation$n / percent_variation$BasePairs) * 100)

write.csv(percent_variation, "data-percent variation.csv")

saveRDS(percent_variation, "pervar_dna_heatmaps.rds")







######################## SAVE PERCENT VARIATION DATA TABLES
# how many genes contain variants for each population (max=30, all genes)
percent_variation |> count(population, sort = FALSE)

# mean percentage of variation (based on n genes containing variants)
percent_variation_means_by_pop_genes = percent_variation |>
    group_by(population) |>
    summarise_at(vars(pervar), list(name = mean))
#saveRDS(dnaPerVar.tablemean, "data-genePerVarMean.rds")





######################## TIDY
rm(list=ls())
gc()
