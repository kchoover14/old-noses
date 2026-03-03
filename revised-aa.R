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
output_dir = "fastas-aa-clean"
dir.create(output_dir, showWarnings = FALSE)

# read all fasta files at once
fasta_files = list.files("fastas-aa", pattern = "\\.fa$", full.names = TRUE)

# format fasta labels to gene_population_chr_start_end
# neandertal camelCase must come before tolower() to preserve distinction between 3 types
lapply(fasta_files, function(f) {
    aa = readAAStringSet(f)
    labels = names(aa)

    # Easy fixes
    labels = gsub("[:\\-]", "_", labels)
    labels = tolower(labels)
    labels = gsub("humref\\.fa", "ref", labels)

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

    names(aa) = labels
    writeXStringSet(aa, file.path(output_dir, basename(f)))
})


# combine all cleaned fastas into single dataframe and save
aa_clean = bind_rows(lapply(list.files(output_dir, pattern = "\\.fa$", full.names = TRUE), function(f) {
    aa = readAAStringSet(f)
    data.frame(label = names(aa), sequence = as.character(aa))
}))
#saveRDS(aa_clean, "data-aa_clean.rds")






######################## VARIANT CALLING - AMINO ACIDS
mismatches_aa <- function(query, ref) {
    result <- pwalign::pairwiseAlignment(ref, query, substitutionMatrix = "BLOSUM62") |>
        pwalign::mismatchTable()
    if (nrow(result) == 0) {
        return(data.frame(ID=names(query), Reference_AA=NA, Sample_AA=NA, Pos=NA))
    }
    result |>
        mutate(ID=names(query),
               Pos=PatternStart,
               Reference_AA=as.character(PatternSubstring),
               Sample_AA=as.character(SubjectSubstring)) |>
        select(ID, Reference_AA, Sample_AA, Pos)
}

fasta_files_clean = list.files("fastas-aa-clean", pattern = "\\.fa$", full.names = TRUE)
aa_list = setNames(lapply(fasta_files_clean, readAAStringSet),
                    gsub("-all with humref.fa", "", basename(fasta_files_clean)))

df = bind_rows(lapply(aa_list, function(aa) {
    bind_rows(lapply(seq_along(aa[-1]), function(i) mismatches_aa(aa[i+1], aa[1])))
}))

# check that all labels have the expected five columns
table(lengths(strsplit(df$ID, "_")))

# split into five columns
df = df |> separate(ID, into = c("gene", "population", "chr", "start", "end"), sep = "_")

#saveRDS(df, "data-dna_variants.rds")





######################## PERCENT VARIATION
# how many variants are present on each gene
# calculate gene length from start/end, inclusive count
percent_variation = df |>
    mutate(BasePairs = as.numeric(end) - as.numeric(start) + 1) |>
    dplyr::count(gene, population, chr, start, BasePairs, sort = FALSE)
percent_variation$pervar = ((percent_variation$n / percent_variation$BasePairs) * 100)

saveRDS(percent_variation, "pervar_aa_heatmaps.rds")




######################## SAVE PERCENT VARIATION DATA TABLES

# how many genes contain variants for each population (max=30, all genes)
df |> count(population, sort = FALSE)

# mean percentage of variation (based on n genes containing variants)
#Shows that when this population varies, this is how much it varies on average
percent_variation_means_by_pop_proteins = percent_variation |>
    group_by(population) |>
    summarise_at(vars(pervar), list(name = mean))
#saveRDS(dnaPerVar.tablemean, "data-genePerVarMean.rds")

# mean percentage of variation (based on 30 genes total)
#Shows this population's overall variation rate across 30 genes
n_genes = 30 #total possible number of genes
percent_variation_means_by_all_proteins = percent_variation |>
    group_by(population) |>
    summarise(pervar = sum(pervar) / n_genes)

#saveRDS(dnaPerVar.table30, "data-genePerVarby30.rds")





######################## TIDY
rm(list=ls())
gc()
