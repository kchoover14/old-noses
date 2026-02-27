######################## LIBRARIES

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpmisc)
library(naniar)
library(Biostrings)
library(pwalign)
library(ape)
library(seqinr)
library(phylogram)
library(dendextend)
library(fitdistrplus)
library(car)

# Bioconductor packages -- install once with:
# install.packages("BiocManager")
# BiocManager::install("Biostrings", update = FALSE, ask = FALSE)
# BiocManager::install("pwalign",    update = FALSE, ask = FALSE)
# BiocManager::install("gdsfmt",     update = FALSE, ask = FALSE)
# BiocManager::install("SNPRelate",  update = FALSE, ask = FALSE)

options(width = 250)
options(max.print = .Machine$integer.max)


######################## VARIANT CALLING -- DNA SUBSTITUTIONS

or1a1.dna  = readDNAStringSet("fastas-DNA/OR1A1-all with humref.fa")
or1c1.dna  = readDNAStringSet("fastas-DNA/OR1C1-all with humref.fa")
or2a25.dna = readDNAStringSet("fastas-DNA/OR2A25-all with humref.fa")
or2b11.dna = readDNAStringSet("fastas-DNA/OR2B11-all with humref.fa")
or2c1.dna  = readDNAStringSet("fastas-DNA/OR2C1-all with humref.fa")
or2j2.dna  = readDNAStringSet("fastas-DNA/OR2J2-all with humref.fa")
or2j3.dna  = readDNAStringSet("fastas-DNA/OR2J3-all with humref.fa")
or2w1.dna  = readDNAStringSet("fastas-DNA/OR2W1-all with humref.fa")
or4e2.dna  = readDNAStringSet("fastas-DNA/OR4E2-all with humref.fa")
or4q3.dna  = readDNAStringSet("fastas-DNA/OR4Q3-all with humref.fa")
or5a1.dna  = readDNAStringSet("fastas-DNA/OR5A1-all with humref.fa")
or5an1.dna = readDNAStringSet("fastas-DNA/OR5AN1-all with humref.fa")
or5k1.dna  = readDNAStringSet("fastas-DNA/OR5K1-all with humref.fa")
or5p3.dna  = readDNAStringSet("fastas-DNA/OR5P3-all with humref.fa")
or6p1.dna  = readDNAStringSet("fastas-DNA/OR6P1-all with humref.fa")
or7c1.dna  = readDNAStringSet("fastas-DNA/OR7C1-all with humref.fa")
or7d4.dna  = readDNAStringSet("fastas-DNA/OR7D4-all with humref.fa")
or8b3.dna  = readDNAStringSet("fastas-DNA/OR8B3-all with humref.fa")
or8d1.dna  = readDNAStringSet("fastas-DNA/OR8D1-all with humref.fa")
or8k3.dna  = readDNAStringSet("fastas-DNA/OR8K3-all with humref.fa")
or10a6.dna = readDNAStringSet("fastas-DNA/OR10A6-all with humref.fa")
or10g3.dna = readDNAStringSet("fastas-DNA/OR10G3-all with humref.fa")
or10g4.dna = readDNAStringSet("fastas-DNA/OR10G4-all with humref.fa")
or10g7.dna = readDNAStringSet("fastas-DNA/OR10G7-all with humref.fa")
or10j5.dna = readDNAStringSet("fastas-DNA/OR10J5-all with humref.fa")
or11a1.dna = readDNAStringSet("fastas-DNA/OR11A1-all with humref.fa")
or51e1.dna = readDNAStringSet("fastas-DNA/OR51E1-all with humref.fa")
or51e2.dna = readDNAStringSet("fastas-DNA/OR51E2-all with humref.fa")
or51l1.dna = readDNAStringSet("fastas-DNA/OR51L1-all with humref.fa")
or56a4.dna = readDNAStringSet("fastas-DNA/OR56A4-all with humref.fa")

# Function adapted from:
# https://stackoverflow.com/questions/48402460/identifying-amino-acid-substitutions-from-local-alignments-in-r
mismatches.dna = function(query, ref) {
  pairwiseAlignment(ref, query, substitutionMatrix = "BLOSUM100") |>
    pwalign::mismatchTable() |>
    dplyr::mutate(
      ID            = names(query),
      Pos           = PatternStart,
      Reference_DNA = as.character(PatternSubstring),
      Sample_DNA    = as.character(SubjectSubstring)
    ) |>
    dplyr::select(ID, Reference_DNA, Sample_DNA, Pos)
}

sink("revised-substitutions-dna-raw data.txt", append = FALSE)
bind_rows(lapply(seq_along(or1a1.dna[-1]),  function(i) mismatches.dna(or1a1.dna[i+1],  or1a1.dna[1])))
bind_rows(lapply(seq_along(or1c1.dna[-1]),  function(i) mismatches.dna(or1c1.dna[i+1],  or1c1.dna[1])))
bind_rows(lapply(seq_along(or2a25.dna[-1]), function(i) mismatches.dna(or2a25.dna[i+1], or2a25.dna[1])))
bind_rows(lapply(seq_along(or2b11.dna[-1]), function(i) mismatches.dna(or2b11.dna[i+1], or2b11.dna[1])))
bind_rows(lapply(seq_along(or2c1.dna[-1]),  function(i) mismatches.dna(or2c1.dna[i+1],  or2c1.dna[1])))
bind_rows(lapply(seq_along(or2j2.dna[-1]),  function(i) mismatches.dna(or2j2.dna[i+1],  or2j2.dna[1])))
bind_rows(lapply(seq_along(or2j3.dna[-1]),  function(i) mismatches.dna(or2j3.dna[i+1],  or2j3.dna[1])))
bind_rows(lapply(seq_along(or2w1.dna[-1]),  function(i) mismatches.dna(or2w1.dna[i+1],  or2w1.dna[1])))
bind_rows(lapply(seq_along(or4e2.dna[-1]),  function(i) mismatches.dna(or4e2.dna[i+1],  or4e2.dna[1])))
bind_rows(lapply(seq_along(or4q3.dna[-1]),  function(i) mismatches.dna(or4q3.dna[i+1],  or4q3.dna[1])))
bind_rows(lapply(seq_along(or5a1.dna[-1]),  function(i) mismatches.dna(or5a1.dna[i+1],  or5a1.dna[1])))
bind_rows(lapply(seq_along(or5an1.dna[-1]), function(i) mismatches.dna(or5an1.dna[i+1], or5an1.dna[1])))
bind_rows(lapply(seq_along(or5k1.dna[-1]),  function(i) mismatches.dna(or5k1.dna[i+1],  or5k1.dna[1])))
bind_rows(lapply(seq_along(or5p3.dna[-1]),  function(i) mismatches.dna(or5p3.dna[i+1],  or5p3.dna[1])))
bind_rows(lapply(seq_along(or6p1.dna[-1]),  function(i) mismatches.dna(or6p1.dna[i+1],  or6p1.dna[1])))
bind_rows(lapply(seq_along(or7c1.dna[-1]),  function(i) mismatches.dna(or7c1.dna[i+1],  or7c1.dna[1])))
bind_rows(lapply(seq_along(or7d4.dna[-1]),  function(i) mismatches.dna(or7d4.dna[i+1],  or7d4.dna[1])))
bind_rows(lapply(seq_along(or8b3.dna[-1]),  function(i) mismatches.dna(or8b3.dna[i+1],  or8b3.dna[1])))
bind_rows(lapply(seq_along(or8d1.dna[-1]),  function(i) mismatches.dna(or8d1.dna[i+1],  or8d1.dna[1])))
bind_rows(lapply(seq_along(or8k3.dna[-1]),  function(i) mismatches.dna(or8k3.dna[i+1],  or8k3.dna[1])))
bind_rows(lapply(seq_along(or10a6.dna[-1]), function(i) mismatches.dna(or10a6.dna[i+1], or10a6.dna[1])))
bind_rows(lapply(seq_along(or10g3.dna[-1]), function(i) mismatches.dna(or10g3.dna[i+1], or10g3.dna[1])))
bind_rows(lapply(seq_along(or10g4.dna[-1]), function(i) mismatches.dna(or10g4.dna[i+1], or10g4.dna[1])))
bind_rows(lapply(seq_along(or10g7.dna[-1]), function(i) mismatches.dna(or10g7.dna[i+1], or10g7.dna[1])))
bind_rows(lapply(seq_along(or10j5.dna[-1]), function(i) mismatches.dna(or10j5.dna[i+1], or10j5.dna[1])))
bind_rows(lapply(seq_along(or11a1.dna[-1]), function(i) mismatches.dna(or11a1.dna[i+1], or11a1.dna[1])))
bind_rows(lapply(seq_along(or51e1.dna[-1]), function(i) mismatches.dna(or51e1.dna[i+1], or51e1.dna[1])))
bind_rows(lapply(seq_along(or51e2.dna[-1]), function(i) mismatches.dna(or51e2.dna[i+1], or51e2.dna[1])))
bind_rows(lapply(seq_along(or51l1.dna[-1]), function(i) mismatches.dna(or51l1.dna[i+1], or51l1.dna[1])))
bind_rows(lapply(seq_along(or56a4.dna[-1]), function(i) mismatches.dna(or56a4.dna[i+1], or56a4.dna[1])))
sink()


######################## VARIANT CALLING -- AMINO ACID SUBSTITUTIONS

or1a1.aa  = readAAStringSet("fastas-aa/or1a1-aa-all with humref.fa")
or1c1.aa  = readAAStringSet("fastas-aa/or1c1-aa-all with humref.fa")
or2a25.aa = readAAStringSet("fastas-aa/or2a25-aa-all with humref.fa")
or2b11.aa = readAAStringSet("fastas-aa/or2b11-aa-all with humref.fa")
or2c1.aa  = readAAStringSet("fastas-aa/or2c1-aa-all with humref.fa")
or2j2.aa  = readAAStringSet("fastas-aa/or2j2-aa-all with humref.fa")
or2j3.aa  = readAAStringSet("fastas-aa/or2j3-aa-all with humref.fa")
or2w1.aa  = readAAStringSet("fastas-aa/or2w1-aa-all with humref.fa")
or4e2.aa  = readAAStringSet("fastas-aa/or4e2-aa-all with humref.fa")
or4q3.aa  = readAAStringSet("fastas-aa/or4q3-aa-all with humref.fa")
or5a1.aa  = readAAStringSet("fastas-aa/or5a1-aa-all with humref.fa")
or5an1.aa = readAAStringSet("fastas-aa/or5an1-aa-all with humref.fa")
or5k1.aa  = readAAStringSet("fastas-aa/or5k1-aa-all with humref.fa")
or5p3.aa  = readAAStringSet("fastas-aa/or5p3-aa-all with humref.fa")
or6p1.aa  = readAAStringSet("fastas-aa/or6p1-aa-all with humref.fa")
or7c1.aa  = readAAStringSet("fastas-aa/or7c1-aa-all with humref.fa")
or7d4.aa  = readAAStringSet("fastas-aa/or7d4-aa-all with humref.fa")
or8b3.aa  = readAAStringSet("fastas-aa/or8b3-aa-all with humref.fa")
or8d1.aa  = readAAStringSet("fastas-aa/or8d1-aa-all with humref.fa")
or8k3.aa  = readAAStringSet("fastas-aa/or8k3-aa-all with humref.fa")
or10a6.aa = readAAStringSet("fastas-aa/or10a6-aa-all with humref.fa")
or10g3.aa = readAAStringSet("fastas-aa/or10g3-aa-all with humref.fa")
or10g4.aa = readAAStringSet("fastas-aa/or10g4-aa-all with humref.fa")
or10g7.aa = readAAStringSet("fastas-aa/or10g7-aa-all with humref.fa")
or10j5.aa = readAAStringSet("fastas-aa/or10j5-aa-all with humref.fa")
or11a1.aa = readAAStringSet("fastas-aa/or11a1-aa-all with humref.fa")
or51e1.aa = readAAStringSet("fastas-aa/or51e1-aa-all with humref.fa")
or51e2.aa = readAAStringSet("fastas-aa/or51e2-aa-all with humref.fa")
or51l1.aa = readAAStringSet("fastas-aa/or51l1-aa-all with humref.fa")
or56a4.aa = readAAStringSet("fastas-aa/or56a4-aa-all with humref.fa")

mismatches.aa = function(query, ref) {
  pairwiseAlignment(ref, query, substitutionMatrix = "BLOSUM100",
                    gapOpening = 3, gapExtension = 1) |>
    pwalign::mismatchTable() |>
    dplyr::mutate(
      ID           = names(query),
      Pos          = PatternStart,
      Reference_AA = as.character(PatternSubstring),
      Sample_AA    = as.character(SubjectSubstring)
    ) |>
    dplyr::select(ID, Reference_AA, Sample_AA, Pos)
}

sink("revised-substitutions-amino acid-raw data.txt", append = FALSE)
bind_rows(lapply(seq_along(or1a1.aa[-1]),  function(i) mismatches.aa(or1a1.aa[i+1],  or1a1.aa[1])))
bind_rows(lapply(seq_along(or1c1.aa[-1]),  function(i) mismatches.aa(or1c1.aa[i+1],  or1c1.aa[1])))
bind_rows(lapply(seq_along(or2a25.aa[-1]), function(i) mismatches.aa(or2a25.aa[i+1], or2a25.aa[1])))
bind_rows(lapply(seq_along(or2b11.aa[-1]), function(i) mismatches.aa(or2b11.aa[i+1], or2b11.aa[1])))
bind_rows(lapply(seq_along(or2c1.aa[-1]),  function(i) mismatches.aa(or2c1.aa[i+1],  or2c1.aa[1])))
bind_rows(lapply(seq_along(or2j2.aa[-1]),  function(i) mismatches.aa(or2j2.aa[i+1],  or2j2.aa[1])))
bind_rows(lapply(seq_along(or2j3.aa[-1]),  function(i) mismatches.aa(or2j3.aa[i+1],  or2j3.aa[1])))
bind_rows(lapply(seq_along(or2w1.aa[-1]),  function(i) mismatches.aa(or2w1.aa[i+1],  or2w1.aa[1])))
bind_rows(lapply(seq_along(or4e2.aa[-1]),  function(i) mismatches.aa(or4e2.aa[i+1],  or4e2.aa[1])))
bind_rows(lapply(seq_along(or4q3.aa[-1]),  function(i) mismatches.aa(or4q3.aa[i+1],  or4q3.aa[1])))
bind_rows(lapply(seq_along(or5a1.aa[-1]),  function(i) mismatches.aa(or5a1.aa[i+1],  or5a1.aa[1])))
bind_rows(lapply(seq_along(or5an1.aa[-1]), function(i) mismatches.aa(or5an1.aa[i+1], or5an1.aa[1])))
bind_rows(lapply(seq_along(or5k1.aa[-1]),  function(i) mismatches.aa(or5k1.aa[i+1],  or5k1.aa[1])))
bind_rows(lapply(seq_along(or5p3.aa[-1]),  function(i) mismatches.aa(or5p3.aa[i+1],  or5p3.aa[1])))
bind_rows(lapply(seq_along(or6p1.aa[-1]),  function(i) mismatches.aa(or6p1.aa[i+1],  or6p1.aa[1])))
bind_rows(lapply(seq_along(or7c1.aa[-1]),  function(i) mismatches.aa(or7c1.aa[i+1],  or7c1.aa[1])))
bind_rows(lapply(seq_along(or7d4.aa[-1]),  function(i) mismatches.aa(or7d4.aa[i+1],  or7d4.aa[1])))
bind_rows(lapply(seq_along(or8b3.aa[-1]),  function(i) mismatches.aa(or8b3.aa[i+1],  or8b3.aa[1])))
bind_rows(lapply(seq_along(or8d1.aa[-1]),  function(i) mismatches.aa(or8d1.aa[i+1],  or8d1.aa[1])))
bind_rows(lapply(seq_along(or8k3.aa[-1]),  function(i) mismatches.aa(or8k3.aa[i+1],  or8k3.aa[1])))
bind_rows(lapply(seq_along(or10a6.aa[-1]), function(i) mismatches.aa(or10a6.aa[i+1], or10a6.aa[1])))
bind_rows(lapply(seq_along(or10g3.aa[-1]), function(i) mismatches.aa(or10g3.aa[i+1], or10g3.aa[1])))
bind_rows(lapply(seq_along(or10g4.aa[-1]), function(i) mismatches.aa(or10g4.aa[i+1], or10g4.aa[1])))
bind_rows(lapply(seq_along(or10g7.aa[-1]), function(i) mismatches.aa(or10g7.aa[i+1], or10g7.aa[1])))
bind_rows(lapply(seq_along(or10j5.aa[-1]), function(i) mismatches.aa(or10j5.aa[i+1], or10j5.aa[1])))
bind_rows(lapply(seq_along(or11a1.aa[-1]), function(i) mismatches.aa(or11a1.aa[i+1], or11a1.aa[1])))
bind_rows(lapply(seq_along(or51e1.aa[-1]), function(i) mismatches.aa(or51e1.aa[i+1], or51e1.aa[1])))
bind_rows(lapply(seq_along(or51e2.aa[-1]), function(i) mismatches.aa(or51e2.aa[i+1], or51e2.aa[1])))
bind_rows(lapply(seq_along(or51l1.aa[-1]), function(i) mismatches.aa(or51l1.aa[i+1], or51l1.aa[1])))
bind_rows(lapply(seq_along(or56a4.aa[-1]), function(i) mismatches.aa(or56a4.aa[i+1], or56a4.aa[1])))
sink()


######################## DATA WRANGLING

# DNA
dnasubs.raw = read_table("revised-substitutions-dna-raw data.txt", skip = 1,
                         col_names = c("Number", "GenePop2", "Ref", "Alt", "Position"))
dnasubs.raw   = dplyr::select(dnasubs.raw, -Number)
dnasubs1      = replace_with_na(dnasubs.raw, replace = list(GenePop2 = "Reference_AA"))
dnasubs2      = replace_with_na(dnasubs1,    replace = list(Ref      = "Sample_AA"))
dnasubs3      = replace_with_na(dnasubs2,    replace = list(Alt      = "Pos"))
dnasubs       = na.omit(dnasubs3)
dnasubs$GenePop = dnasubs$GenePop2
dnasubs       = tidyr::separate(dnasubs, GenePop2, c("Gene", "Population"), "-")
write.csv(dnasubs, "revised-substitutions-dna-wrangled data.csv", quote = FALSE, row.names = FALSE)

# Amino acids
aasubs.raw = read_table("revised-substitutions-amino acid-raw data.txt", skip = 1,
                        col_names = c("Number", "GenePop2", "Ref", "Alt", "Position"))
aasubs.raw   = dplyr::select(aasubs.raw, -Number)
aasubs1      = replace_with_na(aasubs.raw, replace = list(GenePop2 = "Reference_AA"))
aasubs2      = replace_with_na(aasubs1,    replace = list(Ref      = "Sample_AA"))
aasubs3      = replace_with_na(aasubs2,    replace = list(Alt      = "Pos"))
aasubs       = na.omit(aasubs3)
aasubs$GenePop = aasubs$GenePop2
aasubs       = tidyr::separate(aasubs, GenePop2, c("Gene", "Population"), "-")
write.csv(aasubs, "revised-substitutions-amino acids-wrangled data.csv", quote = FALSE, row.names = FALSE)


######################## PERCENT VARIATION CALCULATION

# DNA -- count variants per gene per population, then compute percent variation
# dnasubs is kept in environment from wrangling above; no re-read from csv needed
dna.counts = dplyr::count(dnasubs, Gene, Population)

GeneVar.GeneCounts = dplyr::count(dna.counts, Population)
write.table(GeneVar.GeneCounts, "revised-data-nGeneVar.txt", quote = FALSE, row.names = FALSE)

dna.pervar = dplyr::mutate(dna.counts, BasePairs = case_when(
  Gene == 'or10a6' ~ 945,  Gene == 'or10g3' ~ 942,  Gene == 'or10g4' ~ 936,
  Gene == 'or10g7' ~ 936,  Gene == 'or10j5' ~ 930,  Gene == 'or11a1' ~ 948,
  Gene == 'or1a1'  ~ 930,  Gene == 'or1c1'  ~ 945,  Gene == 'or2a25' ~ 933,
  Gene == 'or2b11' ~ 954,  Gene == 'or2c1'  ~ 939,  Gene == 'or2j2'  ~ 939,
  Gene == 'or2j3'  ~ 936,  Gene == 'or2w1'  ~ 963,  Gene == 'or4q3'  ~ 942,
  Gene == 'or51l1' ~ 948,  Gene == 'or56a4' ~ 1098, Gene == 'or5a1'  ~ 948,
  Gene == 'or5an1' ~ 936,  Gene == 'or5k1'  ~ 927,  Gene == 'or5p3'  ~ 936,
  Gene == 'or6p1'  ~ 954,  Gene == 'or7c1'  ~ 963,  Gene == 'or7d4'  ~ 939,
  Gene == 'or8b3'  ~ 942,  Gene == 'or8d1'  ~ 927,  Gene == 'or8k3'  ~ 939,
  Gene == 'or4e2'  ~ 942,  Gene == 'or51e1' ~ 957,  Gene == 'or51e2' ~ 957
))
dna.pervar$PerVar = (dna.pervar$n / dna.pervar$BasePairs) * 100
write.table(dna.pervar, "revised-dna.pervar.csv", quote = FALSE, row.names = FALSE)

dnaPerVar.tablemean = dna.pervar |>
  dplyr::group_by(Population) |>
  dplyr::summarise_at(vars(PerVar), list(name = mean))
write.table(dnaPerVar.tablemean, "revised-data-GenePerVarMean.txt", quote = FALSE, row.names = FALSE)

dnaPerVar.table30 = dna.pervar |>
  dplyr::group_by(Population) |>
  dplyr::summarise_at(vars(PerVar), list(name = sum))
dnaPerVar.table30$dnaPerVar = dnaPerVar.table30$name / 30
write.table(dnaPerVar.table30, "revised-data-GenePerVarby30.txt", quote = FALSE, row.names = FALSE)

# Amino acids -- aasubs kept in environment from wrangling above
aa.counts = dplyr::count(aasubs, Gene, Population)

nProtVar.table = dplyr::count(aa.counts, Population)
write.table(nProtVar.table, "revised-data-nProtVar.txt", quote = FALSE, row.names = FALSE)

aa.pervar = dplyr::mutate(aa.counts, BasePairs = case_when(
  Gene == 'or10a6' ~ (945/3),  Gene == 'or10g3' ~ (942/3),  Gene == 'or10g4' ~ (936/3),
  Gene == 'or10g7' ~ (936/3),  Gene == 'or10j5' ~ (930/3),  Gene == 'or11a1' ~ (948/3),
  Gene == 'or1a1'  ~ (930/3),  Gene == 'or1c1'  ~ (945/3),  Gene == 'or2a25' ~ (933/3),
  Gene == 'or2b11' ~ (954/3),  Gene == 'or2c1'  ~ (939/3),  Gene == 'or2j2'  ~ (939/3),
  Gene == 'or2j3'  ~ (936/3),  Gene == 'or2w1'  ~ (963/3),  Gene == 'or4q3'  ~ (942/3),
  Gene == 'or51l1' ~ (948/3),  Gene == 'or56a4' ~ (1098/3), Gene == 'or5a1'  ~ (948/3),
  Gene == 'or5an1' ~ (936/3),  Gene == 'or5k1'  ~ (927/3),  Gene == 'or5p3'  ~ (936/3),
  Gene == 'or6p1'  ~ (954/3),  Gene == 'or7c1'  ~ (963/3),  Gene == 'or7d4'  ~ (939/3),
  Gene == 'or8b3'  ~ (942/3),  Gene == 'or8d1'  ~ (927/3),  Gene == 'or8k3'  ~ (939/3),
  Gene == 'or4e2'  ~ (942/3),  Gene == 'or51e1' ~ (957/3),  Gene == 'or51e2' ~ (957/3)
))
aa.pervar$PerVar = (aa.pervar$n / aa.pervar$BasePairs) * 100
write.table(aa.pervar, "revised-aa.pervar.csv", quote = FALSE, row.names = FALSE)

aaPerVar.tablemean = aa.pervar |>
  dplyr::group_by(Population) |>
  dplyr::summarise_at(vars(PerVar), list(name = mean))
write.table(aaPerVar.tablemean, "revised-data-ProtPerVarMean.txt", quote = FALSE, row.names = FALSE)

aaPerVar.table30 = aa.pervar |>
  dplyr::group_by(Population) |>
  dplyr::summarise_at(vars(PerVar), list(name = sum))
aaPerVar.table30$aaPerVar = aaPerVar.table30$name / 30
write.table(aaPerVar.table30, "revised-data-ProtPerVarby30.txt", quote = FALSE, row.names = FALSE)


######################## FIGURE 1 -- PERCENT VARIATION HEAT MAPS

pop.levels = c(
  'Ancient Neandertal Altai', 'Ancient Neandertal Chagyrskaya',
  'Ancient Neandertal Vindija', 'Ancient Denisovan', 'Ancient Ust Ishim',
  'Africa Esan', 'Africa Gambian', 'Africa Luhya', 'Africa Mende', 'Africa Yoruban',
  'Americas Black USA', 'Americas Black Caribbean',
  'Europe British', 'Europe Finn', 'Europe Spanish', 'Europe Tuscan',
  'Americas White USA', 'SW Asia Bengali', 'SW Asia Gujarati', 'SW Asia Punjabi',
  'SW Asia Tamil', 'SW Asia Telugu', 'Asia Dai', 'Asia Han Beijing',
  'Asia Han Southern', 'Asia Japanese', 'Asia Kinh',
  'Americas Columbian', 'Americas Mexican', 'Americas Puerto Rican', 'Americas Peruvian'
)

# dna.pervar is in environment from percent variation section above
dna.pervar$Gene       = as.factor(dna.pervar$Gene)
dna.pervar$Population = factor(dna.pervar$Population, levels = pop.levels)

pa = ggplot(dna.pervar, aes(x = Population, y = Gene, fill = PerVar)) +
  geom_tile(color = 'White', linewidth = 0.1) +
  scale_fill_viridis() +
  coord_equal() +
  xlab("") +
  ylab("Percent of Gene with DNA Variants") +
  theme_classic() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        text         = element_text(size = 14, family = "sans"))

# aa.pervar is in environment from percent variation section above
aa.pervar$Gene       = as.factor(aa.pervar$Gene)
aa.pervar$Population = factor(aa.pervar$Population, levels = pop.levels)

pb = ggplot(aa.pervar, aes(x = Population, y = Gene, fill = PerVar)) +
  geom_tile(color = 'White', linewidth = 0.1) +
  scale_fill_viridis() +
  coord_equal() +
  xlab("") +
  ylab("Percent of Amino Acids with Protein Variants") +
  theme_classic() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        text         = element_text(size = 14, family = "sans"))

png("revised-fig1.png", width = 45, height = 20, res = 400, units = "in")
plot_grid(pa, pb, nrow = 1, labels = c("A", "B"))
dev.off()


######################## FIGURE 2 -- FUNCTIONAL VARIATION (2016 VCF)

den2016 = read_excel("oldnoses-orfuncdata-2016-all.xlsx", sheet = "denisova-2016")
den2016$Odor         = as.factor(den2016$Odor)
den2016$OR           = as.factor(den2016$OR)
den2016$Human_log    = log(den2016$Human)
den2016$Denisova_log = log(den2016$Denisova)

nean2016 = read_excel("oldnoses-orfuncdata-2016-all.xlsx", sheet = "neandertal-2016")
nean2016$Odor           = as.factor(nean2016$Odor)
nean2016$OR             = as.factor(nean2016$OR)
nean2016$Human_log      = log(nean2016$Human)
nean2016$Neandertal_log = log(nean2016$Neandertal)

# Panel A: Neandertal regression
pa = ggplot(nean2016, aes(x = Human_log, y = Neandertal_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  xlim(-1, 2.5) + ylim(-1, 2.5) +
  labs(x = "Human Response(ln)", y = "Neandertal Response(ln)") +
  theme_classic(base_size = 20)

# Panel B: Denisovan regression
pb = ggplot(den2016, aes(x = Human_log, y = Denisova_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  xlim(-1, 2.5) + ylim(-1, 2.5) +
  labs(x = "Human Response(ln)", y = "Denisovan Response(ln)") +
  theme_classic(base_size = 20)

vcf2016long = read_excel("oldnoses-orfuncdata-2016-all.xlsx", sheet = "VCF2016-long")
vcf2016long$OR       = as.factor(vcf2016long$OR)
vcf2016long$Lineage3 = ordered(vcf2016long$Lineage, levels = c("Human", "Denisova", "Neandertal"))
vcf2016long$Lineage3 = gsub("Denisova", "Denisovan", vcf2016long$Lineage3)
vcf2016long$Response_log = log(vcf2016long$Response)

# Panel C: boxplot by lineage
pc = ggplot(vcf2016long, aes(x = Lineage3, y = Response_log, fill = Lineage3)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8, outlier.size = 2) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  scale_fill_brewer(palette = "PuOr") +
  labs(x = "", y = "Response(ln)") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

# Panel D: by OR
pd = ggplot(vcf2016long, aes(Response_log, Lineage3)) +
  geom_jitter(aes(col = OR)) +
  geom_smooth(aes(col = OR), method = "lm", se = FALSE) +
  scale_color_viridis_d() +
  labs(x = "Response by OR", y = "") +
  theme_classic(base_size = 20)

png("revised-fig2.png", width = 14, height = 10, res = 300, units = "in")
row1 = plot_grid(pa, pc, labels = c("a", "c"), nrow = 1)
row2 = plot_grid(pb, pd, labels = c("b", "d"), nrow = 1)
plot_grid(row1, row2, nrow = 2)
dev.off()


######################## PHYLOGENETIC ANALYSIS -- CONCATENATED MSA AND CLADOGRAM

or1a1  = readAAStringSet("fastas-aa/or1a1-aa-all with humref.fa")
or2a25 = readAAStringSet("fastas-aa/or2a25-aa-all with humref.fa")
or2c1  = readAAStringSet("fastas-aa/or2c1-aa-all with humref.fa")
or2j2  = readAAStringSet("fastas-aa/or2j2-aa-all with humref.fa")
or2j3  = readAAStringSet("fastas-aa/or2j3-aa-all with humref.fa")
or4e2  = readAAStringSet("fastas-aa/or4e2-aa-all with humref.fa")
or4q3  = readAAStringSet("fastas-aa/or4q3-aa-all with humref.fa")
or5a1  = readAAStringSet("fastas-aa/or5a1-aa-all with humref.fa")
or5an1 = readAAStringSet("fastas-aa/or5an1-aa-all with humref.fa")
or5k1  = readAAStringSet("fastas-aa/or5k1-aa-all with humref.fa")
or8k3  = readAAStringSet("fastas-aa/or8k3-aa-all with humref.fa")
or10g4 = readAAStringSet("fastas-aa/or10g4-aa-all with humref.fa")
or51e1 = readAAStringSet("fastas-aa/or51e1-aa-all with humref.fa")
or51l1 = readAAStringSet("fastas-aa/or51l1-aa-all with humref.fa")
or1c1  = readAAStringSet("fastas-aa/or1c1-aa-all with humref.fa")
or2b11 = readAAStringSet("fastas-aa/or2b11-aa-all with humref.fa")
or2w1  = readAAStringSet("fastas-aa/or2w1-aa-all with humref.fa")
or5p3  = readAAStringSet("fastas-aa/or5p3-aa-all with humref.fa")
or6p1  = readAAStringSet("fastas-aa/or6p1-aa-all with humref.fa")
or7c1  = readAAStringSet("fastas-aa/or7c1-aa-all with humref.fa")
or7d4  = readAAStringSet("fastas-aa/or7d4-aa-all with humref.fa")
or8b3  = readAAStringSet("fastas-aa/or8b3-aa-all with humref.fa")
or8d1  = readAAStringSet("fastas-aa/or8d1-aa-all with humref.fa")
or10a6 = readAAStringSet("fastas-aa/or10a6-aa-all with humref.fa")
or10g3 = readAAStringSet("fastas-aa/or10g3-aa-all with humref.fa")
or10g7 = readAAStringSet("fastas-aa/or10g7-aa-all with humref.fa")
or10j5 = readAAStringSet("fastas-aa/or10j5-aa-all with humref.fa")
or11a1 = readAAStringSet("fastas-aa/or11a1-aa-all with humref.fa")
or51e2 = readAAStringSet("fastas-aa/or51e2-aa-all with humref.fa")
or56a4 = readAAStringSet("fastas-aa/or56a4-aa-all with humref.fa")

uber.msa = AAStringSet(paste0(or1a1, or2a25, or2c1, or2j2, or2j3, or4e2,
                              or4q3, or5a1, or5an1, or5k1, or8k3, or10g4,
                              or51e1, or51l1, or1c1, or2b11, or2w1, or5p3,
                              or6p1, or7c1, or7d4, or8b3, or8d1, or10a6,
                              or10g3, or10g7, or10j5, or11a1, or51e2, or56a4))
names(uber.msa) = names(or1a1)
writeXStringSet(uber.msa, "revised-uber-msa.fa")

msf.ubertree          = read.alignment("revised-uber-msa.fa", format = "fasta")
distalign.ubertree    = dist.alignment(msf.ubertree, matrix = "identity")
distalign.nj.ubertree = nj(distalign.ubertree)

png("revised-figs2-cladogram.png", width = 10, height = 10, res = 300, units = "in")
plot(distalign.nj.ubertree, main = "", family = "sans", font = 3, ps = 12, cex = 1.5)
dev.off()


######################## SUPPLEMENTAL -- DISTRIBUTION CHECK

odors2016 = read_excel("oldnoses-orfuncdata-2016-all.xlsx", sheet = "VCF2016-long")
descdist(odors2016$Response, discrete = FALSE)
odors2016.lnorm = fitdist(odors2016$Response, "lnorm")
plot(odors2016.lnorm)
odors2016$Response_log = log(odors2016$Response)
descdist(odors2016$Response_log, discrete = FALSE)
odors2016.norm = fitdist(odors2016$Response_log, "norm")
plot(odors2016.norm)

odors2013 = read_excel("oldnoses-orfuncdata-2013-all.xlsx", sheet = "VCF2013-long")
descdist(odors2013$Response, discrete = FALSE)
odors2013.lnorm = fitdist(odors2013$Response, "lnorm")
plot(odors2013.lnorm)
odors2013$Response_log = log(odors2013$Response)
descdist(odors2013$Response_log, discrete = FALSE)
odors2013.norm = fitdist(odors2013$Response_log, "norm")
plot(odors2013.norm)


######################## SUPPLEMENTAL -- 2016 ACTIVE ORs REGRESSION

den2016.act = read_excel("oldnoses-orfuncdata-2016-active.xlsx", sheet = "denisova-2016")
den2016.act$Odor         = as.factor(den2016.act$Odor)
den2016.act$OR           = as.factor(den2016.act$OR)
den2016.act$Human_log    = log(den2016.act$Human)
den2016.act$Denisova_log = log(den2016.act$Denisova)

nean2016.act = read_excel("oldnoses-orfuncdata-2016-active.xlsx", sheet = "neandertal-2016")
nean2016.act$Odor           = as.factor(nean2016.act$Odor)
nean2016.act$OR             = as.factor(nean2016.act$OR)
nean2016.act$Human_log      = log(nean2016.act$Human)
nean2016.act$Neandertal_log = log(nean2016.act$Neandertal)

ps1 = ggplot(den2016.act, aes(x = Human_log, y = Denisova_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  xlim(-1, 2) + ylim(-1, 2) +
  labs(x = "Human Response(ln)", y = "Denisova Response(ln)") +
  theme_classic(base_size = 20)

ps2 = ggplot(nean2016.act, aes(x = Human_log, y = Neandertal_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  xlim(-1, 2) + ylim(-1, 2) +
  labs(x = "Human Response(ln)", y = "Neandertal Response(ln)") +
  theme_classic(base_size = 20)

png("revised-figs-2016-active-regression.png", width = 12, height = 6, res = 300, units = "in")
plot_grid(ps1, ps2, nrow = 1)
dev.off()


######################## SUPPLEMENTAL -- 2013 FUNCTIONAL VARIATION

alt2013 = read_excel("oldnoses-orfuncdata-2013-all.xlsx", sheet = "alt-2013")
alt2013$Odor      = as.factor(alt2013$Odor)
alt2013$OR        = as.factor(alt2013$OR)
alt2013$Human_log = log(alt2013$Human)
alt2013$Altai_log = log(alt2013$Altai)

den2013 = read_excel("oldnoses-orfuncdata-2013-all.xlsx", sheet = "den-2013")
den2013$Odor         = as.factor(den2013$Odor)
den2013$OR           = as.factor(den2013$OR)
den2013$Human_log    = log(den2013$Human)
den2013$Denisova_log = log(den2013$Denisova)

ust2013 = read_excel("oldnoses-orfuncdata-2013-all.xlsx", sheet = "ust-2013")
ust2013$Odor          = as.factor(ust2013$Odor)
ust2013$OR            = as.factor(ust2013$OR)
ust2013$Human_log     = log(ust2013$Human)
ust2013$Ust_Ishim_log = log(ust2013$Ust_Ishim)

ps3 = ggplot(ust2013, aes(x = Human_log, y = Ust_Ishim_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  xlim(-3, 3.5) + ylim(-3, 3.5) +
  labs(x = "Human Response(ln)", y = "Ust'-Ishim Response(ln)") +
  theme_classic(base_size = 20)

ps4 = ggplot(den2013, aes(x = Human_log, y = Denisova_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  xlim(-3, 3.5) + ylim(-3, 3.5) +
  labs(x = "Human Response(ln)", y = "Denisova Response(ln)") +
  theme_classic(base_size = 20)

ps5 = ggplot(alt2013, aes(x = Human_log, y = Altai_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  xlim(-3, 3.5) + ylim(-3, 3.5) +
  labs(x = "Human Response(ln)", y = "Altai Response(ln)") +
  theme_classic(base_size = 20)

vcf2013long = read_excel("oldnoses-orfuncdata-2013-all.xlsx", sheet = "VCF2013-long")
vcf2013long$OR      = as.factor(vcf2013long$OR)
vcf2013long$Lineage = ordered(vcf2013long$Lineage,
                               levels = c("Human", "Ust_Ishim", "Denisova", "Altai"))
vcf2013long$Lineage      = gsub("Denisova", "Denisovan", vcf2013long$Lineage)
vcf2013long$Response_log = log(vcf2013long$Response)

ps6 = ggplot(vcf2013long, aes(x = Lineage, y = Response_log, fill = Lineage)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8, outlier.size = 2) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  scale_fill_brewer(palette = "PuOr") +
  labs(x = "", y = "Response(ln)") +
  theme_classic(base_size = 20) +
  theme(legend.position  = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ps7 = ggplot(vcf2013long, aes(Response_log, Lineage)) +
  geom_jitter(aes(col = OR)) +
  geom_smooth(aes(col = OR), method = "lm", se = FALSE) +
  scale_color_viridis_d() +
  scale_y_discrete(labels = c("ALT", "DEN", "UST", "HUM")) +
  labs(x = "Response(ln) by OR", y = "") +
  theme_classic(base_size = 20)

png("revised-figs-2013-functional-variation.png", width = 20, height = 12, res = 300, units = "in")
row1 = plot_grid(ps3, ps4, ps5, labels = c("a", "b", "c"), nrow = 1)
row2 = plot_grid(ps6, ps7, labels = c("d", "e"), nrow = 1)
plot_grid(row1, row2, nrow = 2, rel_heights = c(1, 1))
dev.off()


######################## SUPPLEMENTAL -- 2013 ACTIVE ORs REGRESSION

den2013.act = read_excel("oldnoses-orfuncdata-2013-active.xlsx", sheet = "den-2013")
den2013.act$Odor         = as.factor(den2013.act$Odor)
den2013.act$OR           = as.factor(den2013.act$OR)
den2013.act$Human_log    = log(den2013.act$Human)
den2013.act$Denisova_log = log(den2013.act$Denisova)

alt2013.act = read_excel("oldnoses-orfuncdata-2013-active.xlsx", sheet = "alt-2013")
alt2013.act$Odor      = as.factor(alt2013.act$Odor)
alt2013.act$OR        = as.factor(alt2013.act$OR)
alt2013.act$Human_log = log(alt2013.act$Human)
alt2013.act$Altai_log = log(alt2013.act$Altai)

ust2013.act = read_excel("oldnoses-orfuncdata-2013-active.xlsx", sheet = "ust-2013")
ust2013.act$Odor          = as.factor(ust2013.act$Odor)
ust2013.act$OR            = as.factor(ust2013.act$OR)
ust2013.act$Human_log     = log(ust2013.act$Human)
ust2013.act$Ust_Ishim_log = log(ust2013.act$Ust_Ishim)

ps8 = ggplot(alt2013.act, aes(x = Human_log, y = Altai_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  labs(x = "Human Response(ln)", y = "Altai 2013 Response(ln)") +
  theme_classic(base_size = 20)

ps9 = ggplot(den2013.act, aes(x = Human_log, y = Denisova_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  labs(x = "Human Response(ln)", y = "Denisova 2013 Response(ln)") +
  theme_classic(base_size = 20)

ps10 = ggplot(ust2013.act, aes(x = Human_log, y = Ust_Ishim_log)) +
  geom_point(aes(colour = OR)) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               method = "lm", label.x = "left", label.y = "top", size = 6) +
  scale_color_viridis_d() +
  labs(x = "Human Response(ln)", y = "Ust'-Ishim 2013 Response(ln)") +
  theme_classic(base_size = 20)

png("revised-figs-2013-active-regression.png", width = 18, height = 6, res = 300, units = "in")
plot_grid(ps8, ps9, ps10, nrow = 1)
dev.off()
