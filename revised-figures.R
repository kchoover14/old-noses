######################## LIBRARIES
# load libraries
library(dplyr) #data wrangling
library(ggplot2) #plotting
library(viridis) #color-blind friendly palette
library(ggpubr) #f1 ggarrange
library(readxl) #read excel files of functional data
library(Biostrings) #read aa sequence
library(ape) #alignment
library(seqinr) #read aa strings
library(phylogram) #create tree

######################## HEAT MAPS
########### Population names
pop_names <- c(
    acb = "Americas Black Caribbean",
    asw = "Americas Black USA",
    clm = "Americas Columbian",
    mxl = "Americas Mexican",
    pel = "Americas Peruvian",
    pur = "Americas Puerto Rican",
    ceu = "Americas White USA",
    chb = "Asia Han Beijing",
    cdx = "Asia Dai",
    jpt = "Asia Japanese",
    khv = "Asia Kinh",
    chs = "Asia Han Southern",
    gbr = "Europe British",
    fin = "Europe Finn",
    ibs = "Europe Spanish",
    tsi = "Europe Tuscan",
    beb = "SW Asia Bengali",
    gih = "SW Asia Gujarati",
    pjl = "SW Asia Punjabi",
    stu = "SW Asia Tamil",
    itu = "SW Asia Telugu",
    gwd = "Africa Gambian",
    lwk = "Africa Luhya",
    msl = "Africa Mende",
    yri = "Africa Yoruban",
    esn = "Africa Esan",
    neandertalAltai = "Ancient Neandertal Altai",
    neandertalChagyrskaya = "Ancient Neandertal Chagyrskaya",
    neandertalVindija = "Ancient Neandertal Vindija",
    denisova = "Ancient Denisovan",
    ustishim = "Ancient Ust Ishim"
)

########### DNA
# read data
genvar = readRDS("pervar_dna_heatmaps.rds")

# add revised population names to data
genvar$population_abbr <- genvar$population
genvar$population <- pop_names[genvar$population]
genvar = relocate(genvar, population_abbr, .before = population)

#make factors
genvar$gene = factor(genvar$gene)

genvar$population <- factor(genvar$population, levels = c(
    "Ancient Neandertal Altai", "Ancient Neandertal Chagyrskaya",
    "Ancient Neandertal Vindija", "Ancient Denisovan",
    "Ancient Ust Ishim", "Africa Esan", "Africa Gambian",
    "Africa Luhya", "Africa Mende", "Africa Yoruban",
    "Americas Black USA", "Americas Black Caribbean",
    "Europe British", "Europe Finn", "Europe Spanish",
    "Europe Tuscan", "Americas White USA", "SW Asia Bengali",
    "SW Asia Gujarati", "SW Asia Punjabi", "SW Asia Tamil",
    "SW Asia Telugu", "Asia Dai", "Asia Han Beijing",
    "Asia Han Southern", "Asia Japanese",
    "Asia Kinh", "Americas Columbian",
    "Americas Mexican", "Americas Puerto Rican", "Americas Peruvian"))

# order genes on y-axis
genvar <- genvar |>
    arrange(desc(chr), desc(start)) |>
    mutate(gene = factor(gene, unique(gene)))

########### AA
# read data
aavar = readRDS("pervar_aa_heatmaps.rds")

# add revised population names to data
aavar$population_abbr <- aavar$population
aavar$population <- pop_names[aavar$population]
aavar = relocate(aavar, population_abbr, .before = population)

# make factors
aavar$population <- factor(aavar$population, levels = c(
    "Ancient Neandertal Altai", "Ancient Neandertal Chagyrskaya",
    "Ancient Neandertal Vindija", "Ancient Denisovan",
    "Ancient Ust Ishim", "Africa Esan", "Africa Gambian",
    "Africa Luhya", "Africa Mende", "Africa Yoruban",
    "Americas Black USA", "Americas Black Caribbean",
    "Europe British", "Europe Finn", "Europe Spanish",
    "Europe Tuscan", "Americas White USA", "SW Asia Bengali",
    "SW Asia Gujarati", "SW Asia Punjabi", "SW Asia Tamil",
    "SW Asia Telugu", "Asia Dai", "Asia Han Beijing",
    "Asia Han Southern", "Asia Japanese",
    "Asia Kinh", "Americas Columbian",
    "Americas Mexican", "Americas Puerto Rican", "Americas Peruvian"))

# order genes on y-axis
aavar <- aavar |>
    arrange(desc(chr), desc(start)) |>
    mutate(gene = factor(gene, unique(gene)))

########### Plots
# genetic variation heatmap
a= ggplot(genvar, aes(x=population, y=gene, fill=pervar)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    xlab("") +
    ylab("Percent of Gene with DNA Variants") +
    geom_tile(color='White', size=0.1) +
    scale_fill_viridis() +
    coord_equal() +
    theme(legend.title = element_blank())+
    theme(text=element_text(size=35,  family="sans"))

# amino acid variation heatmap
b = ggplot(aavar, aes(x=population, y=gene, fill=pervar)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    xlab("") +
    ylab("Percent of Gene with DNA Variants") +
    geom_tile(color='White', size=0.1) +
    scale_fill_viridis() +
    coord_equal() +
    theme(legend.title = element_blank())+
    theme(text=element_text(size=35,  family="sans"))

#arrange plots and save
png("revised-f1-heatmap.png", width=45, height=20,
    res=400, units = "in")
ggarrange(a, ggplot() + theme_void(), b,
          nrow=1, widths = c(4, -1, 4),
          labels = c("A", "", "B"),
          font.label = list(size = 35, face = "bold"))
dev.off()



######################## FUNCTIONAL VARIATION
########### Row 1
den2016 = read_excel("oldnoses-orfuncdata-2016-all.xlsx", sheet = "denisova-2016")
den2016$Odor = factor(den2016$Odor)
den2016$OR = factor(den2016$OR)
den2016$Human_log = log(den2016$Human)
den2016$Denisova_log = log(den2016$Denisova)

nean2016 = read_excel("oldnoses-orfuncdata-2016-all.xlsx", sheet = "neandertal-2016")
nean2016$Odor = factor(nean2016$Odor)
nean2016$OR = factor(nean2016$OR)
nean2016$Human_log = log(nean2016$Human)
nean2016$Neandertal_log = log(nean2016$Neandertal)

# Panel A
den2016_graph = ggplot(den2016, aes(x=Human_log, y=Denisova_log))+
    geom_point(aes(colour = OR)) +
    geom_smooth(method="lm", se= FALSE, col="black") +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             method = "pearson",
             p.accuracy = 0.001, r.accuracy = 0.01,
             label.x = -1, label.y = 2.4, size=6) +
    stat_regline_equation(label.x = -1, label.y = 2.1, size=6) +
    scale_color_viridis_d() +
    xlim(-1, 2.5) +
    ylim(-1, 2.5) +
    labs(x = "Human Response(ln)",
         y = "Denisovan Response(ln)") +
    theme(text = element_text(family = "sans", face = "bold"))+
    theme_pubr(base_size = 20,legend= c(.9,.4))

# Panel B
nean2016_graph = ggplot(nean2016, aes(x=Human_log, y=Neandertal_log))+
    geom_point(aes(colour = OR)) +
    geom_smooth(method="lm", se= FALSE, col="black") +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             method = "pearson",
             p.accuracy = 0.001, r.accuracy = 0.01,
             label.x = -1, label.y = 2.4, size=6) +
    stat_regline_equation(label.x = -1, label.y = 2.1, size=6) +
    xlim(-1, 2.5) +
    ylim(-1, 2.5) +
    labs(x = "Human Response(ln)",
         y = "Neandertal Response(ln)") +
    scale_color_viridis_d()+
    theme(text = element_text(family = "sans", face = "bold"))+
    theme_pubr(base_size = 20,legend= c(.8,.3))

########### Row 2
# Panel C
vcf2016long = read_excel("oldnoses-orfuncdata-2016-all.xlsx", sheet = "VCF2016-long")
vcf2016long$OR = factor(vcf2016long$OR)
vcf2016long$Lineage3 = ordered(vcf2016long$Lineage,
                               levels = c("Human", "Denisova", "Neandertal"))
vcf2016long$Lineage3 <- gsub("Denisova", "Denisovan",vcf2016long$Lineage3)
vcf2016long$Response_log = log(vcf2016long$Response)

boxplot2016_graph=ggplot(vcf2016long, aes(x=Lineage3, y=Response_log, fill=Lineage3)) +
    geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=2)+
    geom_jitter(shape=16, position=position_jitter(0.2))+
    labs(x="", y = "Response(ln)") +
    scale_fill_brewer(palette="PuOr")+
    theme(text = element_text(family = "sans", face = "bold"))+
    theme_pubr(base_size = 20, legend = "none")

#Panel D
byOR2016_graph = ggplot(vcf2016long, aes(Response_log, Lineage3)) +
    geom_jitter(aes(col=OR)) +
    geom_smooth(aes(col=OR), method="lm", se=F)+
    scale_color_viridis_d() +
    labs(x = "Response by OR",
         y = "") +
    theme(text = element_text(family = "sans", face = "bold"))+
    theme_pubr(base_size = 20)

########### PLOT
png("revised-f2-funcvar.png", width=45, height=20,
    res=400, units = "in")
row1 <- ggarrange(nean2016_graph, boxplot2016_graph,
                  labels = c("A", "C"), nrow = 1,
                  font.label = list(size = 35, face = "bold"))

row2 <- ggarrange(den2016_graph, byOR2016_graph,
                  labels = c("B", "D"), nrow = 1,
                  font.label = list(size = 35, face = "bold"))

ggarrange(row1, row2, nrow = 2)
dev.off()





######################## EVOLUTIONARY PLOTS
########### Read aa
or1a1=readAAStringSet(file = "fastas-aa-clean-noref/or1a1-aa-all.fa", format = "fasta")
or2a25=readAAStringSet(file = "fastas-aa-clean-noref/or2a25-aa-all.fa", format = "fasta")
or2c1=readAAStringSet(file = "fastas-aa-clean-noref/or2c1-aa-all.fa", format = "fasta")
or2j2=readAAStringSet(file = "fastas-aa-clean-noref/or2j2-aa-all.fa", format = "fasta")
or2j3=readAAStringSet(file = "fastas-aa-clean-noref/or2j3-aa-all.fa", format = "fasta")
or4e2=readAAStringSet(file = "fastas-aa-clean-noref/or4e2-aa-all.fa", format = "fasta")
or4q3=readAAStringSet(file = "fastas-aa-clean-noref/or4q3-aa-all.fa", format = "fasta")
or5a1=readAAStringSet(file = "fastas-aa-clean-noref/or5a1-aa-all.fa", format = "fasta")
or5an1=readAAStringSet(file = "fastas-aa-clean-noref/or5an1-aa-all.fa", format = "fasta")
or5k1=readAAStringSet(file = "fastas-aa-clean-noref/or5k1-aa-all.fa", format = "fasta")
or8k3=readAAStringSet(file = "fastas-aa-clean-noref/or8k3-aa-all.fa", format = "fasta")
or10g4=readAAStringSet(file = "fastas-aa-clean-noref/or10g4-aa-all.fa", format = "fasta")
or51e1=readAAStringSet(file = "fastas-aa-clean-noref/or51e1-aa-all.fa", format = "fasta")
or51l1=readAAStringSet(file = "fastas-aa-clean-noref/or51l1-aa-all.fa", format = "fasta")
or1c1=readAAStringSet(file = "fastas-aa-clean-noref/or1c1-aa-all.fa", format = "fasta")
or2b11=readAAStringSet(file = "fastas-aa-clean-noref/or2b11-aa-all.fa", format = "fasta")
or2w1=readAAStringSet(file = "fastas-aa-clean-noref/or2w1-aa-all.fa", format = "fasta")
or5p3=readAAStringSet(file = "fastas-aa-clean-noref/or5p3-aa-all.fa", format = "fasta")
or6p1=readAAStringSet(file = "fastas-aa-clean-noref/or6p1-aa-all.fa", format = "fasta")
or7c1=readAAStringSet(file = "fastas-aa-clean-noref/or7c1-aa-all.fa", format = "fasta")
or7d4=readAAStringSet(file = "fastas-aa-clean-noref/or7d4-aa-all.fa", format = "fasta")
or8b3=readAAStringSet(file = "fastas-aa-clean-noref/or8b3-aa-all.fa", format = "fasta")
or8d1=readAAStringSet(file = "fastas-aa-clean-noref/or8d1-aa-all.fa", format = "fasta")
or10a6=readAAStringSet(file = "fastas-aa-clean-noref/or10a6-aa-all.fa", format = "fasta")
or10g3=readAAStringSet(file = "fastas-aa-clean-noref/or10g3-aa-all.fa", format = "fasta")
or10g7=readAAStringSet(file = "fastas-aa-clean-noref/or10g7-aa-all.fa", format = "fasta")
or10j5=readAAStringSet(file = "fastas-aa-clean-noref/or10j5-aa-all.fa", format = "fasta")
or11a1=readAAStringSet(file = "fastas-aa-clean-noref/or11a1-aa-all.fa", format = "fasta")
or51e2=readAAStringSet(file = "fastas-aa-clean-noref/or51e2-aa-all.fa", format = "fasta")
or56a4=readAAStringSet(file = "fastas-aa-clean-noref/or56a4-aa-all.fa", format = "fasta")


# concatenate strings
uber.msa=AAStringSet(paste0(or1a1,or2a25,or2c1,or2j2,or2j3,or4e2,
                            or4q3,or5a1,or5an1,or5k1,or8k3,or10g4,
                            or51e1,or51l1,or1c1,or2b11,or2w1,or5p3,
                            or6p1,or7c1,or7d4,or8b3,or8d1,or10a6,
                            or10g3,or10g7,or10j5,or11a1,or51e2,or56a4))

# label strings with names
names(uber.msa)= names(or1a1)

# save strings
writeXStringSet(uber.msa,"uber-msa.fa")

########### Cladogram

distalign.nj.ubertree = nj(distalign.ubertree)
png(filename="revised-cladogram.png",
    units="in",
    width=10,
    height=10,
    res=300)
plot(distalign.nj.ubertree, main = "", family = "sans", font=3, ps=12,cex=1.5)
dev.off()



########### Phylogram
tree.dendro = read.tree(text = "
                       (Hominids,
                       (Hominin,
                       (Hominina,
                       (Homo,
                       (Neandertal, Denisovan, Human,
                       (Hybrids))))))
                       ;")
png(filename="revised-phylotree.png",
    units="in",
    width=10,
    height=10,
    res=300)
plot(tree.dendro, use.edge.length = FALSE, family = "sans", font=3, ps=12,cex=1.5)
nodelabels("14-12 mya", 9, frame = "r", bg = "beige", adj = 0, family = "sans", font=3, ps=12,cex=1.5)
nodelabels("10-8 mya", 10, frame = "r", bg = "beige", adj = 0, family = "sans", font=3, ps=12,cex=1.5)
nodelabels("8-4 mya", 11, frame = "r", bg = "beige", adj = 0, family = "sans", font=3, ps=12,cex=1.5)
nodelabels("2.5-2 mya", 12, frame = "r", bg = "beige", adj = 0, family = "sans", font=3, ps=12,cex=1.5)
nodelabels("0.8 mya", 13, frame = "r", bg = "beige", adj = 0, family = "sans", font=3, ps=12,cex=1.5)
dev.off()



######################## TIDY
rm(list=ls())
gc()
