# This script creates plots used in Supplemental Information Figure 2

library(readxl); library(ggplot2); library(ggpubr)

#ROW 1
alt2013 = read_excel("oldnoses-orfuncdata-2013-all.xlsx", sheet = "alt-2013")
alt2013$Odor = factor(alt2013$Odor)
alt2013$OR = factor(alt2013$OR)
alt2013$Human_log = log(alt2013$Human)
alt2013$Altai_log = log(alt2013$Altai)
den2013 = read_excel("oldnoses-orfuncdata-2013-all.xlsx", sheet = "den-2013")
den2013$Odor = factor(den2013$Odor)
den2013$OR = factor(den2013$OR)
den2013$Human_log = log(den2013$Human)
den2013$Denisova_log = log(den2013$Denisova)
ust2013 = read_excel("oldnoses-orfuncdata-2013-all.xlsx", sheet = "ust-2013")
ust2013$Odor = factor(ust2013$Odor)
ust2013$OR = factor(ust2013$OR)
ust2013$Human_log = log(ust2013$Human)
ust2013$Ust_Ishim_log = log(ust2013$Ust_Ishim)

#Panel A
ust2013.graph = ggplot(ust2013, aes(x=Human_log, y=Ust_Ishim_log))+
  geom_point(aes(colour=OR)) +
  geom_smooth(method="lm", se= FALSE, col="black") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           method = "pearson",
           p.accuracy = 0.001, r.accuracy = 0.01,
           label.x = -3, label.y = 3.5, size = 6) +
  stat_regline_equation(label.x = -3, label.y = 2.8, size = 6) +
  labs(x = "Human Response(ln)",
       y = "Ust'-Ishim Response(ln)") +
  scale_color_viridis_d() +
  xlim(-3,3.5)+
  ylim(-3,3.5)+
  theme_pubr(base_size = 20,legend= "right")

#Panel B
den2013.graph = ggplot(den2013, aes(x=Human_log, y=Denisova_log))+
  geom_point(aes(colour = OR)) +
  geom_smooth(method="lm", se= FALSE, col="black") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           method = "pearson",
           p.accuracy = 0.001, r.accuracy = 0.01,
           label.x = -3, label.y = 3.5, size = 6) +
  stat_regline_equation(label.x = -3, label.y = 2.8, size = 6) +
  labs(x = "Human Response(ln)",
       y = "Denisova Response(ln)") +
  scale_color_viridis_d() +
  xlim(-3,3.5)+
  ylim(-3,3.5)+
  theme_pubr(base_size = 20,legend= "right")

#Panel C
alt2013.graph = ggplot(alt2013, aes(x=Human_log, y=Altai_log))+
  geom_point(aes(colour = OR)) +
  geom_smooth(method="lm", se= FALSE, col="black") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           method = "pearson",
           p.accuracy = 0.001, r.accuracy = 0.01,
           label.x = -3, label.y = 3.5, size = 6) +
  stat_regline_equation(label.x = -3, label.y = 2.8, size = 6) +
  labs(x = "Human Response(ln)",
       y = "Altai Response(ln)") +
  scale_color_viridis_d() +
  xlim(-3,3.5)+
  ylim(-3,3.5)+
  theme_pubr(base_size = 20,legend= "right")

#ROW 2
vcf2013long = read_excel("oldnoses-orfuncdata-2013-all.xlsx", sheet = "VCF2013-long")
vcf2013long$OR = factor(vcf2013long$OR)
vcf2013long$Lineage = ordered(vcf2013long$Lineage,
  levels = c("Human", "Ust_Ishim", "Denisova", "Altai"))
vcf2013long$Lineage <- gsub("Denisova", "Denisovan",vcf2013long$Lineage)
vcf2013long$Response_log = log(vcf2013long$Response)

#Panel D
boxplot2013.graph1 = ggplot(vcf2013long, aes(x=Lineage, y=Response_log, fill=Lineage)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=2)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(x="", y = "Response(ln)") +
  scale_fill_brewer(palette="PuOr")+
  theme_pubr(base_size = 20, legend = "none")

boxplot2013.graph <- boxplot2013.graph1 +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

#Panel E
byOR2013.graph = ggplot(vcf2013long, aes(Response_log, Lineage)) +
  geom_jitter(aes(col=OR)) +
  geom_smooth(aes(col=OR), method="lm", se=F)+
  scale_color_viridis_d() +
  labs(x = "Response(ln) by OR",
       y = "") +
  scale_y_discrete(labels = c("ALT", "DEN","UST", "HUM"))+
  theme_pubr(base_size = 20, legend="right")

#PLOT
library(cowplot)
tiff("Rscript-oldnoses-2016-func var figure.tiff", width = 20,
     height=12, res=300, units="in",pointsize=12, type="cairo")
row1 = plot_grid(alt2013.graph, den2013.graph, ust2013.graph,
  labels = c('a', 'b', 'd'),nrow=1)
row2 = plot_grid(boxplot2013.graph, byOR2013.graph,
  labels = c('d', 'e'), nrow=1)
plot_grid(row1, row2, nrow=2, rel_heights = c(1,1))
dev.off()
