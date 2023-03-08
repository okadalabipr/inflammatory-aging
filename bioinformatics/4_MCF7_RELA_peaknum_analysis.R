library(dplyr)
library(ggplot2)
library(rtracklayer)
library(NSM3)

setwd("NFkB_aging_paper/")

#List of genes registered in gencode----
gtf <- readGFF("ref_file/gencode.v40.annotation.gtf")
gtf_gene <- subset(gtf, gtf$type == "gene" & gtf$gene_type == "protein_coding")
gtf_gene <- gtf_gene$gene_name %>% as.data.frame() %>% dplyr::rename("Gene.Name" = 1) %>% distinct(Gene.Name, .keep_all=TRUE)
rm(gtf)

#Loading RELA ChIP peak files
rela_siIkBa <- read.table("source_file/siIkBa_120_R2_peaks.annotatePeaks.txt", sep = "\t") 
colnames(rela_siIkBa) <- rela_siIkBa[1,] 
rela_siIkBa <- rela_siIkBa[-1,]
rela_siIkBa <- dplyr::select(rela_siIkBa, "Gene Name", "Distance to TSS")
colnames(rela_siIkBa) <- c("Gene.Name","Distance.to.TSS")
rela_siIkBa <- dplyr::inner_join(gtf_gene, rela_siIkBa,by = "Gene.Name")
rela_siIkBa$Distance.to.TSS <- as.numeric(rela_siIkBa$Distance.to.TSS)
rela_siIkBa_5000 <- dplyr::filter(rela_siIkBa, abs(rela_siIkBa$Distance.to.TSS) < 5000)

#load genelist
rela_target <- rela_siIkBa_5000$Gene.Name %>% unique() %>% as.data.frame() %>% dplyr::rename("Gene.Name" = 1) 
pam1 <- read.csv("pam1.csv", col.names = "Gene.Name") %>% dplyr::inner_join(rela_target, by = "Gene.Name")
pam2 <- read.csv("pam2.csv", col.names = "Gene.Name") %>% dplyr::inner_join(rela_target, by = "Gene.Name")
pam3 <- read.csv("pam3.csv", col.names = "Gene.Name") %>% dplyr::inner_join(rela_target, by = "Gene.Name")
pam4 <- read.csv("pam4.csv", col.names = "Gene.Name") %>% dplyr::inner_join(rela_target, by = "Gene.Name")

#count of RELA peaks
rela_peak_num <- function(x){
  y <- dplyr::inner_join(x, rela_siIkBa, by = "Gene.Name")
  ttmp <- rep(0, nrow(x))
  for (i in 1:nrow(x)){
    a <- dplyr::filter(y, y$Gene.Name == x[i,])
    ttmp[i] <- as.numeric(nrow(a))
  }
  x$num <- ttmp
  x <- dplyr::filter(x, x$num > 0)
  return(x)
}

pam1 <- rela_peak_num(pam1)
pam2 <- rela_peak_num(pam2)
pam3 <- rela_peak_num(pam3)
pam4 <- rela_peak_num(pam4)
rela_target <- rela_peak_num(rela_target)

df <- data.frame(
  label = c(rep("RELA_target", length(rela_target$num)),
            rep("cluster1", length(pam1$num)), rep("cluster2", length(pam2$num)),
            rep("cluster3", length(pam3$num)), rep("cluster4", length(pam4$num))),
  num = c(rela_target$num, pam1$num, pam2$num, pam3$num, pam4$num)
)

df2 <- group_by(df, label) %>% 
  summarise_all(list(mean = mean, 
                     sd = sd, 
                     se = ~sd/sqrt(length(.))))

df2$factor <- factor(df2$label, levels = c("RELA_target", "cluster1", "cluster2", "cluster3", "cluster4"))

#Bar plot of average RELA peaks 
g <- ggplot(df2, aes(x=factor, y = mean, fill=label))+
  geom_bar(stat = "identity", colour = "black", lwd = 0.7, width = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .2)+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#51A39D", "#814374", "#B7695C", "#805f2b", "#d6ecfa"))+
  scale_y_continuous(breaks = seq(1.5, 0, by = -0.5), limits = c(0, round(max(df2$mean)+max(df2$se), digits = 2)), expand = c(0, 0))+
  labs(x="", y="")+
  theme(axis.text.x=element_text(size=20, face = "bold", colour = "black",angle = 45, hjust = 1),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1, face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
plot(g)

ggsave(file = "num_rela_peaks.png", plot = g, dpi = 600, width = 3.6, height = 4.8)

#Steel-Dwass test 
df$num <- as.factor(df$num)
df$label <- as.factor(df$label)
test <- pSDCFlig(df$num, df$label, method="Asymptotic")
test <- data.frame(label = test$labels, p.val = test$p.val)

