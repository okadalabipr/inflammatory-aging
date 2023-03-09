library(ggplot2)
library(dplyr)
library(NSM3)

setwd("inframmatory-aging/bioinfomatics/")
dir.create("ATACseq_analysis/")

#Analysis of chromatin accessibility in open chromatin regions near each cluster (fig3B)----
getdata <- function(n){
  d <- read.csv(n, sep = "\t")
  d$chr <- d$Merged_Peak
  d$start <- d$Merged_Peak
  d$end <- d$Merged_Peak
  d$chr <- gsub("_.*","",d$chr)
  d$start <- gsub("chr1_|chr2_|chr3_|chr4_|chr5_|chr6_|chr7_|chr8_|chr9_|chr10_|chr11_|chr12_|chr13_|chr14_|chr15_|chr16_|chr17_|chr18_|chr19_|chr20_|chr21_|chr22_|chrX_|chrY_", "", d$start)
  d$start <- gsub("_.*", "", d$start)
  d$end <- gsub("c.*_", "", d$end)
  d <- dplyr::filter(d, abs(d$Distance.to.TSS) < 5000)
  return(d)
}

#load annotation file
cs <- getdata("path to your Ctrl and IkBaKD peak comparison file obtained from GENCODE ATAC-seq pipeline (https://github.com/ENCODE-DCC/atac-seq-pipeline)")
cct <- getdata("path to your Ctrl and TNF peak comparison file obtained from GENCODE ATAC-seq pipeline (https://github.com/ENCODE-DCC/atac-seq-pipeline)")
cst <- getdata("path to your Ctrl and IkBaKDTNF peak comparison file obtained from GENCODE ATAC-seq pipeline (https://github.com/ENCODE-DCC/atac-seq-pipeline)")

#load genelist
pam1 <- read.csv("pam1.csv", col.names = "Gene.Name") #Output of 1_MCF7_RNAseq_analysis.R
pam2 <- read.csv("pam2.csv", col.names = "Gene.Name")
pam3 <- read.csv("pam3.csv", col.names = "Gene.Name")
pam4 <- read.csv("pam4.csv", col.names = "Gene.Name")

chrom_ac <- function(x,y){
  pam1_ac <- dplyr::inner_join(x, pam1, by = "Gene.Name")
  pam2_ac <- dplyr::inner_join(x, pam2, by = "Gene.Name")
  pam3_ac <- dplyr::inner_join(x, pam3, by = "Gene.Name")
  pam4_ac <- dplyr::inner_join(x, pam4, by = "Gene.Name")
  
  vec_pam1 <- pam1_ac$log2FoldChange
  vec_pam2 <- pam2_ac$log2FoldChange
  vec_pam3 <- pam3_ac$log2FoldChange
  vec_pam4 <- pam4_ac$log2FoldChange
  
  df <- data.frame(
    label = c(rep("cluster1", length(vec_pam1)), rep("cluster2", length(vec_pam2)),
              rep("cluster3", length(vec_pam3)), rep("cluster4", length(vec_pam4))),
    FC = c(vec_pam1, vec_pam2, vec_pam3, vec_pam4)
  )
  
  df$ID <- rep(y, nrow(df))
  
  return(df)
  
}

cs_df <- chrom_ac(cs, 1)
cct_df <- chrom_ac(cct, 2)
cst_df <- chrom_ac(cst, 3)

all <- rbind(cs_df, cct_df, cst_df)

all_s <- group_by(all, ID, label) %>% 
  summarise_at(vars(FC), list(mean = ~mean(.), 
                              sd = ~sd(.), 
                              se = ~sd(.)/sqrt(length(.))))

all_s$ID <- as.character(all_s$ID)

#Bar plot of chromatin accessibility(fig.3B)
g <- ggplot(all_s, aes(x = label, y = mean,group = ID, fill = ID)) +
  geom_bar(stat = "identity",colour = "black",lwd = 0.7, width = 0.9, position = "dodge") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),position = position_dodge(0.9), width = .2)+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#CCCCCC", "#CCCCCC","#B22222" ))+
  scale_y_continuous(breaks = c(-0.1, 0, 0.1, 0.2), limits = c(min(all_s$mean)-max(all_s$se), max(all_s$mean)+max(all_s$se)) )+
  scale_x_discrete(limit=c("cluster1","cluster2","cluster3","cluster4"))+
  labs(x="", y="")+
  theme(axis.text.x=element_text(size=15, face = "bold", colour = "black",angle = 45,hjust = 1),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1,face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
plot(g)
ggsave(file = "all_chrom_accessibility_5000.png", plot = g, dpi = 600, width = 6, height = 7)

#Steel-Dwass test 
all$label <- paste0(all$label,"-",all$ID)
all$FC <- as.factor(all$FC)
all$label <- as.factor(all$label)
test <- pSDCFlig(all$FC, all$label, method="Asymptotic")
test <- data.frame(label = test$labels , p.val = test$p.val)

#Create MAplot(fig.S14B)
getdata <- function(n){
  d <- read.csv(n, sep = "\t")
  d$chr <- d$Merged_Peak
  d$start <- d$Merged_Peak
  d$end <- d$Merged_Peak
  d$chr <- gsub("_.*", "", d$chr)
  d$start <- gsub("chr1_|chr2_|chr3_|chr4_|chr5_|chr6_|chr7_|chr8_|chr9_|chr10_|chr11_|chr12_|chr13_|chr14_|chr15_|chr16_|chr17_|chr18_|chr19_|chr20_|chr21_|chr22_|chrX_|chrY_", "", d$start)
  d$start <- gsub("_.*", "", d$start)
  d$end <- gsub("c.*_", "", d$end)
  return(d)
}

cst <- getdata("path to your Ctrl and IkBaKDTNF peak comparison file obtained from GENCODE ATAC-seq pipeline (https://github.com/ENCODE-DCC/atac-seq-pipeline)")
cst$log2baseMean <- log2(cst$baseMean)

ttmp <- rep(1, nrow(cst))
ttmp[cst$log2FoldChange > 0 & cst$padj < 0.05] <- 2
cst$threshold <- as.factor(ttmp)

cols <- c("#808080", "#BC0000")

g <- ggplot(data=cst, aes(x=log2baseMean, y=log2FoldChange, colour=threshold)) +
  geom_point(size=1.75, alpha = 0.3) +
  theme_classic()+
  theme(legend.position = "none") +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(breaks = seq(floor(max(cst$log2baseMean)), floor(min(cst$log2baseMean)), by = -2))+
  scale_y_continuous(breaks = seq(floor(max(cst$log2FoldChange)), floor(min(cst$log2FoldChange)), by = -2))+
  scale_fill_manual(values=cols) + 
  scale_color_manual(values=cols) +
  theme(axis.text.x=element_text(size=20, face = "bold", colour = "black"),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1, face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
plot(g)
ggsave(file = "MAplot_opencst.png", plot = g, dpi = 600, width = 8, height = 7)


#Extraction of open chromatin regions under IkBaKD+TNF condition within 10,000 bp from TSS
cst_o <- dplyr::filter(cst,cst$log2FoldChange > 0, cst$padj < 0.05, abs(cst$Distance.to.TSS) < 10000) 

pam1 <- read.csv("pam1.csv", col.names = "Gene.Name") 
pam2 <- read.csv("pam2.csv", col.names = "Gene.Name")
pam3 <- read.csv("pam3.csv", col.names = "Gene.Name")
pam4 <- read.csv("pam4.csv", col.names = "Gene.Name")

pam1 <- dplyr::inner_join(cst_o,pam1, by = "Gene.Name")
write.table(pam1[,20:22], "ATACseq_analysis/pam1_near.bed", col.names = F, row.names = F, quote = F, sep = "\t")
pam2 <- dplyr::inner_join(cst_o,pam2, by = "Gene.Name")
write.table(pam2[,20:22], "ATACseq_analysis/pam2_near.bed", col.names = F, row.names = F, quote = F, sep = "\t")
pam3 <- dplyr::inner_join(cst_o,pam3, by = "Gene.Name")
write.table(pam3[,20:22], "ATACseq_analysis/pam3_near.bed", col.names = F, row.names = F, quote = F, sep = "\t")
pam4 <- dplyr::inner_join(cst_o,pam4, by = "Gene.Name")
write.table(pam4[,20:22], "ATACseq_analysis/pam4_near.bed", col.names = F, row.names = F, quote = F, sep = "\t")
