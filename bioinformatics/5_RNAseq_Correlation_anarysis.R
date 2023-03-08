library(cummeRbund)
library(biomaRt)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("NFkB_aging_paper/")

# Basic function to convert mouse to human gene names----
mouse_human_genes <- read.csv("ref_file/MGI_mouse_human_gene_Correspondence_table.csv")
convert_mouse_to_human <- function(gene_list){
  df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
  for(gene in gene_list2){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      temp <- c(gene,human_genes)
      df <- rbind(df,temp)
    }
  }
  colnames(df) <- c("mice_gene_name","human_gene_name")
  return(df)
}

#Correspondence table between Ensembl geneID and Uniprot symbol----
gtf <- readGFF("ref_file/gencode.v40.annotation.gtf")
gtf_gene <- subset(gtf, gtf$type == "gene" & gtf$gene_type == "protein_coding")
gtf_gene <- gtf_gene[, c("gene_id","gene_name")]
gtf_gene <- gtf_gene %>% distinct(gene_name, .keep_all=TRUE)
colnames(gtf_gene) <- c("Geneid", "gene_name")
rm(gtf)

# FC calculation for each gene in MCF7 (Ctrl vs IkBaKD+TNF)------
mcf7 <- read.table("source_file/salmon_MCF7_48h_gene_count.tsv", 
                   header=T, stringsAsFactors=F,sep="\t")
mcf7 <- dplyr::select(mcf7, "gene_name", "MCF7_CONTsi_1", "MCF7_CONTsi_2", "MCF7_CONTsi_3",
                      "MCF7_IKBAsi_TNF_1", "MCF7_IKBAsi_TNF_2", "MCF7_IKBAsi_TNF_3")

header <- c("gene_name","ctrl_rep1","ctrl_rep2","ctrl_rep3",
            "siIkBa_TNF_rep1","siIkBa_TNF_rep2","siIkBa_TNF_rep3")

colnames(mcf7) <- header
mcf7 <- dplyr::inner_join(mcf7,gtf_gene,by = "gene_name") %>% dplyr::select(-"Geneid")

mcf7 <- dplyr::filter(mcf7, (ctrl_rep1 >= 10 & ctrl_rep2 >= 10 & ctrl_rep3 >= 10)
                      | (siIkBa_TNF_rep1 >= 10 & siIkBa_TNF_rep2 >= 10 & siIkBa_TNF_rep3 >= 10))

mcf7 <- mcf7[!duplicated(mcf7$gene_name), ]
row.names(mcf7) <- mcf7$gene_name
mcf7 <- dplyr::select(mcf7, -"gene_name")

group <- data.frame(con = factor(c(rep("siCtrl", 3), rep("siIkBa_TNF", 3))))

dds <- DESeqDataSetFromMatrix(countData = round(mcf7), colData = group, design = ~ con)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds)

res <- results(dds, contrast = c("con", "siIkBa_TNF", "siCtrl"))
res$gene_id <- row.names(res)
res <- res[!is.na(res$padj), ]
mcf7 <- as.data.frame(res)
mcf7 <- mcf7[, c(7, 2, 6)]
colnames(mcf7) <- c("human_gene_name", "mcf7_cst_log2fc", "mcf7_cst_padj")


# FC calculation for each gene in mouse heart tissue (young vs aged)-----
cnt <- read.table("source_file/salmon_mouse_heart_gene_count.tsv", 
                  header=T, stringsAsFactors=F,sep="\t")
cnt <- dplyr::select(cnt, -"gene_id")

header <- c("gene_name","old_rep1","old_rep2","old_rep3","old_rep4","old_rep5",
            "young_rep1","young_rep2","young_rep3","young_rep4","young_rep5")
colnames(cnt) <- header

cnt <- dplyr::filter(cnt, (old_rep1 >= 5 & old_rep2 >= 5 & old_rep3 >= 5 & old_rep4 >= 5 & old_rep5 >= 5) 
                     | (young_rep1 >= 5 & young_rep2 >= 5 & young_rep3 >= 5 & young_rep4 >= 5 & young_rep5 >= 5))

cnt <- cnt[!duplicated(cnt$gene_name),]
row.names(cnt) <- cnt$gene_name
cnt <- dplyr::select(cnt, -"gene_name")

group <- data.frame(con = factor(c(rep("old", 5), rep("young", 5))))

dds <- DESeqDataSetFromMatrix(countData = round(cnt), colData = group, design = ~ con)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds)

res <- results(dds, contrast = c("con", "old", "young"))
res$gene_id <- row.names(res)
res <- res[!is.na(res$padj), ]
mouse_ya <- as.data.frame(res)

mouse_ya$mice_gene_name <- row.names(mouse_ya)
gene_list <- mouse_ya$mice_gene_name
gene_list2 <- mouse_ya$mice_gene_name

gene_list <- convert_mouse_to_human(gene_list2)

mouse_ya <- dplyr::inner_join(gene_list, mouse_ya, by = "mice_gene_name") 
mouse_ya <- mouse_ya[, c(2, 4, 8)] %>% dplyr::rename("human_gene_name" = 1, "mice_log2fc" = 2, "mice_padj" = 3)

#Integration of MCF7 and mouse heart tissue RNA-seq data
h_m <- dplyr::inner_join(mcf7, mouse_ya, by = "human_gene_name") %>% dplyr::select(c("human_gene_name", "mcf7_cst_log2fc", "mice_log2fc" ))

h_m_list <- h_m$human_gene_name

sp <- read.csv("ref_file/senescence_related_gene.csv") #read senescense gene list 

h_m_list[!h_m_list %in% sp$x] <- ""
h_m$gene_name <- h_m_list

pam1 <- read.csv("pam1.csv",col.names = "human_gene_name") #load genelist
pam2 <- read.csv("pam2.csv",col.names = "human_gene_name")
pam3 <- read.csv("pam3.csv",col.names = "human_gene_name")
pam4 <- read.csv("pam4.csv",col.names = "human_gene_name")

pam1 <- dplyr::inner_join(pam1,h_m,by = "human_gene_name")
pam2 <- dplyr::inner_join(pam2,h_m,by = "human_gene_name")
pam3 <- dplyr::inner_join(pam3,h_m,by = "human_gene_name")
pam4 <- dplyr::inner_join(pam4,h_m,by = "human_gene_name")

all <- rbind(pam1,pam2,pam3,pam4)

#Correlation plots for each cluster (fig.3G)
g_pam1 <- ggplot(pam1, aes(x = mice_log2fc,y = mcf7_cst_log2fc,label = gene_name))+
  geom_point(alpha = 0.7,size = 3,colour = "#51A39D")+
  xlab(NULL)+
  ylab(NULL)+
  geom_text_repel(max.overlaps = 500,size = 7)+
  scale_x_continuous(breaks = c(-2,0,2,4), limits = c(round(min(all$mice_log2fc), digits = 2) - 0.01, round(max(all$mice_log2fc), digits = 2)+0.01))+
  scale_y_continuous(breaks = c(0,4,8,12), limits = c(round(min(all$mcf7_cst_log2fc), digits = 2) - 0.01, round(max(all$mcf7_cst_log2fc), digits = 2)+0.01))+
  theme_classic()+
  theme(axis.text.x=element_text(size=20, face = "bold", colour = "black"),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1,face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

plot(g_pam1)

cor.test(pam1$mice_log2fc, pam1$mcf7_cst_log2fc, method="spearman", exact=FALSE)
ggsave(file = "cluster1_Mice.png", plot = g_pam1, dpi = 600, width = 7, height = 4.8)


g_pam2 <- ggplot(pam2, aes(x = mice_log2fc,y = mcf7_cst_log2fc,label = gene_name))+
  geom_point(alpha = 0.7,size = 3,colour = "#814374")+
  xlab(NULL)+
  ylab(NULL)+
  geom_text_repel(max.overlaps = 500,size = 7)+
  scale_x_continuous(breaks = c(-2,0,2,4), limits = c(round(min(all$mice_log2fc), digits = 2) - 0.01, round(max(all$mice_log2fc), digits = 2)+0.01))+
  scale_y_continuous(breaks = c(0,4,8,12), limits = c(round(min(all$mcf7_cst_log2fc), digits = 2) - 0.01, round(max(all$mcf7_cst_log2fc), digits = 2)+0.01))+
  theme_classic()+
  theme(axis.text.x=element_text(size=20, face = "bold", colour = "black"),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1,face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

plot(g_pam2)

cor.test(pam2$mice_log2fc, pam2$mcf7_cst_log2fc, method="spearman", exact=FALSE)
ggsave(file = "cluster2_Mice.png", plot = g_pam2, dpi = 600, width = 7, height = 4.8)


g_pam3 <- ggplot(pam3, aes(x = mice_log2fc,y = mcf7_cst_log2fc,label = gene_name))+
  geom_point(alpha = 0.7,size = 3,colour = "#B7695C")+
  xlab(NULL)+
  ylab(NULL)+
  geom_text_repel(max.overlaps = 500,size = 7)+
  scale_x_continuous(breaks = c(-2,0,2,4), limits = c(round(min(all$mice_log2fc), digits = 2) - 0.01, round(max(all$mice_log2fc), digits = 2)+0.01))+
  scale_y_continuous(breaks = c(0,4,8,12), limits = c(round(min(all$mcf7_cst_log2fc), digits = 2) - 0.01, round(max(all$mcf7_cst_log2fc), digits = 2)+0.01))+
  theme_classic()+
  theme(axis.text.x=element_text(size=20, face = "bold", colour = "black"),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1,face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

plot(g_pam3)

cor.test(pam3$mice_log2fc, pam3$mcf7_cst_log2fc, method="spearman", exact=FALSE)
ggsave(file = "cluster3_Mice.png", plot = g_pam3, dpi = 600, width = 7, height = 4.8)


g_pam4 <- ggplot(pam4, aes(x = mice_log2fc,y = mcf7_cst_log2fc,label = gene_name))+
  geom_point(alpha = 0.7,size = 3,colour = "#9c7f17")+
  xlab(NULL)+
  ylab(NULL)+
  geom_text_repel(max.overlaps = 500,size = 7)+
  scale_x_continuous(breaks = c(-2,0,2,4), limits = c(round(min(all$mice_log2fc), digits = 2) - 0.01, round(max(all$mice_log2fc), digits = 2)+0.01))+
  scale_y_continuous(breaks = c(0,4,8,12), limits = c(round(min(all$mcf7_cst_log2fc), digits = 2) - 0.01, round(max(all$mcf7_cst_log2fc), digits = 2)+0.01))+
  theme_classic()+
  theme(axis.text.x=element_text(size=20, face = "bold", colour = "black"),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1,face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

plot(g_pam4)

cor.test(pam4$mice_log2fc, pam4$mcf7_cst_log2fc, method="spearman", exact=FALSE)
ggsave(file = "cluster4_Mice.png", plot = g_pam4, dpi = 600, width = 7, height = 4.8)