library(cummeRbund)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(VennDiagram)
library(NSM3)

setwd("inframmatory-aging/bioinfomatics/")
dir.create("mouse_analysis/")

#List of genes registered in gencode (mouse)
gtf <- readGFF("path to your gencode annotation file obtained from https://www.gencodegenes.org/mouse/release_M31.html")
gtf_gene <- subset(gtf, gtf$type == "gene" & gtf$gene_type == "protein_coding")
gtf_gene <- gtf_gene$gene_name %>% as.data.frame() %>% dplyr::rename("gene_name" = 1) %>% distinct(gene_name, .keep_all=TRUE)
rm(gtf)

#load Rela peak file 
rela_annotate <- read.csv("path to your Mouse Rela ChIP-peak file obtained (Details are described in Method)",sep = "\t")
rela_target <- dplyr::filter(rela_annotate,abs(rela_annotate$Distance.to.TSS) < 500) %>% dplyr::rename("gene_name" = 16)
rela_target <- dplyr::inner_join(rela_target,gtf_gene, by = "gene_name") 
rela_target <- rela_target$gene_name %>% unique() %>% as.data.frame()
colnames(rela_target) <- "gene_id"

#Classification by rate of increase in gene expression with aging (mouse heart)-----
cnt <- read.table("path to your Mouse Heart RNA-seq count data obtained from nf-core pipeline (https://nf-co.re/rnaseq/3.5)", header=T, stringsAsFactors=F,sep="\t")
cnt <- dplyr::select(cnt, -"gene_id")

header <- c("gene_name", "old_rep1", "old_rep2", "old_rep3", "old_rep4", "old_rep5",
            "young_rep1", "young_rep2", "young_rep3", "young_rep4", "young_rep5")
colnames(cnt) <- header

cnt <- dplyr::inner_join(gtf_gene, cnt, by = "gene_name")
cnt <- cnt[!duplicated(cnt$gene_name), ]

cnt <- dplyr::filter(cnt, (old_rep1 >= 5 & old_rep2 >= 5 & old_rep3 >= 5 & old_rep4 >= 5 & old_rep5 >= 5) 
                     | (young_rep1 >= 5 & young_rep2 >= 5 & young_rep3 >= 5 & young_rep4 >= 5 & young_rep5 >= 5))

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
heart <- as.data.frame(res)

rela_target_heart <- dplyr::inner_join(heart, rela_target, by = "gene_id")

DEG_all <- dplyr::filter(rela_target_heart, rela_target_heart$padj < 0.01, rela_target_heart$log2FoldChange > 0)
DEG0.5 <- dplyr::filter(DEG_all, DEG_all$log2FoldChange <= 0.5)
DEG1.0 <- dplyr::filter(DEG_all, 0.5 < DEG_all$log2FoldChange, DEG_all$log2FoldChange <= 1.0)
DEG1.0_over <- dplyr::filter(DEG_all, DEG_all$log2FoldChange > 1.0)
rela_potential_nochange <- setdiff(rela_target_heart$gene_id, DEG_all$gene_id) %>% as.data.frame() %>% dplyr::rename("gene_id" = 1)

#venn plot (fig.S17B)
gene_list <- list(rela_target_heart$gene_id, DEG_all$gene_id)
venn.diagram(gene_list, filename = "mouse_analysis/venn.png", imagetype="png", ext.text = F, lty = rep("blank", 2),
             cex = 1, fontface = "bold", cat.pos = c(180, 180), cat.fontface = "bold", cat.cex = 2,
             category.names = c("", ""), fill = c("#add8e6", "#b34b40"))

DEG_all <- DEG_all[order(DEG_all$log2FoldChange, decreasing = T), ]
DEG_all$rank <- seq(1, nrow(DEG_all))

ttmp <- rep(1, nrow(DEG_all))
ttmp[DEG_all$log2FoldChange > 0 & DEG_all$log2FoldChange < 0.5] <- 2
ttmp[DEG_all$log2FoldChange > 0.5 & DEG_all$log2FoldChange < 1.0] <- 3
ttmp[DEG_all$log2FoldChange > 1.0] <- 4
DEG_all$threshold <- as.factor(ttmp)

#bar plot (fig.S17B)
g <- ggplot(DEG_all, aes(x = rank, y = log2FoldChange,fill = threshold))+
  geom_bar(stat = "identity")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  scale_x_reverse(limits = c(nrow(DEG_all)+1, 0), breaks = c(nrow(DEG1.0_over) + nrow(DEG1.0), nrow(DEG1.0_over)), expand = c(0, 0))+
  scale_y_continuous(breaks = c(0, 0.5, 1.0, 2.0, 4.0, 6.0), expand = c(0, 0))+
  scale_fill_manual(values = c("#f15c5c", "#b84a39", "#6f3826"))+
  geom_segment(aes(x = nrow(DEG1.0_over) + nrow(DEG1.0), y = 0.5, xend = nrow(DEG_all)+1, yend = 0.5, colour = "segment"),linetype="dashed", lwd = 1.2,colour = "black")+
  geom_segment(aes(x = nrow(DEG1.0_over), y = 1, xend = nrow(DEG_all)+1, yend = 1, colour = "segment"),linetype="dashed",lwd = 1.2,colour = "black")+
  theme(axis.text.x=element_text(size=20, face = "bold", colour = "black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1, face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
plot(g)
ggsave(file = "mouse_analysis/DEG_rank.png", plot = g, dpi = 600, width = 7, height = 5)

#count of RELA peaks----
rela_annotate <- rela_annotate[, c(2,3,4,10,16)] %>% dplyr::rename("gene_id" = 5) %>% dplyr::filter(abs(rela_annotate$Distance.to.TSS) < 10000)

peak_count <- function(x){
  x_list <- x$gene_id %>% unique()
  x_score <- NULL
  for (i in 1:length(x_list)){
    name <- x_list[i]
    tmp <- dplyr::filter(rela_annotate, rela_annotate$gene_id == name)
    score <- nrow(tmp)
    tmp <- c(name, score) %>% t() %>% as.data.frame()
    x_score <- rbind(x_score, tmp)
  }
  colnames(x_score) <-  c("gene_id", "num")
  x_score$num <- as.numeric(x_score$num)
  return(x_score)
}

rela_potential_nochange_score <- peak_count(rela_potential_nochange)
DEG_0.5_score <- peak_count(DEG0.5)
DEG_1.0_score <- peak_count(DEG1.0)
DEG_1.0_over_score <- peak_count(DEG1.0_over)

df <- data.frame(
  label = c(rep("nochange", length(rela_potential_nochange_score$num)), rep("~0.5",length(DEG_0.5_score$num)),
            rep("0.5~1.0", length(DEG_1.0_score$num)), rep("1.0~", length(DEG_1.0_over_score$num))),
  num = c(rela_potential_nochange_score$num, DEG_0.5_score$num, DEG_1.0_score$num, DEG_1.0_over_score$num))

df2 <- group_by(df, label) %>% 
  summarise_all(list(mean = mean, 
                     sd = sd, 
                     se = ~sd/sqrt(length(.))))

df2$factor <- factor(df2$label, levels = c("nochange","~0.5","0.5~1.0","1.0~"))

#Bar plot of average Rela peaks (fig.3H)
g <- ggplot(df2, aes(x=factor, y = mean, fill=label))+
  geom_bar(stat = "identity",colour = "black",lwd = 0.7,width = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),width = .2)+
  theme_classic(base_size = 20)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#f15c5c","#b84a39","#6f3826","#d6ecfa"))+
  scale_y_continuous(breaks = seq(floor(max(df2$mean)), 0, by = -1), limits = c(0, round(max(df2$mean), digits = 1)+round(max(df2$se), digits = 1) ), expand = c(0, 0) )+
  labs(x="", y="")+
  theme(axis.text.x=element_text(size=20, face = "bold", colour = "black",angle = 45,hjust = 1),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1,face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
plot(g)
ggsave(file = "mouse_analysis/num_rela.png", plot = g, dpi = 600, width = 3.6, height = 4.8)

#Steel-Dwass test 
df$num <- as.factor(df$num)
df$label <- as.factor(df$label)
test <- pSDCFlig(df$num, df$label, method="Asymptotic")
test <- data.frame(label = test$labels, p.val = test$p.val)

#KEGG enrichment analysis-----
#Correspondence table between Ensembl geneID and Uniprot symbol
gtf <- readGFF("path to your gencode annotation file obtained from https://www.gencodegenes.org/mouse/release_M31.html")
gtf_gene <- subset(gtf, gtf$type == "gene" & gtf$gene_type == "protein_coding")
gtf_gene <- gtf_gene[, c("gene_id","gene_name")]
gtf_gene$gene_id <- gsub("\\..*", "", gtf_gene$gene_id)
gtf_gene <- gtf_gene %>% distinct(gene_name, .keep_all=TRUE)
colnames(gtf_gene) <- c("gene_id2", "gene_id")
rm(gtf)

pre_enrich <- function(x){
  x <- dplyr::inner_join(x, gtf_gene, by = "gene_id") 
  x <- bitr(x$gene_id2, fromType="ENSEMBL", toType="ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID
  return(x)
}

R.utils::setOption("clusterProfiler.download.method", "auto")

DEG0.5 <- pre_enrich(DEG0.5)
DEG1.0 <- pre_enrich(DEG1.0)
DEG1.0_over <- pre_enrich(DEG1.0_over)

cls <- list("~0.5" = DEG0.5, "0.5~1.0" = DEG1.0, "1.0~" = DEG1.0_over)

comp_ego <- compareCluster(cls, fun = "enrichKEGG",
                           organism = "mmu",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1)

comp_ego <- as.data.frame(comp_ego)

DEG0.5_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "~0.5" & comp_ego$p.adjust < 0.05) %>% head(3)
DEG1.0_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "0.5~1.0" & comp_ego$p.adjust < 0.05) %>% head(3)
DEG1.0_over_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "1.0~" & comp_ego$p.adjust < 0.05) %>% head(3)

term_list <- c(DEG0.5_comp$Description,DEG1.0_comp$Description,DEG1.0_over_comp$Description) %>% unique() %>% as.data.frame() %>% dplyr::rename("Description" = 1) 
comp_ego <- dplyr::inner_join(comp_ego,term_list,by = "Description")
DEG0.5_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "~0.5") %>% dplyr::select(c(3,6)) %>% dplyr::rename("Low" = 2)
DEG0.5_comp$Low <- -log10(DEG0.5_comp$Low)
DEG1.0_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "0.5~1.0") %>% dplyr::select(c(3,6)) %>% dplyr::rename("Medium" = 2)
DEG1.0_comp$Medium <- -log10(DEG1.0_comp$Medium)
DEG1.0_over_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "1.0~") %>% dplyr::select(c(3,6)) %>% dplyr::rename("High" = 2)
DEG1.0_over_comp$High <- -log10(DEG1.0_over_comp$High)

comp_ego <- dplyr::full_join(DEG0.5_comp,DEG1.0_comp, by = "Description")
comp_ego <- dplyr::full_join(comp_ego,DEG1.0_over_comp, by = "Description")
row.names(comp_ego) <- comp_ego$Description
comp_ego <- comp_ego[,-1]
comp_ego[is.na(comp_ego)] <- 0

#KEGG heatmap (fig.3I)
cols <- c("#fef9fb","#C96363","#b22222")
scale <- c(0, 5 ,10)
col_fun = colorRamp2(scale, cols)
lgd = Legend(col_fun = col_fun, title = expression(paste("-", {Log[10]}, "(pvalue)", sep="")), at = scale, 
             labels = c("","",""), direction = "horizontal",legend_width = unit(2, "cm"),)
set.seed(921)
pdf("mouse_analysis/KEGG_heatmap.pdf",width = 6)
Heatmap(comp_ego, name = "-log10(p)", cluster_columns = F,cluster_rows = T,height = unit(7, "cm"),
        rect_gp = gpar(col = "#1a1a1a", lty = 1, lwd = 2),
        row_names_max_width = max_text_width(rownames(comp_ego)),
        col = col_fun,
        clustering_method_rows = "ward.D2",
        cluster_row_slices = FALSE,
        row_title = NULL,
        show_heatmap_legend = F)
draw(lgd,x = unit(0.2, "npc"), y = unit(0.79, "npc"))
dev.off()
