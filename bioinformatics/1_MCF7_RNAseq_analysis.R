library(biomaRt)
library(rtracklayer)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(clusterProfiler)

setwd("NFkB_aging_paper/")

StripVer <- function(ID){
  N <- regexpr("\\.", ID)[[1]]
  return(substr(ID, 1, N-1))
}

#Correspondence table between Ensembl geneID and Uniprot symbol----
gtf <- readGFF("ref_file/gencode.v40.annotation.gtf")
gtf_gene <- subset(gtf, gtf$type == "gene" & gtf$gene_type == "protein_coding")
gtf_gene <- gtf_gene[, c("gene_id", "gene_name")]
gtf_gene <- gtf_gene %>% distinct(gene_name, .keep_all=TRUE)
colnames(gtf_gene) <- c("Geneid", "gene_name")
rm(gtf)

#DEG detection ------
mcf7 <- read.table("source_file/salmon_MCF7_48h_gene_count.tsv", 
                   header=T, stringsAsFactors=F, sep="\t")
mcf7 <- dplyr::select(mcf7, "gene_name", "MCF7_CONTsi_TNF_1", "MCF7_CONTsi_TNF_2", "MCF7_CONTsi_TNF_3",
                      "MCF7_IKBAsi_TNF_1", "MCF7_IKBAsi_TNF_2", "MCF7_IKBAsi_TNF_3")

header <- c("gene_name","ctrl_rep1","ctrl_rep2","ctrl_rep3",
            "siIkBa_rep1","siIkBa_rep2","siIkBa_rep3")

colnames(mcf7) <- header
mcf7 <- dplyr::inner_join(mcf7,gtf_gene,by = "gene_name") 
mcf7 <- mcf7[!duplicated(mcf7$gene_name),] %>% dplyr::select(-Geneid)

mcf7 <- dplyr::filter(mcf7, (ctrl_rep1 >= 10 & ctrl_rep2 >= 10 & ctrl_rep3 >= 10) 
                      |(siIkBa_rep1 >= 10 & siIkBa_rep2 >= 10 & siIkBa_rep3 >= 10))

row.names(mcf7) <- mcf7$gene_name

mcf7 <- dplyr::select(mcf7, -"gene_name")

group <- data.frame(con = factor(c(rep("siCtrl", 3), rep("siIkBa", 3))))

dds <- DESeqDataSetFromMatrix(countData = round(mcf7), colData = group, design = ~ con)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds)

res <- results(dds, contrast = c("con", "siIkBa", "siCtrl"))
res$gene_id <- row.names(res)
res <- res[!is.na(res$padj),]
mcf7 <- as.data.frame(res)

fc_c <- 0.4 #foldchange criteria
padj_c <- 0.01 #pvalue criteria

DEG <- dplyr::filter(mcf7, mcf7$log2FoldChange > fc_c, mcf7$padj < padj_c)
DEG_down <- dplyr::filter(mcf7,mcf7$log2FoldChange < -fc_c, mcf7$padj < padj_c)
write.csv(DEG$human_gene_name,"DEG_up.csv",row.names = F)
write.csv(DEG_down$human_gene_name,"DEG_down.csv",row.names = F)

DEG <- dplyr::select(DEG, "gene_id") %>% dplyr::rename("human_gene_name" = 1)
DEG_down <- dplyr::select(DEG_down, "gene_id") %>% dplyr::rename("human_gene_name" = 1)

#volcanoplot (fig.S14A)
a <- mcf7
a_list <- a$gene_id

sp <- read.csv("ref_file/senescence_related_gene.csv") #read senescense gene list 

ttmp <- rep(1, nrow(a))
ttmp[a$log2FoldChange > fc_c & a$padj < padj_c] <- 2
ttmp[a$log2FoldChange < -fc_c & a$padj < padj_c] <- 3
a$threshold <- as.factor(ttmp)

a_list[!a_list %in% sp$x] <- ""
a$gene_id <- a_list

for(i in 1:nrow(a)){
  if (a[i,8] == 1){
    a[i,7] <- ""
  }
}

cols <- c("#808080", "#BC0000","#0000BC")
g <- ggplot(data=a, aes(x=log2FoldChange, y=-log10(padj), colour=threshold, label = gene_id)) +
  geom_point(size=1.75) +
  geom_text_repel(max.overlaps = 1000000, size = 4)+
  theme_classic()+
  theme(legend.position = "none") +
  xlab(NULL) +
  ylab(NULL) +
  geom_hline(yintercept = -log10(padj_c),linetype="dashed",lwd = 1)+
  geom_vline(xintercept = c(fc_c,-fc_c),linetype="dashed",lwd = 1)+
  scale_fill_manual(values=cols) + scale_color_manual(values=cols)+
  theme(axis.text.x=element_text(size=20, face = "bold", colour = "black"),
        axis.text.y=element_text(size=20, face = "bold", colour = "black"),
        axis.title=element_text(size=1,face="bold", colour = "black"),
        axis.line = element_line(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
plot(g)
ggsave(file = "volcanoplot.png", plot = g, dpi = 600, width = 10, height = 10)


# clustering by expression pattern -----
mcf7_tc <- read.table("source_file/salmon_MCF7_timecourse_gene_TPM.tsv", 
                      header=T, stringsAsFactors=F,sep="\t")
mcf7_tc <- dplyr::select(mcf7_tc, "gene_name", 
                         "siIkBa_0_min_rep1", "siIkBa_0_min_rep2", "siIkBa_15_min_rep1", "siIkBa_15_min_rep2", 
                         "siIkBa_30_min_rep1", "siIkBa_30_min_rep2", "siIkBa_45_min_rep1", "siIkBa_45_min_rep2", 
                         "siIkBa_60_min_rep1", "siIkBa_60_min_rep2", "siIkBa_75_min_rep1", "siIkBa_75_min_rep2", 
                         "siIkBa_90_min_rep1", "siIkBa_90_min_rep2", "siIkBa_105_min_rep1", "siIkBa_105_min_rep2", 
                         "siIkBa_120_min_rep1", "siIkBa_120_min_rep2", "siIkBa_135_min_rep1", "siIkBa_135_min_rep2",
                         "siIkBa_150_min_rep1", "siIkBa_150_min_rep2")
mcf7_tc <- mcf7_tc[!duplicated(mcf7_tc$gene_name), ]
mcf7_tc_genelist <- dplyr::select(mcf7_tc,"gene_name")

mcf7_tc <- mcf7_tc[,-1]

min_value <- apply(mcf7_tc,1,function(x) if(all(x==0)) 0 else min(x[x>0])) #Calculation of the minimum value excluding 0
min_value <- min_value[min_value > 0] %>% min()

mcf7_tc <- log2(mcf7_tc + min_value/2) #Add half of the minimum to the whole, then logarithmize by 2
mcf7_tc <- cbind(mcf7_tc_genelist,mcf7_tc)

mcf7_tc$siIkBa_0 <- (mcf7_tc$siIkBa_0_min_rep1 + mcf7_tc$siIkBa_0_min_rep2)/2
mcf7_tc$siIkBa_15 <- (mcf7_tc$siIkBa_15_min_rep1 + mcf7_tc$siIkBa_15_min_rep2)/2
mcf7_tc$siIkBa_30 <- (mcf7_tc$siIkBa_30_min_rep1 + mcf7_tc$siIkBa_30_min_rep2)/2
mcf7_tc$siIkBa_45 <- (mcf7_tc$siIkBa_45_min_rep1 + mcf7_tc$siIkBa_45_min_rep2)/2
mcf7_tc$siIkBa_60 <- (mcf7_tc$siIkBa_60_min_rep1 + mcf7_tc$siIkBa_60_min_rep2)/2
mcf7_tc$siIkBa_75 <- (mcf7_tc$siIkBa_75_min_rep1 + mcf7_tc$siIkBa_75_min_rep2)/2
mcf7_tc$siIkBa_90 <- (mcf7_tc$siIkBa_90_min_rep1 + mcf7_tc$siIkBa_90_min_rep2)/2
mcf7_tc$siIkBa_105 <- (mcf7_tc$siIkBa_105_min_rep1 + mcf7_tc$siIkBa_105_min_rep2)/2
mcf7_tc$siIkBa_120 <- (mcf7_tc$siIkBa_120_min_rep1 + mcf7_tc$siIkBa_120_min_rep2)/2
mcf7_tc$siIkBa_135 <- (mcf7_tc$siIkBa_135_min_rep1 + mcf7_tc$siIkBa_135_min_rep2)/2
mcf7_tc$siIkBa_150 <- (mcf7_tc$siIkBa_150_min_rep1 + mcf7_tc$siIkBa_150_min_rep2)/2

mcf7_tc <- mcf7_tc[c(1,24:ncol(mcf7_tc))] %>% dplyr::rename("human_gene_name" = 1)

mcf7_tc <- dplyr::inner_join(DEG,mcf7_tc, by = "human_gene_name")
row.names(mcf7_tc) <- mcf7_tc$human_gene_name
mcf7_tc <- dplyr::select(mcf7_tc, -"human_gene_name")

z_scored_DEG <- apply(mcf7_tc, 1, scale) %>% t() %>% as.data.frame()
z_scored_DEG <- as.data.frame(z_scored_DEG)
colnames(z_scored_DEG) <- c("0","15","30","45","60","75","90","105","120","135","150")
z_scored_DEG[is.na(z_scored_DEG)] <- 0

set.seed(921)
pa <- cluster::pam(z_scored_DEG, k = 4)

#Heatmap (fig.3A)
cols <- c("#483D8B", "#EEEEEE", "#CD2626")
scale <- c(-2, 0 ,2)
col_fun = colorRamp2(scale, cols)
lgd = Legend(col_fun = col_fun, title = "z-score", at = scale, 
             labels = c("", "", ""), direction = "horizontal", legend_width = unit(2, "cm"),)
pdf("mcf7_timecourse.pdf", width = 3.5)
Heatmap(z_scored_DEG, name = "z_score", cluster_columns = F, cluster_rows = F, height = unit(13, "cm"),
        col = col_fun,
        split = paste0("pam", pa$clustering),
        show_row_names = F,
        show_heatmap_legend = F,
        border_gp = gpar(col = "#1a1a1a", lty = 1, lwd = 2), )
draw(lgd, x = unit(0.87, "npc"), y = unit(0.95, "npc"))
dev.off()

#Obtain gene list for each cluster
extract_cluster_genes <- function(x){
  d <- names(which(pa$clustering == x)) %>% Vectorize(StripVer)()
  d <- data.frame(gene_id = d)
  d$gene_name <- rownames(d)
  d <- as.data.frame(d$gene_name)
  d <- as.data.frame(d[order(d$`d$gene_name`, decreasing = F), ])
  colnames(d) <- "gene_name"
  return(d)
}

pam1 <- extract_cluster_genes(1)
pam2 <- extract_cluster_genes(2)
pam3 <- extract_cluster_genes(3)
pam4 <- extract_cluster_genes(4)

write.csv(pam1, "pam1.csv", row.names = F)
write.csv(pam2, "pam2.csv", row.names = F)
write.csv(pam3, "pam3.csv", row.names = F)
write.csv(pam4, "pam4.csv", row.names = F)

#Calculate the average z-score
z_scored_DEG$gene_name <- rownames(z_scored_DEG)

cal_a_z <- function(pam){
  z_scored_pam <- dplyr::inner_join(pam, z_scored_DEG, by = "gene_name")
  rownames(z_scored_pam) <- z_scored_pam$gene_name
  z_scored_pam <- dplyr::select(z_scored_pam, -"gene_name")
  z_scored_pam_mean <- as.data.frame(colMeans(z_scored_pam))
  z_scored_pam_mean <- t(z_scored_pam_mean)
  return(z_scored_pam_mean)
}

z_scored_DEG_pam1_mean <- cal_a_z(pam1)
z_scored_DEG_pam2_mean <- cal_a_z(pam2)
z_scored_DEG_pam3_mean <- cal_a_z(pam3)
z_scored_DEG_pam4_mean <- cal_a_z(pam4)

z_mean <- rbind(z_scored_DEG_pam1_mean, z_scored_DEG_pam2_mean, 
                z_scored_DEG_pam3_mean, z_scored_DEG_pam4_mean)
z_mean <- as.data.frame(z_mean)
rownames(z_mean) <- c("cluster1","cluster2","cluster3","cluster4")
z_mean <- as.matrix(z_mean) %>% melt()
colnames(z_mean) <- c("cluster","time","z_score")

#Lineplot (fig.3A)
g <- ggplot(z_mean, aes(x = time, y = z_score, color = cluster))+
  geom_point( size = 3)+
  xlab("")+
  ylab("")+
  scale_color_manual(values = c("black","black","black","black"))+
  scale_y_continuous(breaks = c(-1,0,1),limits = c(-1.5,1.5) )+
  stat_smooth(method = 'loess',formula = y ~ x,se = TRUE, size = 2)+
  theme_classic()+
  facet_wrap(. ~ cluster, nrow = 4)+
  theme(strip.text.x = element_text(size=15, face="bold"),
        axis.line = element_line(size = 1.2, colour = "black"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")
plot(g)
ggsave(file = "z_mean.png", plot = g, dpi = 600, width = 3, height = 10)

#KEGG enrichment analysis----
R.utils::setOption("clusterProfiler.download.method", "auto")
colnames(gtf_gene) <- c("gene_id", "gene_name")
gtf_gene$gene_id <- gsub("\\..*", "", gtf_gene$gene_id)

pre_enrich <- function(x){
  x <- dplyr::inner_join(x, gtf_gene, by = "gene_name") 
  x <- bitr(x$gene_id, fromType="ENSEMBL", toType="ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
  return(x)
}

cluster1 <- pre_enrich(pam1)
cluster2 <- pre_enrich(pam2)
cluster3 <- pre_enrich(pam3)
cluster4 <- pre_enrich(pam4)

cls <- list(cluster1 = cluster1, cluster2 = cluster2, cluster3 = cluster3, cluster4 = cluster4)

comp_ego <- compareCluster(cls, fun = "enrichKEGG",
                           organism = "hsa",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1)

comp_ego <- as.data.frame(comp_ego)

cluster1_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "cluster1" & comp_ego$p.adjust < 0.05 & comp_ego$Count >= 5)
cluster2_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "cluster2" & comp_ego$p.adjust < 0.05 & comp_ego$Count >= 5)
cluster3_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "cluster3" & comp_ego$p.adjust < 0.05 & comp_ego$Count >= 5)
cluster4_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "cluster4" & comp_ego$p.adjust < 0.05 & comp_ego$Count >= 5)

term_list <- c(cluster1_comp$Description, cluster2_comp$Description, cluster3_comp$Description, cluster4_comp$Description) %>% 
  unique() %>% as.data.frame() %>% dplyr::rename("Description" = 1)

comp_ego <- dplyr::inner_join(comp_ego, term_list, by = "Description")

cluster1_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "cluster1") %>% dplyr::select("Description","pvalue") %>% dplyr::rename("cluster1" = 2)
cluster1_comp$cluster1 <- -log10(cluster1_comp$cluster1)
cluster2_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "cluster2") %>% dplyr::select("Description","pvalue") %>% dplyr::rename("cluster2" = 2)
cluster2_comp$cluster2 <- -log10(cluster2_comp$cluster2)
cluster3_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "cluster3") %>% dplyr::select("Description","pvalue") %>% dplyr::rename("cluster3" = 2)
cluster3_comp$cluster3 <- -log10(cluster3_comp$cluster3)
cluster4_comp <- dplyr::filter(comp_ego,comp_ego$Cluster == "cluster4") %>% dplyr::select("Description","pvalue") %>% dplyr::rename("cluster4" = 2)
cluster4_comp$cluster4 <- -log10(cluster4_comp$cluster4)

comp_ego <- dplyr::full_join(cluster1_comp,cluster2_comp, by = "Description") %>%
  dplyr::full_join(., cluster3_comp, by = "Description") %>% 
  dplyr::full_join(., cluster4_comp, by = "Description")

row.names(comp_ego) <- comp_ego$Description
comp_ego <- dplyr::select(comp_ego, -"Description")
comp_ego[is.na(comp_ego)] <- 0

#KEGG Heatmap (fig.S16C)
cols <- c("#fef9fb","#C96363","#b22222")
scale <- c(0, 2 ,4)
col_fun = colorRamp2(scale, cols)
lgd = Legend(col_fun = col_fun, title = expression(paste("-", {Log[10]}, "(pvalue)", sep="")), at = scale, 
             labels = scale, direction = "horizontal",legend_width = unit(2, "cm"),)

pdf("KEGG_heatmap.pdf",width = 6)
Heatmap(comp_ego, name = "-log10(p)", cluster_columns = F,cluster_rows = T, height = unit(6, "cm"),
        rect_gp = gpar(col = "#1a1a1a", lty = 1, lwd = 2),
        row_names_max_width = max_text_width(rownames(comp_ego)),
        col = col_fun,
        clustering_method_rows = "ward.D2",
        cluster_row_slices = FALSE,
        row_title = NULL,
        show_heatmap_legend = F)
draw(lgd,x = unit(0.15, "npc"), y = unit(0.77, "npc"))
dev.off()
