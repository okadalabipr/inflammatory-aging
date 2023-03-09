#!/bin/bash

RELA_ChIP_siCtrl_120 = "path to your RELA ChIP-seq bigWig file in Ctrl condition obtained from nfcore ChIP-seq pipeline (https://nf-co.re/chipseq/1.2.2)"
RELA_ChIP_siIkBa_120 = "path to your RELA ChIP-seq bigWig file in Ctrl condition obtained from nfcore ChIP-seq pipeline (https://nf-co.re/chipseq/1.2.2)"
BED = "path to your bed file of open chromatin regions near each cluster obtained from 2_MCF7_ATACseq_analysis.R"

#fig.3C
computeMatrix reference-point \
    -S ${RELA_ChIP_siCtrl_120} ${RELA_ChIP_siIkBa_120} \
    -R ${BED} \
    -a 1500 -b 1500 --referencePoint center -o RelA_120_matrix.mat.gz -p 8

plotHeatmap -m RelA_120_matrix.mat.gz --heatmapHeight 15 --whatToShow "heatmap and colorbar" --colorMap coolwarm --colorList darkslateblue,lightgray,firebrick --heatmapWidth 7 -out RelA_120_Heatmap.png
plotProfile -m RelA_120_matrix.mat.gz --plotType se --plotWidth 7 --colors "#51A39D" "#814374" "#B7695C" "#805f2b" -out RelA_120_plot.png


