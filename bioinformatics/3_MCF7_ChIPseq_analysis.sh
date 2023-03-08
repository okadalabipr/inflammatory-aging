#!/bin/bash

#fig.3C
computeMatrix reference-point -S ${R_home_directory}/NFkB_aging_paper/source_file/RELA_ChIP_siCtrl_120.bigWig ${R_home_directory}/NFkB_aging_paper/source_file/RELA_ChIP_siIkBa_120.bigWig -R ${R_home_directory}/NFkB_aging_paper/ATACseq_analysis/*_near.bed -a 1500 -b 1500 --referencePoint center -o ${R_home_directory}/${R_home_directory}/NFkB_aging_paper/RelA_120_matrix.mat.gz -p 8
plotHeatmap -m ${R_home_directory}/NFkB_aging_paper/RelA_120_matrix.mat.gz --heatmapHeight 15 --whatToShow "heatmap and colorbar" --colorMap coolwarm --colorList darkslateblue,lightgray,firebrick --heatmapWidth 7 -out ${R_home_directory}/NFkB_aging_paper/RelA_120_Heatmap.png
plotProfile -m ${R_home_directory}/NFkB_aging_paper/RelA_120_matrix.mat.gz --plotType se --plotWidth 7 --colors "#51A39D" "#814374" "#B7695C" "#805f2b" -out ${R_home_directory}/NFkB_aging_paper/RelA_120_plot.png



