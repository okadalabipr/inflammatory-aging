# Bioinformatic analysis of NF-κB dynamics on inflammatory aging
## Introduction
This repository contains the source code for the sequence analysis.  
Primary analysis was performed using an established pipeline.  
* RNA-seq  -> https://nf-co.re/rnaseq/3.5  
* ChIP-seq -> https://nf-co.re/chipseq/1.2.2  
* ATAC-seq -> https://github.com/ENCODE-DCC/atac-seq-pipeline

## Requirements 
### for R
R version 4.1.3
| Package Name                                                   | Version                                                                                                  | 
| -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- | 
| biomaRt	                                                     | 2.50.3                                                                                                   |
| circlize	                                                     | 0.4.15                                                                                                   |
| clusterProfiler                                                | 4.2.2                                                                                                    |
| ComplexHeatmap	                                             | 2.10.0                                                                                                   |
| cummeRbund	                                                 | 2.36.0                                                                                                   |
| DESeq2	                                                     | 1.34.0                                                                                                   |
| dplyr	                                                         | 1.0.10                                                                                                   |
| ggplot2	                                                     | 3.3.6                                                                                                    |
| ggrepel	                                                     | 0.9.1                                                                                                    |
| NSM3	                                                         | 1.17                                                                                                     |
| reshape2	                                                     | 1.4.4                                                                                                    |
| rtracklayer	                                                 | 1.54.0                                                                                                   |
| VennDiagram	                                                 | 1.7.3                                                                                                    |

### for shell script
| Package Name                                                   | Version                                                                                                  | 
| -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- | 
| deeptools	                                                     | 3.5.1                                                                                                    |
| Homer 	                                                     | 4.11                                                                                                     |

## Brief usage instruction
Please add reference files in [`ref_file/`](./ref_file/). Required files are listed in each script.  
To reproduce results, please perform the analysis in the following order.  
* 1_MCF7_RNAseq_analysis.R
* 2_MCF7_ATACseq_analysis.R
* 3_MCF7_ChIPseq_analysis.sh
* 4_MCF7_RELA_peaknum_analysis.R
* 5_RNAseq_Correlation_anarysis.R
* 6_Mouse_Rela_peaknum_analysis.R

## License