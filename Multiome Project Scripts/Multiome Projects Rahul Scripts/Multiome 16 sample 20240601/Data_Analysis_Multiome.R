#------------------------------------------------------------------------------
# Original Author and Code Provided By: Kory R Johnson
# Affiliation: NINDS/NIH Bioinformatics Core
# Contact information: johnsonko@mail.nih.gov
# Date: 20240608

# Edited for use by Rahul Akkem
# Sara Inati Lab 
# Editing: 20230928

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Launch
#------------------------------------------------------------------------------
#cd /data/johnsonko/Analysis_Projects/Lorna_Role/Rajebhosale_Prithviraj/scRNA_Seq_Project_20230522/h5_files
#cd /Users/akkemrr/Desktop/Genomic_scRNA_seq_Analysis


#sinteractive --mem=200g --cpus-per-task=4 --time 36:00:00 --gres=lscratch:200 # augur
#sinteractive --gres=lscratch:5

#------------------------------------------------------------------------------
# Set-up
#------------------------------------------------------------------------------

#module load R Rstudio
#module load Rstudio
#rstudio &

#------------------------------------------------------------------------------
# Initialize
#------------------------------------------------------------------------------

#rm(list=ls())
#options(object.size=Inf)
# setwd("/data/johnsonko/Analysis_Projects/Lorna_Role/Rajebhosale_Prithviraj/scRNA_Seq_Project_20230522/h5_files")

#setwd("/Users/akkemrr/Desktop/Genomic_scRNA_seq_Analysis")


#------------------------------------------------------------------------------
# Install
#------------------------------------------------------------------------------

# install.packages("usethis")
# devtools::install_github("hadley/devtools")

library(devtools)

# install.packages("harmony")
library("Rcpp")
library(harmony)

# BiocManager::install("DropletUtils")
library(DropletUtils)

# BiocManager::install("scDblFinder") 
library(scDblFinder)

#remotes::install_version("Seurat", version = "5.0.3") 
#devtools::install_github('satijalab/seurat-data')
library(Seurat)
library(SeuratData)
library(cowplot)

# remotes::install_github('satijalab/azimuth', version = "0.5.0") # not working now
# devtools::install_github('satijalab/azimuth', version = "0.5.0")
# remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
library(Azimuth)

# BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

# BiocManager::install("ggplot2")
library(ggplot2)

# BiocManager::install("scales")
library(scales)

# BiocManager::install("dplyr") 
library("dplyr")

# BiocManager::install("clustree")
library("clustree")

# BiocManager::install("scater")
library("scater") # no work here and above

# BiocManager::install("Signac")
library(Signac)

# BiocManager::install("EnsDb.Mmusculus.v79")
library(EnsDb.Mmusculus.v79)

# BiocManager::install("GenomicRanges")
library(GenomicRanges)

# BiocManager::install("future")
library(future)

# BiocManager::install("stringr")
library(stringr)

# BiocManager::install("openxlsx")
library(openxlsx)

# BiocManager::install("pheatmap")
library(pheatmap)

'%!in%' <- Negate('%in%') # WHAT IS THIS?

# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)

# BiocManager::install("chromVAR")
library(chromVAR)

# BiocManager::install("JASPAR2020")
library(JASPAR2020)

# BiocManager::install("TFBSTools")
library(TFBSTools)

# BiocManager::install("motifmatchr")
library(motifmatchr)

# devtools::install_github("neurorestore/Augur")
library(Augur)
library(BiocParallel)

# BiocManager::install("openxlsx")
library(openxlsx)

# BiocManager::install("patchwork")
library(patchwork)

#------------------------------------------------------------------------------
# Import & Split
#------------------------------------------------------------------------------

s1356_A <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1356_A_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1356_A.gex<-s1356_A[[1]]
s1356_A.seurat <- CreateSeuratObject(counts = s1356_A.gex, project="s1356_A",min.cells = 3, min.features=200)

s1356_C <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1356_C_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1356_C.gex<-s1356_C[[1]]
s1356_C.seurat <- CreateSeuratObject(counts = s1356_C.gex,project="s1356_C",min.cells = 3, min.features=200)

s1521_4 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1521_4_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1521_4.gex<-s1521_4[[1]]
s1521_4.seurat <- CreateSeuratObject(counts = s1521_4.gex,project="s1521_4",min.cells = 3, min.features=200)

s1521_7 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1521_7_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1521_7.gex<-s1521_7[[1]]
s1521_7.seurat <- CreateSeuratObject(counts = s1521_7.gex,project="s1521_7",min.cells = 3, min.features=200)

s1295_1 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1295_1_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1295_1.gex<-s1295_1[[1]]
s1295_1.seurat <- CreateSeuratObject(counts = s1295_1.gex,project="s1295_1",min.cells = 3, min.features=200)

s1295_3 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1295_3_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1295_3.gex<-s1295_1[[1]]
s1295_3.seurat <- CreateSeuratObject(counts = s1295_3.gex,project="s1295_3",min.cells = 3, min.features=200)

s1717_2 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1717_2_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1717_2.gex<-s1717_2[[1]]
s1717_2.seurat <- CreateSeuratObject(counts = s1717_2.gex,project="s1717_2",min.cells = 3, min.features=200)

s1717_3 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1717_3_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1717_3.gex<-s1717_3[[1]]
s1717_3.seurat <- CreateSeuratObject(counts = s1717_3.gex,project="s1717_3",min.cells = 3, min.features=200)

s1815_2 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1815_2_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1815_2.gex<-s1815_2[[1]]
s1815_2.seurat <- CreateSeuratObject(counts = s1815_2.gex,project="s1815_2",min.cells = 3, min.features=200)

s1815_4 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1815_4_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1815_4.gex<-s1815_4[[1]]
s1815_4.seurat <- CreateSeuratObject(counts = s1815_4.gex,project="s1815_4",min.cells = 3, min.features=200)

s1822_1 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1822_1_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1822_1.gex<-s1822_1[[1]]
s1822_1.seurat <- CreateSeuratObject(counts = s1822_1.gex,project="s1822_1",min.cells = 3, min.features=200)

s1822_2 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/1822_2_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s1822_2.gex<-s1822_2[[1]]
s1822_2.seurat <- CreateSeuratObject(counts = s1822_2.gex,project="s1822_2",min.cells = 3, min.features=200)

s2281_1 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/2281_1_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s2281_1.gex<-s2281_1[[1]]
s2281_1.seurat <- CreateSeuratObject(counts = s2281_1.gex,project="s2281_1",min.cells = 3, min.features=200)

s2281_3 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/2281_3_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s2281_3.gex<-s2281_3[[1]]
s2281_3.seurat <- CreateSeuratObject(counts = s2281_3.gex,project="s2281_3",min.cells = 3, min.features=200)

s2392_1 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/2392_1_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s2392_1.gex<-s2392_1[[1]]
s2392_1.seurat <- CreateSeuratObject(counts = s2392_1.gex,project="s2392_1",min.cells = 3, min.features=200)

s2392_3 <- Read10X_h5("/data/akkemrr/Single_RNA_Analy/2392_3_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
s2392_3.gex<-s2392_3[[1]]
s2392_3.seurat <- CreateSeuratObject(counts = s2392_3.gex,project="s2392_3",min.cells = 3, min.features=200)


Idents(s1356_A.seurat) <- "s1356_A"
Idents(s1356_C.seurat) <- "s1356_C"
Idents(s1521_4.seurat) <- "s1521_4"
Idents(s1521_7.seurat) <- "s1521_7"
Idents(s1295_1.seurat) <- "s1295_1"
Idents(s1295_3.seurat) <- "s1295_3"
Idents(s1717_2.seurat) <- "s1717_2"
Idents(s1717_3.seurat) <- "s1717_3"
Idents(s1815_2.seurat) <- "s1815_2"
Idents(s1815_4.seurat) <- "s1815_4"
Idents(s1822_1.seurat) <- "s1822_1"
Idents(s1822_2.seurat) <- "s1822_2"
Idents(s2281_1.seurat) <- "s2281_1"
Idents(s2281_3.seurat) <- "s2281_3"
Idents(s2392_1.seurat) <- "s2392_1"
Idents(s2392_3.seurat) <- "s2392_3"

s1356_A.seurat$dataset <- "s1356_A"
s1356_C.seurat$dataset <- "s1356_C"
s1521_4.seurat$dataset <- "s1521_4"
s1521_7.seurat$dataset <- "s1521_7"
s1295_1.seurat$dataset <- "s1295_1"
s1295_3.seurat$dataset <- "s1295_3"
s1717_2.seurat$dataset <- "s1717_2"
s1717_3.seurat$dataset <- "s1717_3"
s1815_2.seurat$dataset <- "s1815_2"
s1815_4.seurat$dataset <- "s1815_4"
s1822_1.seurat$dataset <- "s1822_1"
s1822_2.seurat$dataset <- "s1822_2"
s2281_1.seurat$dataset <- "s2281_1"
s2281_3.seurat$dataset <- "s2281_3"
s2392_1.seurat$dataset <- "s2392_1"
s2392_3.seurat$dataset <- "s2392_3"




# Calculate MT stats

s1356_A.seurat <- PercentageFeatureSet(s1356_A.seurat, pattern = "^MT-", col.name = "percent.mt")
s1356_C.seurat <- PercentageFeatureSet(s1356_C.seurat, pattern = "^MT-", col.name = "percent.mt")
s1521_4.seurat <- PercentageFeatureSet(s1521_4.seurat, pattern = "^MT-", col.name = "percent.mt")
s1521_7.seurat <- PercentageFeatureSet(s1521_7.seurat, pattern = "^MT-", col.name = "percent.mt")
s1295_1.seurat <- PercentageFeatureSet(s1295_1.seurat, pattern = "^MT-", col.name = "percent.mt")
s1295_3.seurat <- PercentageFeatureSet(s1295_3.seurat, pattern = "^MT-", col.name = "percent.mt")
s1717_2.seurat <- PercentageFeatureSet(s1717_2.seurat, pattern = "^MT-", col.name = "percent.mt")
s1717_3.seurat <- PercentageFeatureSet(s1717_3.seurat, pattern = "^MT-", col.name = "percent.mt")
s1815_2.seurat <- PercentageFeatureSet(s1815_2.seurat, pattern = "^MT-", col.name = "percent.mt")
s1815_4.seurat <- PercentageFeatureSet(s1815_4.seurat, pattern = "^MT-", col.name = "percent.mt")
s1822_1.seurat <- PercentageFeatureSet(s1822_1.seurat, pattern = "^MT-", col.name = "percent.mt")
s1822_2.seurat <- PercentageFeatureSet(s1822_2.seurat, pattern = "^MT-", col.name = "percent.mt")
s2281_1.seurat <- PercentageFeatureSet(s2281_1.seurat, pattern = "^MT-", col.name = "percent.mt")
s2281_3.seurat <- PercentageFeatureSet(s2281_3.seurat, pattern = "^MT-", col.name = "percent.mt")
s2392_1.seurat <- PercentageFeatureSet(s2392_1.seurat, pattern = "^MT-", col.name = "percent.mt")
s2392_3.seurat <- PercentageFeatureSet(s2392_3.seurat, pattern = "^MT-", col.name = "percent.mt")



# Merge

all.normalized <- merge(s1295_1.seurat, y = c(s1295_3.seurat,s1356_A.seurat,s1356_C.seurat,s1521_4.seurat,s1521_7.seurat,s1717_2.seurat,s1717_3.seurat,s1815_2.seurat,s1815_4.seurat,s1822_1.seurat,s1822_2.seurat,s2281_1.seurat,s2281_3.seurat,s2392_1.seurat,s2392_3.seurat), add.cell.ids=c("s1295_1","s1295_3","s1356_A","s1356_C","s1521_4","s1521_7","s1717_2","s1717_3","s1815_2","s1815_4","s1822_1","s1822_2","s2281_1","s2281_3","s2392_1","s2392_3"), project="Rahul_Akkem_20240604", merge.data=TRUE)

temp <- colnames(all.normalized)

temp1 <- grep("s1295_1",temp)
names(temp1) <- rep("s1295_1",length(temp1))

temp2 <- grep("s1295_3",temp)
names(temp2) <- rep("s1295_3",length(temp2))

temp3 <- grep("s1356_A",temp)
names(temp3) <- rep("s1356_A",length(temp3))

temp4 <- grep("s1356_C",temp)
names(temp4) <- rep("s1356_C",length(temp4))

temp5 <- grep("s1521_4",temp)
names(temp5) <- rep("s1521_4",length(temp5))

temp6 <- grep("s1521_7",temp)
names(temp6) <- rep("s1521_7",length(temp6))

temp7 <- grep("s1717_2",temp)
names(temp7) <- rep("s1717_2",length(temp7))

temp8 <- grep("s1717_3",temp)
names(temp8) <- rep("s1717_3",length(temp8))

temp9 <- grep("s1815_2",temp)
names(temp9) <- rep("s1815_2",length(temp9))

temp10 <- grep("s1815_4",temp)
names(temp10) <- rep("s1815_4",length(temp10))

temp11 <- grep("s1822_1",temp)
names(temp11) <- rep("s1822_1",length(temp11))

temp12 <- grep("s1822_2",temp)
names(temp12) <- rep("s1822_2",length(temp12))

temp13 <- grep("s2281_1",temp)
names(temp13) <- rep("s2281_1",length(temp13))

temp14 <- grep("s2281_3",temp)
names(temp14) <- rep("s2281_3",length(temp14))

temp15 <- grep("s2392_1",temp)
names(temp15) <- rep("s2392_1",length(temp15))

temp16 <- grep("s2392_3",temp)
names(temp16) <- rep("s2392_3",length(temp16))





temp <- names(sort(c(temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12,temp13,temp14,temp15,temp16)))
all.normalized@meta.data$orig.ident <- temp




# Inspect

svg("1_gex_primary_inspection_barplot_number_cells_per_sample.svg")
temp1 <- table(all.normalized@meta.data$orig.ident)
ccc <- hue_pal()(8)
par(mar=c(11,4,4,2))
barplot(as.numeric(temp1),ylab="nCells",names.arg = names(temp1),las = 2,col=ccc) + theme(axis.text.x = element_text(angle = 90))
abline(h=0)
dev.off()
ccc <- hue_pal()(8)
par(mar=c(11,4,4,2))
barplot(as.numeric(temp1),ylab="nCells",names.arg = names(temp1),las = 2,col=ccc) + theme(axis.text.x = element_text(angle = 90))
abline(h=0)







svg("2_gex_primary_inspection_violin_plot.svg")
VlnPlot(all.normalized, features = c("percent.mt","nCount_RNA","nFeature_RNA"), ncol = 3,pt.size=0)
dev.off()
VlnPlot(all.normalized, features = c("percent.mt","nCount_RNA","nFeature_RNA"), ncol = 3,pt.size=0)

plot1 <- FeatureScatter(all.normalized, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all.normalized, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")






svg("3_gex_primary_inspection_scatter_plots.svg")
plot1 + NoLegend() + plot2
dev.off()
plot1 + NoLegend() + plot2





svg("4_gex_primary_inspection_MT_density_plot.svg")
plot(density(as.numeric(all.normalized@meta.data$percent.mt)),main="percent.mt")
abline(v=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(1*sd(as.numeric(all.normalized@meta.data$percent.mt)))),lty=2)
(mean(as.numeric(all.normalized@meta.data$percent.mt))+(1*sd(as.numeric(all.normalized@meta.data$percent.mt))))
abline(v=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(2*sd(as.numeric(all.normalized@meta.data$percent.mt)))),lty=3)
(mean(as.numeric(all.normalized@meta.data$percent.mt))+(2*sd(as.numeric(all.normalized@meta.data$percent.mt))))
abline(v=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(3*sd(as.numeric(all.normalized@meta.data$percent.mt)))),lty=4)
(mean(as.numeric(all.normalized@meta.data$percent.mt))+(3*sd(as.numeric(all.normalized@meta.data$percent.mt))))
abline(v=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(4*sd(as.numeric(all.normalized@meta.data$percent.mt)))),lty=5)
(mean(as.numeric(all.normalized@meta.data$percent.mt))+(4*sd(as.numeric(all.normalized@meta.data$percent.mt))))
abline(v=5,lty=6)
legend("topright",c("1SD","2SD","3SD","4SD","5%"),lwd=2,lty=c(2,3,4,5,6))
dev.off()
plot(density(as.numeric(all.normalized@meta.data$percent.mt)),main="percent.mt")
abline(v=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(1*sd(as.numeric(all.normalized@meta.data$percent.mt)))),lty=2)
(mean(as.numeric(all.normalized@meta.data$percent.mt))+(1*sd(as.numeric(all.normalized@meta.data$percent.mt))))
abline(v=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(2*sd(as.numeric(all.normalized@meta.data$percent.mt)))),lty=3)
(mean(as.numeric(all.normalized@meta.data$percent.mt))+(2*sd(as.numeric(all.normalized@meta.data$percent.mt))))
abline(v=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(3*sd(as.numeric(all.normalized@meta.data$percent.mt)))),lty=4)
(mean(as.numeric(all.normalized@meta.data$percent.mt))+(3*sd(as.numeric(all.normalized@meta.data$percent.mt))))
abline(v=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(4*sd(as.numeric(all.normalized@meta.data$percent.mt)))),lty=5)
(mean(as.numeric(all.normalized@meta.data$percent.mt))+(4*sd(as.numeric(all.normalized@meta.data$percent.mt))))
abline(v=5,lty=6)
legend("topright",c("1SD","2SD","3SD","4SD","5%"),lwd=2,lty=c(2,3,4,5,6))

svg("5_gex_primary_inspection_MT_violin_plot_with_selected_cutoffs.svg")
VlnPlot(all.normalized, features = c("percent.mt"),pt.size = 0) + NoLegend() + geom_hline(yintercept=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(3*sd(as.numeric(all.normalized@meta.data$percent.mt)))),color="black",linetype=c(2)) + theme(axis.text.x = element_text(angle = 90))
dev.off()
VlnPlot(all.normalized, features = c("percent.mt"),pt.size = 0) + NoLegend() + geom_hline(yintercept=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(3*sd(as.numeric(all.normalized@meta.data$percent.mt)))),color="black",linetype=c(2)) + theme(axis.text.x = element_text(angle = 90))

svg("6_gex_primary_inspection_ncount_density_plot.svg")
plot(density(as.numeric(all.normalized@meta.data$nCount_RNA)),main="nCount_RNA")
abline(v=(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(1*sd(as.numeric(all.normalized@meta.data$nCount_RNA)))),lty=2)
(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(1*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(2*sd(as.numeric(all.normalized@meta.data$nCount_RNA)))),lty=3)
(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(2*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(3*sd(as.numeric(all.normalized@meta.data$nCount_RNA)))),lty=4)
(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(3*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nCount_RNA)))),lty=5)
(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))
abline(v=500,lty=1)
legend("topright",c("250","1SD","2SD","3SD","4SD"),lwd=2,lty=c(1,2,3,4,5))
dev.off()
plot(density(as.numeric(all.normalized@meta.data$nCount_RNA)),main="nCount_RNA")
abline(v=(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(1*sd(as.numeric(all.normalized@meta.data$nCount_RNA)))),lty=2)
(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(1*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(2*sd(as.numeric(all.normalized@meta.data$nCount_RNA)))),lty=3)
(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(2*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(3*sd(as.numeric(all.normalized@meta.data$nCount_RNA)))),lty=4)
abline(v=(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nCount_RNA)))),lty=5)
(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))
abline(v=500,lty=1)
legend("topright",c("250","1SD","2SD","3SD","4SD"),lwd=2,lty=c(1,2,3,4,5))

svg("7_gex_primary_inspection_ncount_violin_plot_with_selected_cutoffs.svg")
VlnPlot(all.normalized, features = "nCount_RNA",pt.size = 0) + NoLegend() + geom_hline(yintercept=c(500,(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))),color="black",linetype=c(1,3)) + theme(axis.text.x = element_text(angle = 90))
dev.off()
VlnPlot(all.normalized, features = "nCount_RNA",pt.size = 0) + NoLegend() + geom_hline(yintercept=c(500,(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))),color="black",linetype=c(1,3)) + theme(axis.text.x = element_text(angle = 90))

svg("8_gex_primary_inspection_nfeatures_density_plot.svg")
plot(density(as.numeric(all.normalized@meta.data$nFeature_RNA)),main="nFeature_RNA")
abline(v=(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(1*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),lty=2)
(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(1*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(2*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),lty=3)
(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(2*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(3*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),lty=4)
(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(3*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),lty=5)
(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))
abline(v=250,lty=1)
legend("topright",c("250","1SD","2SD","3SD","4SD"),lwd=2,lty=c(1,2,3,4,5))
dev.off()
plot(density(as.numeric(all.normalized@meta.data$nFeature_RNA)),main="nFeature_RNA")
abline(v=(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(1*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),lty=2)
(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(1*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(2*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),lty=3)
(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(2*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(3*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),lty=4)
(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(3*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))
abline(v=(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),lty=5)
(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))
abline(v=250,lty=1)
legend("topright",c("250","1SD","2SD","3SD","4SD"),lwd=2,lty=c(1,2,3,4,5))

svg("9_gex_primary_inspection_nfeatures_violin_plot_with_selected_cutoffs.svg")
VlnPlot(all.normalized, features = "nFeature_RNA",pt.size = 0) + NoLegend() + geom_hline(yintercept=c(250,mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),color="black",linetype=c(1,3)) + theme(axis.text.x = element_text(angle = 90))
dev.off()
VlnPlot(all.normalized, features = "nFeature_RNA",pt.size = 0) + NoLegend() + geom_hline(yintercept=c(250,mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))),color="black",linetype=c(1,3)) + theme(axis.text.x = element_text(angle = 90))

# Filter

svg("10_gex_pre_hard_filtering_scatterplot_with_selected_cutoffs.svg")
plot1 + geom_hline(yintercept=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(3*sd(as.numeric(all.normalized@meta.data$percent.mt)))),color="black",linetype=2) + NoLegend() + plot2 + geom_hline(yintercept=c(250,(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))),color="black",linetype=2) + geom_vline(xintercept=c(500,(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))),color="black",linetype=2)
dev.off()
plot1 + geom_hline(yintercept=(mean(as.numeric(all.normalized@meta.data$percent.mt))+(3*sd(as.numeric(all.normalized@meta.data$percent.mt)))),color="black",linetype=2) + NoLegend() + plot2 + geom_hline(yintercept=c(250,(mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nFeature_RNA))))),color="black",linetype=2) + geom_vline(xintercept=c(500,(mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nCount_RNA))))),color="black",linetype=2)

# MT Filter = 2.787142
(mean(as.numeric(all.normalized@meta.data$percent.mt))+(3*sd(as.numeric(all.normalized@meta.data$percent.mt))))

# nFeature Filter = 3307.526.   or 5061.235
mean(as.numeric(all.normalized@meta.data$nFeature_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nFeature_RNA)))

# nCount Filter = 8916.1     or 15988.05
mean(as.numeric(all.normalized@meta.data$nCount_RNA))+(4*sd(as.numeric(all.normalized@meta.data$nCount_RNA)))

all.filtered <- subset(all.normalized, subset = nCount_RNA > 500 & nFeature_RNA > 250 & nFeature_RNA < 6209.353 & nCount_RNA < 20259.15 & percent.mt <  16.00667)

temp <- rownames(all.filtered)
temp <- grep("^mt-",temp)
all.filtered <- all.filtered[-c(temp),]

plot1 <- FeatureScatter(all.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

svg("11_gex_post_hard_filtering_scatterplot_with_selected_cutoffs.svg")
plot2
dev.off()
plot2

svg("12_gex_post_hard_filtering_barplot_number_cells_per_sample.svg")
temp1 <- table(all.filtered@meta.data$orig.ident)
ccc <- hue_pal()(8)
par(mar=c(11,4,4,2))
barplot(as.numeric(temp1),ylab="nCells",names.arg = names(temp1),las = 2,col=ccc) + theme(axis.text.x = element_text(angle = 90))
abline(h=0)
dev.off()
ccc <- hue_pal()(8)
par(mar=c(11,4,4,2))
barplot(as.numeric(temp1),ylab="nCells",names.arg = names(temp1),las = 2,col=ccc) + theme(axis.text.x = element_text(angle = 90))
abline(h=0)

# identify and remove doublet cells by sample

all.list <- SplitObject(all.filtered, split.by = "orig.ident")
running.keep <- NULL
for(i in 1:length(all.list)) {
  print(i)
  temp.grab <- all.list[[i]]
  genes.to.use <- rownames(temp.grab)
  DefaultAssay(temp.grab) <- "RNA"
  temp.grab.sce <- as.SingleCellExperiment(temp.grab)
  temp.grab.sce <- temp.grab.sce[genes.to.use,]
  set.seed(123)
  dbl.dens <- computeDoubletDensity(temp.grab.sce, d=100)
  temp.grab$DoubletScore <- log10(dbl.dens+1)
  temp1 <- temp.grab$DoubletScore
  names(temp1) <- colnames(temp.grab.sce)
  temp2 <- as.numeric(temp1)>(mean(as.numeric(temp1))+(2*sd(as.numeric(temp1))))
  names(temp2) <- names(temp1)
  temp3 <- temp2[temp2=="FALSE"]
  running.keep <- c(running.keep,names(temp3))
}

all.curated <- all.filtered[,running.keep]

svg("13_gex_post_doublet_cell_filtering_barplot_number_cells_per_sample.svg")
temp1 <- table(all.curated@meta.data$orig.ident)
ccc <- hue_pal()(8)
par(mar=c(11,4,4,2))
barplot(as.numeric(temp1),ylab="nCells",names.arg = names(temp1),las = 2,col=ccc) + theme(axis.text.x = element_text(angle = 90))
abline(h=0)
dev.off()
temp1 <- table(all.curated@meta.data$orig.ident)
ccc <- hue_pal()(8)
par(mar=c(11,4,4,2))
barplot(as.numeric(temp1),ylab="nCells",names.arg = names(temp1),las = 2,col=ccc) + theme(axis.text.x = element_text(angle = 90))
abline(h=0)

# integrate via harmony

harmony.go <- Seurat::NormalizeData(all.curated,verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(vars.to.regress="percent.mt",verbose = FALSE) %>% 
  RunPCA(pc.genes = harmony.go@var.genes, npcs = 200, verbose = FALSE)
p <- ElbowPlot(harmony.go,ndims=200)
svg("14_pre_harmony_integration_scree_plot.svg")
p
dev.off()
p

harmony.go <- RunPCA(harmony.go, pc.genes = harmony.go@var.genes, npcs = 75, verbose = FALSE) %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:75)

svg("15_post_harmony_integration_umap_plot_with_sample_overlaid.svg")
DimPlot(harmony.go, label=FALSE, group.by="orig.ident")
dev.off()
DimPlot(harmony.go, label=FALSE, group.by="orig.ident")

svg("16_post_harmony_integration_umap_plot_split_out_by_sample.svg",width=20)
DimPlot(harmony.go, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident",ncol=4)
dev.off()
DimPlot(harmony.go, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident",ncol=4)

# perform initial clustering and inspect by clustree

harmony.go <- harmony.go %>% FindNeighbors(harmony.go, reduction = "harmony", dims = 1:75, verbose = FALSE) %>% 
  FindClusters(resolution=c(0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5),save.SNN=TRUE)

svg("17_initial_clustering_clustree_results.svg",width=10,height=25)
clustree(harmony.go,prefix="RNA_snn_res.")
dev.off()
clustree(harmony.go,prefix="RNA_snn_res.")



pdf("18_initial_clustering_umap_plots_per_resolution.pdf")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.04)))))
harmony.go@meta.data$RNA_snn_res.0.04 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.04, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.04")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.05)))))
harmony.go@meta.data$RNA_snn_res.0.05 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.05, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.05")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.06)))))
harmony.go@meta.data$RNA_snn_res.0.06 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.06, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.06")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.07)))))
harmony.go@meta.data$RNA_snn_res.0.07 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.07, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.07")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.08)))))
harmony.go@meta.data$RNA_snn_res.0.08 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.08, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.08")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.09)))))
harmony.go@meta.data$RNA_snn_res.0.09 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.09, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.09")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.1)))))
harmony.go@meta.data$RNA_snn_res.0.1 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.1, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.1")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.2)))))
harmony.go@meta.data$RNA_snn_res.0.2 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.2, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.2")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.3)))))
harmony.go@meta.data$RNA_snn_res.0.3 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.3, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.3")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.4)))))
harmony.go@meta.data$RNA_snn_res.0.4 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.4, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.4")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.5)))))
harmony.go@meta.data$RNA_snn_res.0.5 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.5, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.5")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.6)))))
harmony.go@meta.data$RNA_snn_res.0.6 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.6, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.6")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.7)))))
harmony.go@meta.data$RNA_snn_res.0.7 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.7, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.7")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.8)))))
harmony.go@meta.data$RNA_snn_res.0.8 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.8, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.8")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.0.9)))))
harmony.go@meta.data$RNA_snn_res.0.9 <- factor(x=harmony.go@meta.data$RNA_snn_res.0.9, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.0.9")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1)))))
harmony.go@meta.data$RNA_snn_res.1 <- factor(x=harmony.go@meta.data$RNA_snn_res.1, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1.1)))))
harmony.go@meta.data$RNA_snn_res.1.1 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.1, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.1")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1.2)))))
harmony.go@meta.data$RNA_snn_res.1.2 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.2, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.2")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1.3)))))
harmony.go@meta.data$RNA_snn_res.1.3 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.3, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.3")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1.4)))))
harmony.go@meta.data$RNA_snn_res.1.4 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.4, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.4")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1.5)))))
harmony.go@meta.data$RNA_snn_res.1.5 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.5, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.5")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1.6)))))
harmony.go@meta.data$RNA_snn_res.1.6 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.6, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.6")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1.7)))))
harmony.go@meta.data$RNA_snn_res.1.7 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.7, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.7")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1.8)))))
harmony.go@meta.data$RNA_snn_res.1.8 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.8, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.8")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.1.9)))))
harmony.go@meta.data$RNA_snn_res.1.9 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.9, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.9")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2)))))
harmony.go@meta.data$RNA_snn_res.2 <- factor(x=harmony.go@meta.data$RNA_snn_res.2, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2.1)))))
harmony.go@meta.data$RNA_snn_res.2.1 <- factor(x=harmony.go@meta.data$RNA_snn_res.2.1, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2.1")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2.2)))))
harmony.go@meta.data$RNA_snn_res.2.2 <- factor(x=harmony.go@meta.data$RNA_snn_res.2.2, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2.2")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2.3)))))
harmony.go@meta.data$RNA_snn_res.2.3 <- factor(x=harmony.go@meta.data$RNA_snn_res.2.3, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2.3")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2.4)))))
harmony.go@meta.data$RNA_snn_res.2.4 <- factor(x=harmony.go@meta.data$RNA_snn_res.2.4, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2.4")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2.5)))))
harmony.go@meta.data$RNA_snn_res.2.5 <- factor(x=harmony.go@meta.data$RNA_snn_res.2.5, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2.5")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2.6)))))
harmony.go@meta.data$RNA_snn_res.2.6 <- factor(x=harmony.go@meta.data$RNA_snn_res.2.6, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2.6")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2.7)))))
harmony.go@meta.data$RNA_snn_res.2.7 <- factor(x=harmony.go@meta.data$RNA_snn_res.2.7, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2.7")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2.8)))))
harmony.go@meta.data$RNA_snn_res.2.8 <- factor(x=harmony.go@meta.data$RNA_snn_res.2.8, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2.8")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.2.9)))))
harmony.go@meta.data$RNA_snn_res.2.9 <- factor(x=harmony.go@meta.data$RNA_snn_res.2.9, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.2.9")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.3)))))
harmony.go@meta.data$RNA_snn_res.3 <- factor(x=harmony.go@meta.data$RNA_snn_res.3, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.3")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.3.1)))))
harmony.go@meta.data$RNA_snn_res.3.1 <- factor(x=harmony.go@meta.data$RNA_snn_res.3.1, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.3.1")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.3.2)))))
harmony.go@meta.data$RNA_snn_res.3.2 <- factor(x=harmony.go@meta.data$RNA_snn_res.3.2, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.3.2")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.3.3)))))
harmony.go@meta.data$RNA_snn_res.3.3 <- factor(x=harmony.go@meta.data$RNA_snn_res.3.3, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.3.3")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.3.4)))))
harmony.go@meta.data$RNA_snn_res.3.4 <- factor(x=harmony.go@meta.data$RNA_snn_res.3.4, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.3.4")

my_levels <- as.character(c(0:length(table(as.character(harmony.go@meta.data$RNA_snn_res.3.5)))))
harmony.go@meta.data$RNA_snn_res.3.5 <- factor(x=harmony.go@meta.data$RNA_snn_res.3.5, levels=my_levels)
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.3.5")

dev.off()

svg("19_gex_post_initial_clustering_barplot_number_cells_per_cluster_for_res_1.6.svg")
temp3 <- table(as.character(harmony.go@meta.data$RNA_snn_res.1.6))
temp4 <- as.data.frame(cbind(names(temp3),as.numeric(temp3)))
temp5 <- temp4[,2]
names(temp5) <- temp4[,1]
temp5 <- temp5[as.character(c(0:38))]
ccc <- hue_pal()(length(temp5))
barplot(as.numeric(temp5),col=ccc,xlab="Cluster",ylab="nCells",names.arg = names(temp5),las = 2)
abline(h=0)
dev.off()
temp3 <- table(as.character(harmony.go@meta.data$RNA_snn_res.1.6))
temp4 <- as.data.frame(cbind(names(temp3),as.numeric(temp3)))
temp5 <- temp4[,2]
names(temp5) <- temp4[,1]
temp5 <- temp5[as.character(c(0:38))]
ccc <- hue_pal()(length(temp5))
barplot(as.numeric(temp5),col=ccc,xlab="Cluster",ylab="nCells",names.arg = names(temp5),las = 2)
abline(h=0)

svg("20_gex_post_initial_clustering_umap_for_res_1.6.svg")
harmony.go@meta.data$RNA_snn_res.1.6 <- factor(x=harmony.go@meta.data$RNA_snn_res.1.6, levels=as.character(c(0:38)))
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.6")
dev.off()
DimPlot(harmony.go, label=TRUE,group.by="RNA_snn_res.1.6",repel=T)




# recode

temp <- as.character(unlist(harmony.go@meta.data$RNA_snn_res.1.6))
temp <- ifelse(temp=="0","1",temp)
temp <- ifelse(temp=="3","1",temp)
temp <- ifelse(temp=="2","1",temp)
temp <- ifelse(temp=="10","1",temp)
temp <- ifelse(temp=="12","1",temp)
temp <- ifelse(temp=="5","4",temp)
temp <- ifelse(temp=="26","18",temp)

#c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39")
#harmony.go@meta.data$initial.clusters <- factor(temp,levels=as.character(c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39")))
harmony.go@meta.data$initial.clusters <- factor(temp,levels=as.character(c("1","4","6","7","8","9","11","13","14","15","16","17","18","19","20","21","22","23","24","25","27","28","29","30","31","32","33","34","35","36","37")))

Idents(harmony.go) <- temp











svg("21_gex_post_initial_clustering_umap_for_res_1.6_recoded.svg")
DimPlot(harmony.go, label=TRUE,group.by="initial.clusters")
dev.off()
DimPlot(harmony.go, label=TRUE,group.by="initial.clusters")

# identify and remove doublet clusters
# https://github.com/plger/scDblFinder
# https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/findDoubletClusters.html

DefaultAssay(harmony.go) <- "RNA"
harmony.go.joinedlayers <- JoinLayers(harmony.go) # Had to create a joined layers because of combined seurat harmony wont work in the SingleCellExperiment and will split it in the next line
harmony.go.sce <- as.SingleCellExperiment(harmony.go.joinedlayers)
set.seed(123)
dbl.out <- findDoubletClusters(x=harmony.go.sce, clusters=colData(harmony.go.sce)$initial.clusters,assay.type="counts")
chosen.doublet <- rownames(dbl.out)[isOutlier(as.numeric(dbl.out$num.de), type="lower", log=TRUE)]
chosen.doublet

# doublet clusters = none

svg("22_gex_post_doublet_cluster_filtering_barplot_for_res_1.6_recoded.svg")
temp <- dbl.out$num.de
names(temp) <- rownames(dbl.out)
m <- mean(temp)
s <- sd(temp)
temp <- temp-m
temp <- temp/s
temp1 <- min(c(-3,temp))
temp2 <- max(c(3,temp))
barplot(rev(sort(temp)),xlab="Cluster",ylab="Z-score (Number Differential Genes)",las = 2, cex.names = 0.7,space=0, ylim = c(temp1,temp2),col=c(rep(hue_pal()(2)[1],length(rownames(dbl.out))-length(chosen.doublet)),rep(hue_pal()(2)[2],length(chosen.doublet))),main="Doublet Cluster")
abline(h=0)
box()
legend("topright", legend = c("False","True"), fill = c(hue_pal()(2)[1],hue_pal()(2)[2]))
abline(h=(-2),lty=2)
dev.off()
temp <- dbl.out$num.de
names(temp) <- rownames(dbl.out)
m <- mean(temp)
s <- sd(temp)
temp <- temp-m
temp <- temp/s
temp1 <- min(c(-3,temp))
temp2 <- max(c(3,temp))
barplot(rev(sort(temp)),xlab="Cluster",ylab="Z-score (Number Differential Genes)",las = 2, cex.names = 0.7,space=0, ylim = c(temp1,temp2),col=c(rep(hue_pal()(2)[1],length(rownames(dbl.out))-length(chosen.doublet)),rep(hue_pal()(2)[2],length(chosen.doublet))),main="Doublet Cluster")
abline(h=0)
box()
legend("topright", legend = c("False","True"), fill = c(hue_pal()(2)[1],hue_pal()(2)[2]))
abline(h=(-2),lty=2)

curated.harmony.go <- harmony.go
rm(harmony.go)
curated.harmony.go@meta.data$doubletcluster <- ifelse(curated.harmony.go@meta.data$initial.clusters%in%chosen.doublet,"True","False")

svg("23_gex_post_doublet_cluster_filtering_umap_for_res_1.6_recoded.svg")
my_levels <- as.character(c(0:length(table(as.character(curated.harmony.go@meta.data$initial.clusters)))))
DimPlot(curated.harmony.go,group.by="doubletcluster")
dev.off()
DimPlot(curated.harmony.go,group.by="doubletcluster")

if(length(chosen.doublet)>0) {
  temp1 <- chosen.doublet
  for (i in 1:length(temp1)) {
    temp2 <- curated.harmony.go@meta.data$initial.clusters!=temp1[i]
    curated.harmony.go <- curated.harmony.go[,temp2]
  }
}

# identify doublet cells
# https://github.com/plger/scDblFinder
# https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/computeDoubletDensity.html

DefaultAssay(curated.harmony.go) <- "RNA"
curated.harmony.go.joinedlayers=JoinLayers(curated.harmony.go)
curated.harmony.go.sce <- as.SingleCellExperiment(curated.harmony.go.joinedlayers)
set.seed(123)
dbl.dens <- computeDoubletDensity(curated.harmony.go.sce, d=50)
curated.harmony.go.sce$DoubletScore <- log10(dbl.dens+1)
temp1 <- curated.harmony.go.sce$DoubletScore
names(temp1) <- colnames(curated.harmony.go.sce)
temp1 <- temp1[colnames(curated.harmony.go)]
curated.harmony.go@meta.data$DoubletScore <- as.numeric(temp1)
curated.harmony.go@meta.data$DoubletScoreCoded <- curated.harmony.go@meta.data$DoubletScore>(mean(as.numeric(curated.harmony.go@meta.data$DoubletScore))+(2*sd(as.numeric(curated.harmony.go@meta.data$DoubletScore))))

svg("24_gex_post_doublet_cell_filtering_density_plot_for_res_1.6_recoded.svg")
plot(density(as.numeric(temp1)),main="DoubletScore")
abline(v=(mean(as.numeric(curated.harmony.go@meta.data$DoubletScore))+(2*sd(as.numeric(curated.harmony.go@meta.data$DoubletScore)))),lty=2)
legend("topright", legend = c("2SD"),lty=2)
dev.off()
plot(density(as.numeric(temp1)),main="DoubletScore")
abline(v=(mean(as.numeric(curated.harmony.go@meta.data$DoubletScore))+(2*sd(as.numeric(curated.harmony.go@meta.data$DoubletScore)))),lty=2)
legend("topright", legend = c("2SD"),lty=2)

svg("25_gex_post_doublet_cell_filtering_umap_for_res_1.6_recoded.svg")
DimPlot(curated.harmony.go, group.by="DoubletScoreCoded",label=FALSE,cols=c(hue_pal()(2)))
dev.off()
DimPlot(curated.harmony.go, group.by="DoubletScoreCoded",label=FALSE,cols=c(hue_pal()(2)))

svg("26_gex_post_doublet_cell_filtering_umap_for_res_1.6_recoded_split_by_True_False.svg")
DimPlot(curated.harmony.go, split.by="DoubletScoreCoded", group.by="DoubletScoreCoded",label=FALSE,cols=c(hue_pal()(2)))
dev.off()
DimPlot(curated.harmony.go, split.by="DoubletScoreCoded", group.by="DoubletScoreCoded",label=FALSE,cols=c(hue_pal()(2)))

final.curated.harmony.go <- curated.harmony.go[,curated.harmony.go@meta.data$DoubletScoreCoded==FALSE]
rm(curated.harmony.go)

# perform final clustering and inspect by clustree

final.curated.harmony.go <- final.curated.harmony.go %>% FindNeighbors(final.curated.harmony.go, reduction = "harmony", dims = 1:75, verbose = FALSE) %>% FindClusters(resolution=c(0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5),save.SNN=TRUE)

svg("27_final_clustering_clustree_results.svg",width=10,height=25)
clustree(final.curated.harmony.go,prefix="RNA_snn_res.")
dev.off()
clustree(final.curated.harmony.go,prefix="RNA_snn_res.")

pdf("28_final_clustering_umap_plots_per_resolution.pdf")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.04)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.04 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.04, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.04")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.05)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.05 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.05, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.05")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.06)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.06 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.06, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.06")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.07)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.07 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.07, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.07")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.08)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.08 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.08, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.08")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.09)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.09 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.09, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.09")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.1)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.1 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.1, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.1")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.2)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.2 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.2, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.2")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.3)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.3 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.3, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.3")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.4)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.4 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.4, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.4")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.5)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.5 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.5, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.5")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.6)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.6 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.6, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.6")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.7)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.7 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.7, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.7")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.8)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.8 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.8, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.8")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.0.9)))))
final.curated.harmony.go@meta.data$RNA_snn_res.0.9 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.0.9, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.0.9")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.1)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1.1 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.1, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.1")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.2)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1.2 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.2, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.2")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.3)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1.3 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.3, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.3")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.4)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1.4 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.4, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.4")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.5)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1.5 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.5, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.5")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.6)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1.6 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.6, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.6")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.7)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1.7 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.7, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.7")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.8)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1.8 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.8, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.8")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.9)))))
final.curated.harmony.go@meta.data$RNA_snn_res.1.9 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.9, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.9")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2.1)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2.1 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2.1, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2.1")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2.2)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2.2 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2.2, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2.2")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2.3)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2.3 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2.3, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2.3")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2.4)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2.4 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2.4, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2.4")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2.5)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2.5 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2.5, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2.5")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2.6)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2.6 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2.6, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2.6")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2.7)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2.7 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2.7, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2.7")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2.8)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2.8 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2.8, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2.8")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.2.9)))))
final.curated.harmony.go@meta.data$RNA_snn_res.2.9 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.2.9, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.2.9")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.3)))))
final.curated.harmony.go@meta.data$RNA_snn_res.3 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.3, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.3")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.3.1)))))
final.curated.harmony.go@meta.data$RNA_snn_res.3.1 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.3.1, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.3.1")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.3.2)))))
final.curated.harmony.go@meta.data$RNA_snn_res.3.2 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.3.2, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.3.2")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.3.3)))))
final.curated.harmony.go@meta.data$RNA_snn_res.3.3 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.3.3, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.3.3")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.3.4)))))
final.curated.harmony.go@meta.data$RNA_snn_res.3.4 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.3.4, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.3.4")

my_levels <- as.character(c(0:length(table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.3.5)))))
final.curated.harmony.go@meta.data$RNA_snn_res.3.5 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.3.5, levels=my_levels)
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.3.5")

dev.off()

svg("29_gex_post_final_clustering_barplot_number_cells_per_cluster_for_res_1.7.svg")
temp3 <- table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.7))
temp4 <- as.data.frame(cbind(names(temp3),as.numeric(temp3)))
temp5 <- temp4[,2]
names(temp5) <- temp4[,1]
temp5 <- temp5[as.character(c(0:40))]
ccc <- hue_pal()(length(temp5))
barplot(as.numeric(temp5),col=ccc,xlab="Cluster",ylab="nCells",names.arg = names(temp5),las = 2)
abline(h=0)
dev.off()
temp3 <- table(as.character(final.curated.harmony.go@meta.data$RNA_snn_res.1.7))
temp4 <- as.data.frame(cbind(names(temp3),as.numeric(temp3)))
temp5 <- temp4[,2]
names(temp5) <- temp4[,1]
temp5 <- temp5[as.character(c(0:40))]
ccc <- hue_pal()(length(temp5))
barplot(as.numeric(temp5),col=ccc,xlab="Cluster",ylab="nCells",names.arg = names(temp5),las = 2)
abline(h=0)

svg("30_gex_post_final_clustering_umap_for_res_1.7.svg")
final.curated.harmony.go@meta.data$RNA_snn_res.1.7 <- factor(x=final.curated.harmony.go@meta.data$RNA_snn_res.1.7, levels=as.character(c(0:40))) # changed from 0:29 to 0:40
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.7")
dev.off()
DimPlot(final.curated.harmony.go, label=TRUE,group.by="RNA_snn_res.1.7",repel=TRUE)

# recode to combine some clusters

temp <- as.character(unlist(final.curated.harmony.go@meta.data$RNA_snn_res.1.7))
temp <- ifelse(temp=="0","1",temp)
temp <- ifelse(temp=="2","1",temp)
temp <- ifelse(temp=="13","1",temp)
temp <- ifelse(temp=="27","1",temp)
temp <- ifelse(temp=="11","1",temp)
temp <- ifelse(temp=="8","1",temp)
temp <- ifelse(temp=="12","7",temp)
temp <- ifelse(temp=="25","17",temp)
temp <- ifelse(temp=="3","4",temp)



#final.curated.harmony.go@meta.data$final.clusters <- factor(temp,levels=as.character( c("0","1","5","8","9","10","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","29")))
#Idents(final.curated.harmony.go) <- factor(temp,levels=as.character(c("0","1","5","8","9","10","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","29")))

c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39")
final.curated.harmony.go@meta.data$final.clusters <- factor(temp,levels=as.character( c("1","4","5","6","7","9","10","14","15","16","17","18","19","20","21","22","23","24","26","28","29","30","31","32","33","34","35","36","37")))
Idents(final.curated.harmony.go) <- factor(temp,levels=as.character(c("1","4","5","6","7","9","10","14","15","16","17","18","19","20","21","22","23","24","26","28","29","30","31","32","33","34","35","36","37")))




svg("31_gex_post_final_clustering_umap_for_res_1.7_recoded.svg")
DimPlot(final.curated.harmony.go, label=TRUE,group.by="final.clusters")
dev.off()
DimPlot(final.curated.harmony.go, label=TRUE,group.by="final.clusters")


# perform cluster by cluster pruning, total clusters = 23

to.do <- names((table(Idents(final.curated.harmony.go))))
length(to.do)

to.remove <- NULL

i <- 1
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(3.75),linetype='dashed') + geom_hline(yintercept=c(3.75),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>3.75
temp10[,2] <- temp10[,2]>3.75
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 2
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(3),linetype='dashed') + geom_hline(yintercept=c(3.75),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>3
temp10[,2] <- temp10[,2]>3.75
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 3
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(3),linetype='dashed') + geom_hline(yintercept=c(3),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>3
temp10[,2] <- temp10[,2]>3
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 4
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2.5),linetype='dashed') + geom_hline(yintercept=c(3.5),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2.5
temp10[,2] <- temp10[,2]>3.5
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 5
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(3),linetype='dashed') + geom_hline(yintercept=c(2.5),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>3
temp10[,2] <- temp10[,2]>2.5
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 6
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(5),linetype='dashed') + geom_hline(yintercept=c(4),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>5
temp10[,2] <- temp10[,2]>4
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 7
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2.5),linetype='dashed') + geom_hline(yintercept=c(2.5),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2.5
temp10[,2] <- temp10[,2]>2.5
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 8
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2),linetype='dashed') + geom_hline(yintercept=c(3),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2
temp10[,2] <- temp10[,2]>3
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 9
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(3.5),linetype='dashed') + geom_hline(yintercept=c(2),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>3.5
temp10[,2] <- temp10[,2]>2
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 10
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(3.5),linetype='dashed') + geom_hline(yintercept=c(4),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>3.5
temp10[,2] <- temp10[,2]>4
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 11
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(1),linetype='dashed') + geom_hline(yintercept=c(2),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>1
temp10[,2] <- temp10[,2]>2
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)


i <- 12
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2),linetype='dashed') + geom_hline(yintercept=c(2),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2
temp10[,2] <- temp10[,2]>2
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 13
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2.5),linetype='dashed') + geom_hline(yintercept=c(3),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2.5
temp10[,2] <- temp10[,2]>3
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 14
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(3),linetype='dashed') + geom_hline(yintercept=c(2.5),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>3
temp10[,2] <- temp10[,2]>2.5
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],"_first_pass.svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

#i <- 14
#temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
#temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
#temp4 <- temp3$data[,1]
#temp5 <- temp3$data[,2]
#X <- as.matrix(cbind(temp4,temp5))
#dimnames(X)[[1]] <- colnames(temp2)
#dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
#p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(-1.5),linetype='dashed') + geom_vline(xintercept=c(-0.5),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Supervised Revisit"),sep="",collapse=""))
#p1
#temp10 <- X
#temp10a <- temp10[,1]<(-0.5)
#temp10b <- temp10[,1]>(-1.5)
#temp11 <- temp10a+temp10b
#temp11 <- dimnames(temp10)[[1]][temp11<2]
#temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
#p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
#plot_lst <- vector("list", length = 2)
#plot_lst[[1]] <- p1
#plot_lst[[2]] <- p2
#cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
#svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],"_second_pass.svg"),sep="",collapse=""),width = 14)
#cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
#dev.off()
#to.remove <- c(to.remove,temp11)

i <- 15
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2),linetype='dashed') + geom_hline(yintercept=c(2),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2
temp10[,2] <- temp10[,2]>2
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 16
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(1),linetype='dashed') + geom_hline(yintercept=c(0.75),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>1
temp10[,2] <- temp10[,2]>0.75
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 17
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2.5),linetype='dashed') + geom_hline(yintercept=c(3),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2.5
temp10[,2] <- temp10[,2]>3
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 18
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2),linetype='dashed') + geom_hline(yintercept=c(2.5),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2
temp10[,2] <- temp10[,2]>2.5
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 19
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2),linetype='dashed') + geom_hline(yintercept=c(2.5),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2
temp10[,2] <- temp10[,2]>2.5
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 20
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2),linetype='dashed') + geom_hline(yintercept=c(1.5),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2
temp10[,2] <- temp10[,2]>1.5
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 21
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2),linetype='dashed') + geom_hline(yintercept=c(3),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2
temp10[,2] <- temp10[,2]>3
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 22
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(1.25),linetype='dashed') + geom_hline(yintercept=c(1),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>1.25
temp10[,2] <- temp10[,2]>1
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 23
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2.5),linetype='dashed') + geom_hline(yintercept=c(3),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2.5
temp10[,2] <- temp10[,2]>3
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 24
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(3),linetype='dashed') + geom_hline(yintercept=c(3),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>3
temp10[,2] <- temp10[,2]>3
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 25
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(2.5),linetype='dashed') + geom_hline(yintercept=c(2),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>2.5
temp10[,2] <- temp10[,2]>2
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 26
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(1.5),linetype='dashed') + geom_hline(yintercept=c(2),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>1.5
temp10[,2] <- temp10[,2]>2
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 27
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(3),linetype='dashed') + geom_hline(yintercept=c(3),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>3
temp10[,2] <- temp10[,2]>3
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 28
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(1.5),linetype='dashed') + geom_hline(yintercept=c(1),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>1.5
temp10[,2] <- temp10[,2]>1
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)

i <- 29
temp2 <- final.curated.harmony.go[,(final.curated.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
temp6 <- mean(temp4,trim=0.20)
temp7 <- mean(temp5,trim=0.20)
temp44 <- abs((temp4-temp6)/sd(temp4))
temp55 <- abs((temp5-temp7)/sd(temp5))
X <- as.matrix(cbind(temp44,temp55))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(1),linetype='dashed') + geom_hline(yintercept=c(1),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Zscore-based Outlier Threshold Selection"),sep="",collapse=""))
p1
temp10 <- X
temp10[,1] <- temp10[,1]>1
temp10[,2] <- temp10[,2]>1
temp11 <- dimnames(temp10)[[1]][apply(temp10,1,sum)>0]
temp2@meta.data$outliers <- rownames(temp2@meta.data)%in%temp11
p2 <- DimPlot(temp2, label=FALSE,group.by="outliers") + ggtitle(paste(c("Cluster ",to.do[i],": Outlier Cells To Remove"),sep="",collapse=""))
plot_lst <- vector("list", length = 2)
plot_lst[[1]] <- p1
plot_lst[[2]] <- p2
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
svg(paste(c("32_outlier_cell_curation_for_cluster_",to.do[i],".svg"),sep="",collapse=""),width = 14)
cowplot::plot_grid(plotlist = plot_lst, nrow = 1)
dev.off()
to.remove <- c(to.remove,temp11)


final.pruned.harmony.go <- final.curated.harmony.go
final.pruned.harmony.go@meta.data$pruned <- colnames(final.curated.harmony.go)%in%to.remove

svg("33_final_umap_post_cell_curation_for_res_1.7_recoded.svg")
DimPlot(final.pruned.harmony.go, label=FALSE,group.by="pruned")
dev.off()
DimPlot(final.pruned.harmony.go, label=FALSE,group.by="pruned")

final.pruned.harmony.go <- final.pruned.harmony.go[,final.pruned.harmony.go@meta.data$pruned=="FALSE"]

# Characterize

svg("34_final_umap_by_cluster_v1_pre_cluster_27_split.svg")
DimPlot(final.pruned.harmony.go, label=FALSE,group.by="final.clusters")
dev.off()
DimPlot(final.pruned.harmony.go, label=FALSE,group.by="final.clusters")

svg("34_final_umap_by_cluster_v2_pre_cluster_27_split.svg")
DimPlot(final.pruned.harmony.go, label=TRUE,group.by="final.clusters",repel=TRUE)
dev.off()
DimPlot(final.pruned.harmony.go, label=TRUE,group.by="final.clusters",repel=TRUE)

svg("34_final_umap_by_cluster_v3_pre_cluster_27_split.svg")
DimPlot(final.pruned.harmony.go, label=TRUE,group.by="final.clusters",repel=TRUE) + NoLegend()
dev.off()
DimPlot(final.pruned.harmony.go, label=TRUE,group.by="final.clusters",repel=TRUE) + NoLegend()

# split 27 -> 27 & 30

i <- 21
temp2 <- final.pruned.harmony.go[,(final.pruned.harmony.go@meta.data$final.clusters==to.do[i])]
temp3 <- DimPlot(temp2, label=FALSE,group.by="final.clusters")
temp4 <- temp3$data[,1]
temp5 <- temp3$data[,2]
X <- as.matrix(cbind(temp4,temp5))
dimnames(X)[[1]] <- colnames(temp2)
dimnames(X)[[2]] <- c("UMAP_1","UMAP_2")
p1 <- ggplot(data = as.data.frame(X), aes(x = UMAP_1, y = UMAP_2)) + geom_point() + geom_vline(xintercept=c(4.5),linetype='dashed') + geom_hline(yintercept=c(4.5),linetype='dashed') + ggtitle(paste(c("Cluster ",to.do[i],": Supervised Revisit"),sep="",collapse=""))
p1
temp10 <- X
temp10a <- temp10[,1]>(4.5)
temp10b <- temp10[,2]<(4.5)
temp11 <- temp10a+temp10b
temp12 <- names(temp11[temp11>0])
temp13 <- as.character(final.pruned.harmony.go@meta.data$final.clusters)
temp13[dimnames(final.pruned.harmony.go@meta.data)[[1]]%in%temp12] <- c(rep("30",length(temp12)))
final.pruned.harmony.go@meta.data$final.clusters <- factor(temp13,levels=as.character(c("0","1","5","8","9","10","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","29","30")))

svg("34_final_umap_by_cluster_v1_post_cluster_27_split.svg")
DimPlot(final.pruned.harmony.go, label=FALSE,group.by="final.clusters")
dev.off()
DimPlot(final.pruned.harmony.go, label=FALSE,group.by="final.clusters")

svg("34_final_umap_by_cluster_v2_post_cluster_27_split.svg")
DimPlot(final.pruned.harmony.go, label=TRUE,group.by="final.clusters",repel=TRUE)
dev.off()
DimPlot(final.pruned.harmony.go, label=TRUE,group.by="final.clusters",repel=TRUE)

svg("34_final_umap_by_cluster_v3_post_cluster_27_split.svg")
DimPlot(final.pruned.harmony.go, label=TRUE,group.by="final.clusters",repel=TRUE) + NoLegend()
dev.off()
DimPlot(final.pruned.harmony.go, label=TRUE,group.by="final.clusters",repel=TRUE) + NoLegend()

svg("35_final_umap_by_sample.svg")
DimPlot(final.pruned.harmony.go, label=FALSE, group.by="orig.ident")
dev.off()
DimPlot(final.pruned.harmony.go, label=FALSE, group.by="orig.ident")

svg("36_final_umap_split_by_sample.svg",width=20)
DimPlot(final.pruned.harmony.go, reduction = "umap", split.by = "orig.ident",ncol=4)
dev.off()
DimPlot(final.pruned.harmony.go, reduction = "umap", split.by = "orig.ident",ncol=4)

level.order <- as.character(names(rev(sort(table(final.pruned.harmony.go@meta.data$final.clusters)))))
temp <- table(final.pruned.harmony.go@meta.data$final.clusters)
temp <- as.data.frame(cbind(names(temp),as.numeric(as.character(temp))))
dimnames(temp)[[2]] <- c("Cluster","nCells")
temp$Cluster <- factor(temp$Cluster,levels=level.order)
temp$nCells <- as.numeric(as.character(temp$nCells))
svg("37_final_barplot_cluster_nCells.svg",width=20)
print(ggplot(temp, aes(x = factor(Cluster, level=level.order), y = nCells, fill = Cluster)) +  geom_bar(stat = "identity") + xlab("Cluster")) + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
print(ggplot(temp, aes(x = factor(Cluster, level=level.order), y = nCells, fill = Cluster)) +  geom_bar(stat = "identity") + xlab("Cluster")) + theme(axis.text.x = element_text(angle = 45, hjust=1))

level.order <- as.character(names(rev(sort(table(final.pruned.harmony.go@meta.data$final.clusters)))))
temp <- final.pruned.harmony.go@meta.data
temp <- as.data.frame(temp[,c(1,52)])
temptemp <- apply(temp,1,paste,collapse="_",sep="_")
temptemptemp <- table(temptemp)
temp <- cbind(names(temptemptemp),as.numeric(temptemptemp))
temptemp <- str_split(as.character(temp[,1]),"_")
temptemp <- as.data.frame(temptemp)
temptemp <- t(temptemp)
temp <- cbind(temp,temptemp)
temptemp <- temp[,c(3:4)]
temptemp <- apply(temptemp,1,paste,collapse="_",sep="_")
temp <- cbind(temp,temptemp)
rownames(temp) <- temp[,1]
temp[,1] <- temp[,6]     
temp <- temp[,c(1,5,2)]
dimnames(temp)[[2]] <- c("Sample","Cluster","nCells")
temp <- as.data.frame(temp)
temp$Cluster <- factor(as.character(temp$Cluster),levels=level.order)
temp$nCells <- as.numeric(as.character(temp$nCells))
svg("38_final_barplot_by_cluster_vs_sample_breakdown.svg",width=10)
print(ggplot(temp, aes(x = factor(Cluster, level=level.order), y = nCells, fill = Sample)) + geom_bar(stat = "identity",position="stack") + xlab('Cluster') + ylab('nCells') + theme(axis.text.x = element_text(angle = 45, hjust=1)))
dev.off()
print(ggplot(temp, aes(x = factor(Cluster, level=level.order), y = nCells, fill = Sample)) + geom_bar(stat = "identity",position="stack") + xlab('Cluster') + ylab('nCells') + theme(axis.text.x = element_text(angle = 45, hjust=1)))

svg("39_final_barplot_by_cluster_vs_sample_proportion.svg",width=10)
print(ggplot(temp, aes(x = factor(Cluster, level=level.order), y = nCells, fill = Sample)) + geom_bar(stat = "identity",position="fill") + xlab('Cluster') + ylab('Proportion') + theme(axis.text.x = element_text(angle = 45, hjust=1)))
dev.off()
print(ggplot(temp, aes(x = factor(Cluster, level=level.order), y = nCells, fill = Sample)) + geom_bar(stat = "identity",position="fill") + xlab('Cluster') + ylab('Proportion') + theme(axis.text.x = element_text(angle = 45, hjust=1)))




# perform automated annotation
# https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
# https://seurat.nygenome.org/

available_data <- AvailableData()
available_data[grep("Azimuth", available_data[, 3]), 1:3]

InstallData("humancortexref.SeuratData")
library("humancortexref.SeuratData")

final.annotated.harmony.go <- final.pruned.harmony.go
final.annotated.harmony.go.joinedlayers <- JoinLayers(final.pruned.harmony.go) #Added so that the next line can be processed
final.annotated.harmony.go <- RunAzimuth(final.annotated.harmony.go.joinedlayers, reference = "humancortexref")

temp.cluster <- as.character(final.annotated.harmony.go@meta.data$final.clusters)
temp.cohort <- final.annotated.harmony.go@meta.data$orig.ident
temp.class <- final.annotated.harmony.go@meta.data$predicted.class
temp.subclass <- final.annotated.harmony.go@meta.data$predicted.subclass
final.preds.out.1 <- cbind(colnames(final.annotated.harmony.go),temp.cohort,temp.cluster,temp.class,temp.subclass)
final.preds.out.1 <- as.data.frame(final.preds.out.1)
wb <- createWorkbook("automated_annotations")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", final.preds.out.1, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "40_post_automated_annotations_using_humancortexref.xlsx", overwrite = TRUE)

p1 <- DimPlot(final.annotated.harmony.go, group.by = "predicted.class", label = FALSE, repel = TRUE)
svg("41_post_annotation_umap_with_humancortexref_class_overlaid.svg")
print(p1)
dev.off()
p1

p2 <- DimPlot(final.annotated.harmony.go, group.by = "predicted.subclass", label = TRUE, repel = TRUE)
svg("42_post_annotation_umap_with_humancortexref_subclass_overlaid.svg")
print(p2)
dev.off()
p2

# cluster plots

temp.cx <- cbind(temp.class,temp.subclass)
temp.cx <- apply(temp.cx,1,paste,sep="_",collapse="_")
temp.cx <- str_replace_all(temp.cx," ","_")
final.preds.out.1 <- cbind(final.preds.out.1,temp.cx)

max.pred <- NULL
temp <- sort(unique(final.preds.out.1$temp.cluster))
for (i in 1:length(temp)) {
  svg(paste(c("43_post_annotation_barplot_summary_for_Cluster_",temp[i],".svg"),collapse="",sep=""))
  par(mfrow=c(1,1),mar=c(10, 4, 5, 2))
  barplot(table(final.preds.out.1[final.preds.out.1$temp.cluster==temp[i],6])/sum(table(final.preds.out.1[final.preds.out.1$temp.cluster==temp[i],6])),ylab="Membership Proportion",las = 2,space=0,cex.names = 0.7,main="Azimuth::Cortex") / abline(h=0)
  dev.off()
  temptemp <- table(final.preds.out.1[final.preds.out.1$temp.cluster==temp[i],6])/sum(table(final.preds.out.1[final.preds.out.1$temp.cluster==temp[i],6]))
  temptemp <- names(rev(sort(temptemp))[1])
  max.pred <- c(max.pred,temptemp)
}
names(max.pred) <- temp
write.table(max.pred,"final.annotation.per.cluster.azimuth.cortex.txt",sep="\t")

# cross cluster plot

temp7 <- names(table(as.character(final.preds.out.1$temp.cx)))
sumarize.cluster.types <- NULL
for (i in 1:length(temp)) {
  print(temp[i])
  temp1 <- final.preds.out.1[final.preds.out.1$temp.cluster==temp[i],]
  temp2 <- table(temp1$temp.cx)/sum(table(temp1$temp.cx))
  print(length(temp2))
  if(length(temp2)<length(temp7)) {
    temp3 <- setdiff(as.character(temp7),as.character(names(temp2)))
    temp4 <- rep(0,length(temp3))
    names(temp4) <- temp3
    temp4 <- c(temp2,temp4)
    temp4 <- temp4[temp7]
  } else {
    temp4 <- temp2[temp7]
  }
  sumarize.cluster.types <- rbind(sumarize.cluster.types,as.numeric(temp4))
}
dimnames(sumarize.cluster.types)[[1]] <- temp
dimnames(sumarize.cluster.types)[[2]] <- temp7
ppp <- t(as.matrix(sumarize.cluster.types))
ppp <- ifelse(ppp<0.01,0,ppp)
svg("45_post_annotation_final_heatmap.svg",width=20)
print(pheatmap::pheatmap(ppp,scale="none",fontsize_row=8,display_numbers=TRUE,number_color="white"))
dev.off()
print(pheatmap::pheatmap(ppp,scale="none",fontsize_row=8,display_numbers=TRUE,number_color="white"))

temp1 <- as.character(final.annotated.harmony.go@meta.data$final.clusters)
temp1 <- ifelse(temp1=="0","Glutamatergic_L2/3_IT",temp1)
temp1 <- ifelse(temp1=="1","Glutamatergic_L2/3_IT",temp1)
temp1 <- ifelse(temp1=="5","Glutamatergic_L2/3_IT",temp1)
temp1 <- ifelse(temp1=="8","Non-Neuronal_Astro",temp1)
temp1 <- ifelse(temp1=="9","Non-Neuronal_Oligo",temp1)
temp1 <- ifelse(temp1=="10","Glutamatergic_L5_IT",temp1)
temp1 <- ifelse(temp1=="12","Glutamatergic_L5_IT",temp1)
temp1 <- ifelse(temp1=="13","Non-Neuronal_Micro-PVM",temp1)
temp1 <- ifelse(temp1=="14","Glutamatergic_L5_ET",temp1)
temp1 <- ifelse(temp1=="15","Non-Neuronal_OPC",temp1)
temp1 <- ifelse(temp1=="17","GABAergic_Meis2",temp1)
temp1 <- ifelse(temp1=="18","Glutamatergic_L2/3_IT",temp1)
temp1 <- ifelse(temp1=="19","GABAergic_Sst",temp1)
temp1 <- ifelse(temp1=="20","Non-Neuronal_Astro",temp1)
temp1 <- ifelse(temp1=="21","GABAergic_Lamp5",temp1)
temp1 <- ifelse(temp1=="22","GABAergic_Vip",temp1)
temp1 <- ifelse(temp1=="23","Glutamatergic_L6b",temp1)
temp1 <- ifelse(temp1=="24","GABAergic_Lamp5",temp1)
temp1 <- ifelse(temp1=="25","Glutamatergic_L2/3_IT",temp1)
temp1 <- ifelse(temp1=="26","Glutamatergic_L5/6_NP",temp1)
temp1 <- ifelse(temp1=="27","Non-Neuronal_VLMC",temp1)
temp1 <- ifelse(temp1=="28","Non-Neuronal_Peri",temp1)
temp1 <- ifelse(temp1=="29","Non-Neuronal_Oligo",temp1)
temp1 <- ifelse(temp1=="30","Non-Neuronal_Endo",temp1)
final.annotated.harmony.go@meta.data$final.annotations <- as.factor(temp1)






saveRDS(object = final.curated.harmony.go, file = "final.curated.harmony.go.rds")
saveRDS(object = final.pruned.harmony.go, file = "final.pruned.harmony.go.rds")
saveRDS(object = final.annotated.harmony.go, file = "final.annotated.harmony.go.rds")






# Conserved testing 

final.clusters <- table(final.annotated.harmony.go@meta.data$final.clusters)
final.clusters <- names(final.clusters)
DefaultAssay(final.annotated.harmony.go) <- "RNA"
Idents(final.annotated.harmony.go) <- factor(final.annotated.harmony.go@meta.data$final.clusters,levels=as.character(sort(as.numeric(final.clusters))))
conserved.markers.0 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "0", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.1 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "1", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.5 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "5", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.8 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "8", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.9 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "9", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.10 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "10", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.12 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "12", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.13 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "13", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.14 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "14", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.15 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "15", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.17 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "17", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.18 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "18", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.19 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "19", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.20 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "20", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.21 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "21", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.22 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "22", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.23 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "23", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.24 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "24", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.25 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "25", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.26 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "26", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.27 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "27", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.28 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "28", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.29 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "29", grouping.var="orig.ident", only.pos = TRUE)
conserved.markers.30 <- FindConservedMarkers(final.annotated.harmony.go, ident.1 = "30", grouping.var="orig.ident", only.pos = TRUE)

wb <- createWorkbook("conserved.markers.0")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.0, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_0.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.1")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.1, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_1.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.5")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.5, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_5.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.8")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.8, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_8.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.9")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.9, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_9.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.10")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.10, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_10.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.12")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.12, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_12.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.13")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.13, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_13.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.14")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.14, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_14.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.15")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.15, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_15.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.17")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.17, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_17.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.18")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.18, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_18.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.19")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.19, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_19.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.20")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.20, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_20.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.21")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.21, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_21.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.22")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.22, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_22.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.23")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.23, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_23.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.24")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.24, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_24.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.25")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.25, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_25.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.26")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.26, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_26.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.27")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.27, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_27.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.28")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.28, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_28.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.29")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.29, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_29.xlsx", overwrite = TRUE)

wb <- createWorkbook("conserved.markers.30")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", conserved.markers.30, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "46_conserved_marker_testing_results_Cluster_30.xlsx", overwrite = TRUE)

markers.to.plot <- c(rownames(conserved.markers.0)[1],
                     rownames(conserved.markers.1)[1],
                     rownames(conserved.markers.5)[1],
                     rownames(conserved.markers.8)[1],
                     rownames(conserved.markers.9)[1],
                     rownames(conserved.markers.10)[1],
                     rownames(conserved.markers.12)[1],
                     rownames(conserved.markers.13)[1],
                     rownames(conserved.markers.14)[1],
                     rownames(conserved.markers.15)[1],
                     rownames(conserved.markers.17)[1],
                     rownames(conserved.markers.18)[1],
                     rownames(conserved.markers.19)[1],
                     rownames(conserved.markers.20)[1],
                     rownames(conserved.markers.21)[1],
                     rownames(conserved.markers.22)[1],
                     rownames(conserved.markers.23)[1],
                     rownames(conserved.markers.24)[1],
                     rownames(conserved.markers.25)[1],
                     rownames(conserved.markers.26)[1],
                     rownames(conserved.markers.27)[1],
                     rownames(conserved.markers.28)[1],
                     rownames(conserved.markers.29)[1],
                     rownames(conserved.markers.30)[1])
markers.to.plot <- unique(markers.to.plot)

svg("47_post_conserved_marker_testing_dot_plot.svg",width=12)
DotPlot(final.annotated.harmony.go, features = unique(markers.to.plot), dot.scale = 8) + RotatedAxis()
dev.off()
DotPlot(final.annotated.harmony.go, features = unique(markers.to.plot), dot.scale = 8) + RotatedAxis()

p <- VlnPlot(final.annotated.harmony.go,markers.to.plot,stack=TRUE,sort=FALSE,flip=TRUE,fill.by="ident")
p <- p + xlab("Cluster")
svg("48_post_conserved_marker_testing_violin_plot.svg",width=20)
print(p)
dev.off()
print(p)

# Non-Conserved testing 

non.conserved.markers.0 <- FindMarkers(final.annotated.harmony.go, ident.1 = "0", only.pos = TRUE)
non.conserved.markers.1 <- FindMarkers(final.annotated.harmony.go, ident.1 = "1", only.pos = TRUE)
non.conserved.markers.5 <- FindMarkers(final.annotated.harmony.go, ident.1 = "5", only.pos = TRUE)
non.conserved.markers.8 <- FindMarkers(final.annotated.harmony.go, ident.1 = "8", only.pos = TRUE)
non.conserved.markers.9 <- FindMarkers(final.annotated.harmony.go, ident.1 = "9", only.pos = TRUE)
non.conserved.markers.10 <- FindMarkers(final.annotated.harmony.go, ident.1 = "10", only.pos = TRUE)
non.conserved.markers.12 <- FindMarkers(final.annotated.harmony.go, ident.1 = "12", only.pos = TRUE)
non.conserved.markers.13 <- FindMarkers(final.annotated.harmony.go, ident.1 = "13", only.pos = TRUE)
non.conserved.markers.14 <- FindMarkers(final.annotated.harmony.go, ident.1 = "14", only.pos = TRUE)
non.conserved.markers.15 <- FindMarkers(final.annotated.harmony.go, ident.1 = "15", only.pos = TRUE)
non.conserved.markers.17 <- FindMarkers(final.annotated.harmony.go, ident.1 = "17", only.pos = TRUE)
non.conserved.markers.18 <- FindMarkers(final.annotated.harmony.go, ident.1 = "18", only.pos = TRUE)
non.conserved.markers.19 <- FindMarkers(final.annotated.harmony.go, ident.1 = "19", only.pos = TRUE)
non.conserved.markers.20 <- FindMarkers(final.annotated.harmony.go, ident.1 = "20", only.pos = TRUE)
non.conserved.markers.21 <- FindMarkers(final.annotated.harmony.go, ident.1 = "21", only.pos = TRUE)
non.conserved.markers.22 <- FindMarkers(final.annotated.harmony.go, ident.1 = "22", only.pos = TRUE)
non.conserved.markers.23 <- FindMarkers(final.annotated.harmony.go, ident.1 = "23", only.pos = TRUE)
non.conserved.markers.24 <- FindMarkers(final.annotated.harmony.go, ident.1 = "24", only.pos = TRUE)
non.conserved.markers.25 <- FindMarkers(final.annotated.harmony.go, ident.1 = "25", only.pos = TRUE)
non.conserved.markers.26 <- FindMarkers(final.annotated.harmony.go, ident.1 = "26", only.pos = TRUE)
non.conserved.markers.27 <- FindMarkers(final.annotated.harmony.go, ident.1 = "27", only.pos = TRUE)
non.conserved.markers.28 <- FindMarkers(final.annotated.harmony.go, ident.1 = "28", only.pos = TRUE)
non.conserved.markers.29 <- FindMarkers(final.annotated.harmony.go, ident.1 = "29", only.pos = TRUE)
non.conserved.markers.30 <- FindMarkers(final.annotated.harmony.go, ident.1 = "30", only.pos = TRUE)

wb <- createWorkbook("non.conserved.markers.0")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.0, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_0.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.1")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.1, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_1.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.5")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.5, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_5.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.8")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.8, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_8.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.9")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.9, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_9.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.10")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.10, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_10.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.12")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.12, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_12.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.13")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.13, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_13.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.14")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.14, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_14.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.15")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.15, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_15.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.17")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.17, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_17.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.18")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.18, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_18.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.19")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.19, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_19.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.20")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.20, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_20.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.21")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.21, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_21.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.22")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.22, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_22.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.23")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.23, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_23.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.24")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.24, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_24.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.25")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.25, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_25.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.26")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.26, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_26.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.27")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.27, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_27.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.28")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.28, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_28.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.29")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.29, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_29.xlsx", overwrite = TRUE)

wb <- createWorkbook("non.conserved.markers.30")
addWorksheet(wb, "results")
writeData(wb, sheet = "results", non.conserved.markers.30, colNames = TRUE,rowNames = TRUE)
saveWorkbook(wb, "49_non_conserved_marker_testing_results_Cluster_30.xlsx", overwrite = TRUE)

markers.to.plot <- c(rownames(non.conserved.markers.0)[1],
                     rownames(non.conserved.markers.1)[1],
                     rownames(non.conserved.markers.5)[1],
                     rownames(non.conserved.markers.8)[1],
                     rownames(non.conserved.markers.9)[1],
                     rownames(non.conserved.markers.10)[1],
                     rownames(non.conserved.markers.12)[1],
                     rownames(non.conserved.markers.13)[1],
                     rownames(non.conserved.markers.14)[1],
                     rownames(non.conserved.markers.15)[1],
                     rownames(non.conserved.markers.17)[1],
                     rownames(non.conserved.markers.18)[1],
                     rownames(non.conserved.markers.19)[1],
                     rownames(non.conserved.markers.20)[1],
                     rownames(non.conserved.markers.21)[1],
                     rownames(non.conserved.markers.22)[1],
                     rownames(non.conserved.markers.23)[1],
                     rownames(non.conserved.markers.24)[1],
                     rownames(non.conserved.markers.25)[1],
                     rownames(non.conserved.markers.26)[1],
                     rownames(non.conserved.markers.27)[1],
                     rownames(non.conserved.markers.28)[1],
                     rownames(non.conserved.markers.29)[1],
                     rownames(non.conserved.markers.30)[1])
markers.to.plot <- unique(markers.to.plot)

svg("50_post_non_conserved_marker_testing_dot_plot.svg",width=12)
DotPlot(final.annotated.harmony.go, features = unique(markers.to.plot), dot.scale = 8) + RotatedAxis()
dev.off()
DotPlot(final.annotated.harmony.go, features = unique(markers.to.plot), dot.scale = 8) + RotatedAxis()

p <- VlnPlot(final.annotated.harmony.go,markers.to.plot,stack=TRUE,sort=FALSE,flip=TRUE,fill.by="ident")
p <- p + xlab("Cluster")
svg("51_post_non_conserved_marker_testing_violin_plot.svg",width=20)
print(p)
dev.off()
print(p)

#------------------------------------------------------------------------------
# Generate additional metric plots
#------------------------------------------------------------------------------

Idents(final.annotated.harmony.go) <- final.annotated.harmony.go@meta.data$final.clusters


# add meta column for M and F coding
code.m.f <- final.annotated.harmony.go@meta.data$orig.ident
code.m.f  <- ifelse(code.m.f=="DG_F1_VV","F",code.m.f)
code.m.f  <- ifelse(code.m.f=="DG_F2_VV","F",code.m.f)
code.m.f  <- ifelse(code.m.f=="LL_F1_DG","F",code.m.f)
code.m.f  <- ifelse(code.m.f=="LL_F2_DG","F",code.m.f)
code.m.f <- ifelse(code.m.f=="F",code.m.f,"M")
final.annotated.harmony.go@meta.data$code.m.f <- as.factor(code.m.f)
Idents(final.annotated.harmony.go) <- final.annotated.harmony.go@meta.data$code.m.f

# generate umap plot overlating M vs F
DimPlot(final.annotated.harmony.go, label=FALSE,group.by="code.m.f",cols=c("lightgrey","black"),pt.size=1)

# generate split umap plot overlating M vs F
DimPlot(final.annotated.harmony.go, label=FALSE,split.by="code.m.f",cols=c("lightgrey","black"),pt.size=1)

# generate bar plot describing number cells per M and F class
level.order <- c("F","M")
temp <- table(final.annotated.harmony.go@meta.data$code.m.f)
temp <- as.data.frame(cbind(names(temp),as.numeric(as.character(temp))))
dimnames(temp)[[2]] <- c("Gender","nCells")
temp$Gender <- factor(temp$Gender,levels=level.order)
temp$nCells <- as.numeric(as.character(temp$nCells))
print(ggplot(temp, aes(x = factor(Gender, level=level.order), y = nCells, fill = Gender)) +  geom_bar(stat = "identity") + xlab("Cluster") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_fill_manual(values=c("lightgrey","black")))

# output bar plot stats
table(code.m.f)

# generate stacked bar plot describing number cells per cluster per M and F class
level.order <- c("F","M")
temp <- data.frame(Class=as.character(unlist(final.annotated.harmony.go@meta.data$code.m.f)),Cluster=as.character(unlist(final.annotated.harmony.go@meta.data$final.clusters)))
temptemp <- apply(temp,1,paste,collapse="_",sep="_")
temptemptemp <- table(temptemp)
temp <- cbind(names(temptemptemp),as.numeric(temptemptemp))
temptemp <- str_split(as.character(temp[,1]),"_")
temptemp <- as.data.frame(temptemp)
temptemp <- t(temptemp)
temp <- cbind(temp,temptemp)
rownames(temp) <- temp[,1]
temp <- temp[,-c(1)]
temp <- data.frame(Gender=temp[,2],Cluster=temp[,3],nCells=as.numeric(temp[,1]))
temp$Cluster <- factor(as.character(temp$Cluster),levels=as.character(c(0,1,5,8,9,10,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30)))
temp$Gender <- factor(as.character(temp$Gender),levels=level.order)
print(ggplot(temp, aes(x = factor(Gender, level=level.order), y = nCells, fill = Cluster)) + geom_bar(stat = "identity",position="stack") + xlab('Gender') + ylab('nCells') + theme(axis.text.x = element_text(angle = 45, hjust=1)))

# generate proportion bar plot describing number cells per cluster per M and F class
print(ggplot(temp, aes(x = factor(Gender, level=level.order), y = nCells, fill = Cluster)) + geom_bar(stat = "identity",position="fill") + xlab('Gender') + ylab('nCells') + theme(axis.text.x = element_text(angle = 45, hjust=1)))

# output bar plot stats
temp.f <- sum(temp[,3][temp[,1]=="F"])
temp.m <- sum(temp[,3][temp[,1]=="M"])
temptemp <- ifelse(temp[,1]=="F",temp[,3]/temp.f,temp[,3]/temp.m)
temp <- cbind(temp,temptemp)
dimnames(temp)[[2]][4] <- "Proportion"
print(temp)

# umap with MT overlay

temp <- DimPlot(final.annotated.harmony.go)
temp <- temp$data
percent.mt <- final.annotated.harmony.go@meta.data$percent.mt
names(percent.mt) <- rownames(final.annotated.harmony.go@meta.data)
percent.mt <- percent.mt[rownames(temp)]
percent.nfeatures <- final.annotated.harmony.go@meta.data$nFeature_RNA
names(percent.nfeatures) <- rownames(final.annotated.harmony.go@meta.data)
percent.nfeatures <- percent.nfeatures[rownames(temp)]
percent.ncount <- final.annotated.harmony.go@meta.data$nCount_RNA
names(percent.ncount) <- rownames(final.annotated.harmony.go@meta.data)
percent.ncount <- percent.ncount[rownames(temp)]
doublet.score <- final.annotated.harmony.go@meta.data$DoubletScore
names(doublet.score) <- rownames(final.annotated.harmony.go@meta.data)
doublet.score <- doublet.score[rownames(temp)]
temp <- cbind(temp,percent.mt,percent.nfeatures,percent.ncount,doublet.score)
ggplot(temp, aes(x=umap_1, y=umap_2,colour=percent.mt)) + geom_point(size=1) + scale_colour_gradient(low="black",high="gray")

# umap with percent nFeatures overlay
ggplot(temp, aes(x=umap_1, y=umap_2,colour=percent.nfeatures)) + geom_point(size=1) + scale_colour_gradient(low="black",high="gray")

# umap with percent nCount overlay
ggplot(temp, aes(x=umap_1, y=umap_2,colour=percent.ncount)) + geom_point(size=1) + scale_colour_gradient(low="black",high="gray")

# umap with Doublet Score overlay
ggplot(temp, aes(x=umap_1, y=umap_2,colour=doublet.score)) + geom_point(size=1) + scale_colour_gradient(low="black",high="gray")

#------------------------------------------------------------------------------
# trajectory fitting LL & VV 
#------------------------------------------------------------------------------

# source required libraries

library(remotes)
remotes::install_github("satijalab/seurat-wrappers")
library(monocle3)
library(scales)

# load working object

final.annotated.harmony.go <- readRDS("final.annotated.harmony.go.rds")
Idents(final.annotated.harmony.go) <- final.annotated.harmony.go@meta.data$final.clusters

# add meta column for VV and LL coding

code.vv.ll <- final.annotated.harmony.go@meta.data$orig.ident
code.vv.ll <- ifelse(code.vv.ll=="DG_F1_VV","VV",code.vv.ll)
code.vv.ll <- ifelse(code.vv.ll=="DG_F2_VV","VV",code.vv.ll)
code.vv.ll <- ifelse(code.vv.ll=="DG_M1_VV","VV",code.vv.ll)
code.vv.ll <- ifelse(code.vv.ll=="DG_M2_VV","VV",code.vv.ll)
code.vv.ll <- ifelse(code.vv.ll=="VV",code.vv.ll,"LL")
final.annotated.harmony.go@meta.data$code.vv.ll <- as.factor(code.vv.ll)
Idents(final.annotated.harmony.go) <- final.annotated.harmony.go@meta.data$code.vv.ll

# prep for fit

prep.4.fit  <- final.annotated.harmony.go[,final.annotated.harmony.go@meta.data$final.clusters%in%c(0,1,5,18)]
DefaultAssay(prep.4.fit) <- "RNA"
monocle_object <- SeuratWrappers::as.cell_data_set(prep.4.fit)
monocle_object <- preprocess_cds(monocle_object)
monocle_object <- cluster_cells(cds = monocle_object, reduction_method = "UMAP")

# perform fit

monocle_object <- learn_graph(monocle_object)
monocle_object <- estimate_size_factors(monocle_object)
monocle_object@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(final.annotated.harmony.go[["RNA"]])
monocle_object.cx <- monocle_object 

# plot fit

ccc <- hue_pal()(24)[c(1,2,3,12)]
cx <- plot_cells(monocle_object.cx,color_cells_by = "final.clusters", label_cell_groups=FALSE,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=1,trajectory_graph_segment_size=2,trajectory_graph_color = "black") + scale_colour_manual(values=ccc)
cx

# evaluate root node choice for pseudotime testing

goi <- c("Camk4","Ntng1")
g <- plot_cells(monocle_object.cx,genes=goi,label_cell_groups=FALSE,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=1,trajectory_graph_segment_size=2,trajectory_graph_color = "black")
g

# select root

monocle_object.cx <- order_cells(monocle_object.cx,reduction_method = "UMAP")
cx.fit <- plot_cells(monocle_object.cx,color_cells_by = "pseudotime", label_cell_groups=FALSE,show_trajectory_graph=TRUE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=1,trajectory_graph_segment_size=2,trajectory_graph_color = "black",graph_label_size=5)
cx.fit

# test for and identify genes that explain pseudotime

cx.fit.genes <- graph_test(monocle_object.cx, neighbor_graph="principal_graph")
cx.fit.genes.filtered <- row.names(subset(cx.fit.genes, q_value < 0.05))
write.table(cx.fit.genes.filtered,"cx.fit.genes.filtered.txt",sep="\t")

# plot top gene by q_value then p_value then morans_test_statistic

plot_genes_in_pseudotime(monocle_object.cx[rowData(monocle_object.cx)$gene_short_name %in% c("Sorcs3"),],color_cells_by="pseudotime",min_expr=0.5)

plot_cells(monocle_object.cx,genes=c("Sorcs3"),label_cell_groups=FALSE,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=1,trajectory_graph_segment_size=2,trajectory_graph_color = "black",graph_label_size=5)

# test for and identify depleted vs enriched modules of genes that were identified to explain pseduotime

gene_module_df.cx <- find_gene_modules(monocle_object.cx[cx.fit.genes.filtered,], resolution=c(10^seq(-6,-1)))
write.table(gene_module_df.cx,"gene_module_df.cx.txt",sep="\t")
cell_group_df.cx <- tibble::tibble(cell=row.names(colData(monocle_object.cx)),cell_group=colData(monocle_object.cx)$final.clusters)

# plot module heat map

agg_mat.cx <- aggregate_gene_expression(monocle_object.cx, gene_module_df.cx, cell_group_df.cx)
row.names(agg_mat.cx) <- stringr::str_c("Module ", row.names(agg_mat.cx))
pheatmap::pheatmap(agg_mat.cx,scale="column", clustering_method="ward.D2")

# plot top traj gene for cluster 0

plot_genes_in_pseudotime(monocle_object.cx[rowData(monocle_object.cx)$gene_short_name %in% c("Ahcyl2"),],color_cells_by="pseudotime",min_expr=0.5)

plot_cells(monocle_object.cx,genes=c("Ahcyl2"),label_cell_groups=FALSE,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=1,trajectory_graph_segment_size=2,trajectory_graph_color = "black",graph_label_size=5)

# plot top traj gene for cluster 1

plot_genes_in_pseudotime(monocle_object.cx[rowData(monocle_object.cx)$gene_short_name %in% c("Cadm2"),],color_cells_by="pseudotime",min_expr=0.5)

plot_cells(monocle_object.cx,genes=c("Cadm2"),label_cell_groups=FALSE,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=1,trajectory_graph_segment_size=2,trajectory_graph_color = "black",graph_label_size=5)

# plot top traj gene for cluster 5

plot_genes_in_pseudotime(monocle_object.cx[rowData(monocle_object.cx)$gene_short_name %in% c("Olfm1"),],color_cells_by="pseudotime",min_expr=0.5)

plot_cells(monocle_object.cx,genes=c("Olfm1"),label_cell_groups=FALSE,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=1,trajectory_graph_segment_size=2,trajectory_graph_color = "black",graph_label_size=5)

# plot top traj gene for cluster 18

plot_genes_in_pseudotime(monocle_object.cx[rowData(monocle_object.cx)$gene_short_name %in% c("Sorcs3"),],color_cells_by="pseudotime",min_expr=0.5)

plot_cells(monocle_object.cx,genes=c("Sorcs3"),label_cell_groups=FALSE,show_trajectory_graph=FALSE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=1,trajectory_graph_segment_size=2,trajectory_graph_color = "black",graph_label_size=5)

