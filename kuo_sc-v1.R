setwd('C://Users/lvlah/linux/ac_lab/kuo-colab/sc/')
library(PISCES)
library(Seurat)
library(DESeq2)
library(viper)
library(pheatmap)
library(RColorBrewer)
# singleR
library(SingleR)
library(celldex)
mouse.ref <- MouseRNAseqData()

DESeqGES <- function(test.counts, ref.counts) {
  # make combined matrix
  comb.mat <- cbind(test.counts, ref.counts)
  # make meta data
  meta.df <- data.frame('Group' = c(rep('Test', ncol(test.counts)), rep('Ref', ncol(ref.counts))))
  # make DESeq object
  deseq.obj <- DESeqDataSetFromMatrix(countData = comb.mat, colData = meta.df, design = ~Group)
  deseq.obj$Group <- relevel(deseq.obj$Group, ref = 'Ref')
  deseq.obj <- DESeq(deseq.obj)
  # extract GES
  ges.df <- as.data.frame(results(object = deseq.obj, contrast = c('Group', 'Test', 'Ref'),
                                  cooksCutoff = FALSE, independentFiltering = FALSE))
  ges.vals <- qnorm(p = ges.df$pvalue / 2, lower.tail = FALSE) * sign(ges.df$log2FoldChange)
  names(ges.vals) <- rownames(ges.df)
  # return final vector
  return(ges.vals)
}
pact.col <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))

### initial K0 analysis
###############
k0.counts <- Read10X('KO/filtered_feature_bc_matrix/')
k0.seur <- CreateSeuratObject(counts = k0.counts, project = 'k0')
## filter; normalize; SCT
mt.features <- intersect(mt.genes$mur.symb, rownames(k0.seur))
k0.seur[["percent.mt"]] <- PercentageFeatureSet(object = k0.seur, features = mt.features)
qc.plot <- QCPlots(k0.seur)
qc.plot <- annotate_figure(qc.plot, top = text_grob("KO QC Metrics", size = 18))
jpeg('C://Users/lvlah/linux/ac_lab/kuo-colab/paper_components/fig-s1_b.png', width = 8, height = 4, units = "in", res = 2000)
print(qc.plot)
dev.off()
## subset
k0.seur <- subset(k0.seur, subset = nCount_RNA > 1000 & nFeature_RNA > 1250 & percent.mt < 25)
k0.seur <- SCTransform(k0.seur, vars.to.regress = 'percent.mt', verbose = FALSE)
## make SingleR labels
k0.singleR <- SingleR(k0.seur@assays$RNA@data, ref = mouse.ref, labels = mouse.ref$label.main, de.method="wilcox")
k0.seur@misc[['SingleR']] <- k0.singleR
k0.seur <- k0.seur[, which(k0.singleR$labels == 'Endothelial cells')]
## gene expression cluster
k0.seur <- CorDist(k0.seur)
k0.seur <- RunPCA(k0.seur, verbose = FALSE)
k0.seur <- RunUMAP(k0.seur, dims = 1:30, verbose = FALSE)
k0.seur <- FindNeighbors(k0.seur, dims = 1:30, verbose = FALSE)
k0.seur <- FindClusters(k0.seur, verbose = FALSE)
DimPlot(k0.seur, label = TRUE) + NoLegend()
## make aracne matrix
arac.cells <- sample(colnames(k0.seur), 500)
arac.mat <- CPMTransform(as.matrix(k0.seur@assays$RNA@counts)[, arac.cells])
saveRDS(arac.mat, 'k0_arac.rds')
meta.arac <- MetaCells(as.matrix(k0.seur@assays$RNA@counts), k0.seur@assays$SCT@misc$dist.mat, 
                       num.neighbors = 5, subset = NULL)
saveRDS(meta.arac[[1]][, arac.cells], 'k0_k5-arac.rds')
## save
saveRDS(k0.seur, file = 'k0_seur.rds')
###############

### initial K0F4 analysis
###############
k0f4.counts <- Read10X('KOF4/filtered_feature_bc_matrix/')
k0f4.seur <- CreateSeuratObject(counts = k0f4.counts, project = 'k0f4')
## filter; normalize; SCT
mt.features <- intersect(mt.genes$mur.symb, rownames(k0f4.seur))
k0f4.seur[["percent.mt"]] <- PercentageFeatureSet(object = k0f4.seur, features = mt.features)
qc.plot <- QCPlots(k0f4.seur)
qc.plot <- annotate_figure(qc.plot, top = text_grob("KOF4 QC Metrics", size = 18))
jpeg('C://Users/lvlah/linux/ac_lab/kuo-colab/paper_components/fig-s1_c.png', width = 8, height = 4, units = "in", res = 2000)
print(qc.plot)
dev.off()
k0f4.seur <- subset(k0f4.seur, subset = nCount_RNA > 3000 & nFeature_RNA > 1000 & percent.mt < 12.5)
k0f4.seur <- SCTransform(k0f4.seur, vars.to.regress = 'percent.mt', verbose = FALSE)
## make SingleR labels
k0f4.singleR <- SingleR(k0f4.seur@assays$RNA@data, ref = mouse.ref, labels = mouse.ref$label.main, de.method="wilcox")
k0f4.seur@misc[['SingleR']] <- k0f4.singleR
k0f4.seur <- k0f4.seur[, which(k0f4.singleR$labels == 'Endothelial cells')]
## gene expression cluster
k0f4.seur <- CorDist(k0f4.seur)
k0f4.seur <- RunPCA(k0f4.seur, verbose = FALSE)
k0f4.seur <- RunUMAP(k0f4.seur, dims = 1:30, verbose = FALSE)
k0f4.seur <- FindNeighbors(k0f4.seur, dims = 1:30, verbose = FALSE)
k0f4.seur <- FindClusters(k0f4.seur, verbose = FALSE)
DimPlot(k0f4.seur, label = TRUE) + NoLegend()
## make aracne matrix
arac.cells <- sample(colnames(k0f4.seur), 500)
arac.mat <- CPMTransform(as.matrix(k0f4.seur@assays$RNA@counts)[, arac.cells])
saveRDS(arac.mat, 'k0f4_arac.rds')
meta.arac <- MetaCells(as.matrix(k0f4.seur@assays$RNA@counts), k0f4.seur@assays$SCT@misc$dist.mat, 
                       num.neighbors = 5, subset = NULL)
saveRDS(meta.arac[[1]][, arac.cells], 'k0f4_k5-arac.rds')
saveRDS(k0f4.seur, file = 'k0f4_seur.rds')
###############

### initial WT analysis
###############
wt.counts <- Read10X('WT/filtered_feature_bc_matrix/')
wt.seur <- CreateSeuratObject(counts = wt.counts, project = 'wt')
## filter; normalize; SCT
mt.features <- intersect(mt.genes$mur.symb, rownames(wt.seur))
wt.seur[["percent.mt"]] <- PercentageFeatureSet(object = wt.seur, features = mt.features)
qc.plot <- QCPlots(wt.seur)
qc.plot <- annotate_figure(qc.plot, top = text_grob("WT QC Metrics", size = 18))
jpeg('C://Users/lvlah/linux/ac_lab/kuo-colab/paper_components/fig-s1_a.png', width = 8, height = 4, units = "in", res = 2000)
print(qc.plot)
dev.off()
wt.seur <- subset(wt.seur, subset = nCount_RNA > 1250 & nCount_RNA < 75000 & nFeature_RNA > 1000 & percent.mt < 12.5)
wt.seur <- SCTransform(wt.seur, vars.to.regress = 'percent.mt', verbose = FALSE)
## make SingleR labels
wt.singleR <- SingleR(wt.seur@assays$RNA@data, ref = mouse.ref, labels = mouse.ref$label.main, de.method = "wilcox")
wt.seur@misc[['SingleR']] <- wt.singleR
wt.seur <- wt.seur[, which(wt.singleR$labels == 'Endothelial cells')]
## gene expression cluster
wt.seur <- CorDist(wt.seur)
wt.seur <- RunPCA(wt.seur, verbose = FALSE)
wt.seur <- RunUMAP(wt.seur, dims = 1:30, verbose = FALSE)
wt.seur <- FindNeighbors(wt.seur, dims = 1:30, verbose = FALSE)
wt.seur <- FindClusters(wt.seur, verbose = FALSE)
DimPlot(wt.seur, label = TRUE) + NoLegend()
## make aracne matrix
arac.cells <- sample(colnames(wt.seur), 500)
arac.mat <- CPMTransform(as.matrix(wt.seur@assays$RNA@counts)[, arac.cells])
saveRDS(arac.mat, 'wt_arac.rds')
meta.arac <- MetaCells(as.matrix(wt.seur@assays$RNA@counts), wt.seur@assays$SCT@misc$dist.mat, 
                       num.neighbors = 5, subset = NULL)
saveRDS(meta.arac[[1]][, arac.cells], 'wt_k5-arac.rds')
saveRDS(wt.seur, file = 'wt_seur.rds')
###############

### initial WTF4 analysis
###############
wtf4.counts <- Read10X('WTF4/filtered_feature_bc_matrix/')
wtf4.seur <- CreateSeuratObject(counts = wtf4.counts, project = 'wt')
## filter; normalize; SCT
mt.features <- intersect(mt.genes$mur.symb, rownames(wtf4.seur))
wtf4.seur[["percent.mt"]] <- PercentageFeatureSet(object = wtf4.seur, features = mt.features)
QCPlots(wtf4.seur)
wtf4.seur <- subset(wtf4.seur, subset = nCount_RNA > 1250 & nCount_RNA < 75000 & nFeature_RNA > 1250 & percent.mt < 12.5)
wtf4.seur <- SCTransform(wtf4.seur, vars.to.regress = 'percent.mt', verbose = FALSE)
## make SingleR labels
wtf4.singleR <- SingleR(wtf4.seur@assays$RNA@data, ref = mouse.ref, labels = mouse.ref$label.main, de.method = "wilcox")
wtf4.seur@misc[['SingleR']] <- wtf4.singleR
wtf4.seur <- wtf4.seur[, which(wtf4.singleR$labels == 'Endothelial cells')]
## gene expression cluster
wtf4.seur <- CorDist(wtf4.seur)
wtf4.seur <- RunPCA(wtf4.seur, verbose = FALSE)
wtf4.seur <- RunUMAP(wtf4.seur, dims = 1:30, verbose = FALSE)
wtf4.seur <- FindNeighbors(wtf4.seur, dims = 1:30, verbose = FALSE)
wtf4.seur <- FindClusters(wtf4.seur, verbose = FALSE)
DimPlot(wtf4.seur, label = TRUE) + NoLegend()
## make aracne matrix
arac.cells <- sample(colnames(wtf4.seur), 500)
arac.mat <- CPMTransform(as.matrix(wtf4.seur@assays$RNA@counts)[, arac.cells])
saveRDS(arac.mat, 'wtf4_arac.rds')
meta.arac <- MetaCells(as.matrix(wtf4.seur@assays$RNA@counts), wtf4.seur@assays$SCT@misc$dist.mat, 
                       num.neighbors = 5, subset = NULL)
saveRDS(meta.arac[[1]][, arac.cells], 'wtf4_k5-arac.rds')
saveRDS(wtf4.seur, file = 'wtf4_seur.rds')
###############

### versus WT analysis
###############
k0.seur <- readRDS('k0_seur.rds')
k0f4.seur <- readRDS('k0f4_seur.rds')
wt.seur <- readRDS('wt_seur.rds')
## load networks
k0.net <- readRDS('sc-nets/k0-net_pruned.rds')
k0.meta <- readRDS('sc-nets/k0_k5-net_pruned.rds')
k0f4.net <- readRDS('sc-nets/k0f4-net_pruned.rds')
k0f4.meta <- readRDS('sc-nets/k0f4_k5-net_pruned.rds')
## separate counts
k0.counts <- as.matrix(k0.seur@assays$SCT@counts)
k0f4.counts <- as.matrix(k0f4.seur@assays$SCT@counts)
wt.counts <- as.matrix(wt.seur@assays$SCT@counts)
## filter for shared genes
shared.genes <- intersect(intersect(rownames(k0.counts), rownames(k0f4.counts)), rownames(wt.counts))
k0.counts <- k0.counts[shared.genes,]
k0f4.counts <- k0f4.counts[shared.genes,]
wt.counts <- wt.counts[shared.genes,]
## make ges
k0.ges <- DESeqGES(k0.counts, wt.counts)
saveRDS(k0.ges, 'k0-v-wt_deseq2-ges.rds')
k0f4.ges <- DESeqGES(k0f4.counts, wt.counts)
saveRDS(k0f4.ges, 'k0f4-v-wt_deseq2-ges.rds')
## run viper
ges.mat <- cbind(as.matrix(k0.ges), as.matrix(k0f4.ges))
colnames(ges.mat) <- c('K0', 'K0F4')
vip.mat <- viper(ges.mat, k0.meta)
saveRDS(vip.mat, 'k0-k0f4-v-wt_vip.rds')
## write tables
num.mrs <- 50
k0.mrs <- sort(vip.mat[,1]); k0.mrs <- c(k0.mrs[1:num.mrs], tail(k0.mrs, num.mrs))
write.table(as.data.frame(k0.mrs), file = 'k0-v-wt_vip-mrs.csv', sep = ',', quote = FALSE)
k0f4.mrs <- sort(vip.mat[,2]); k0.mrs <- c(k0f4.mrs[1:num.mrs], tail(k0f4.mrs, num.mrs))
write.table(as.data.frame(k0f4.mrs), file = 'k0f4-v-wt_vip-mrs.csv', sep = ',', quote = FALSE)
## compare master regulators
nes.thresh <- qnorm(1 - (0.05 / nrow(vip.mat)))
k0.mrs <- vip.mat[which(abs(vip.mat[,1]) > nes.thresh), 1]
k0f4.mrs <- vip.mat[which(abs(vip.mat[,2]) > nes.thresh), 2]
## sample-level viper
wt.cpm <- CPMTransform(wt.counts)
k0.cpm <- CPMTransform(k0.counts)
k0f4.cpm <- CPMTransform(k0f4.counts)
wt.mean <- rowMeans(wt.cpm); wt.sd <- apply(wt.cpm, 1, sd)
k0.scGES <- (k0.cpm - wt.mean) / wt.sd
k0f4.scGES <- (k0f4.cpm - wt.mean) / wt.sd
k0.vip <- viper(k0.scGES, k0.meta)
saveRDS(k0.vip, file = 'k0-v-wt_sc-vip.rds')
k0f4.vip <- viper(k0.scGES, k0f4.meta)
saveRDS(k0f4.vip, file = 'k0f4-v-wt_sc-vip.rds')
###############

### versus WT plots
###############
vip.mat <- readRDS('k0-k0f4-v-wt_vip.rds')
## pull out mrs
k0.mrs <- sort(vip.mat[,1])
k0f4.mrs <- sort(vip.mat[,2])
## build mr set
num.mrs <- 25
k0.mrs <- c(k0.mrs[1:num.mrs], tail(k0.mrs, num.mrs))
k0f4.mrs <- c(k0f4.mrs[1:num.mrs], tail(k0f4.mrs, num.mrs))
## build row annotation
k0.unique.mrs <- setdiff(names(k0.mrs), names(k0f4.mrs))
k0f4.unique.mrs <- setdiff(names(k0f4.mrs), names(k0.mrs))
intersect.mrs <- intersect(names(k0.mrs), names(k0f4.mrs))
row.annot <- c(rep('KO', length(k0.unique.mrs)),
               rep('KOF4', length(k0f4.unique.mrs)),
               rep('Intersect', length(intersect.mrs)))
names(row.annot) <- c(k0.unique.mrs, k0f4.unique.mrs, intersect.mrs)
annot.df <- data.frame('MR' = as.factor(row.annot))
row.colors <- ClusterColors(3); names(row.colors) <- c('KO', 'KOF4', 'Intersect')
annot.colors <- list('MR' = row.colors)
## heatmap
mr.set <- c(k0.unique.mrs, k0f4.unique.mrs, intersect.mrs)
plot.mat <- vip.mat[mr.set,]
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/k0-k0f4-v-wt_mrs.jpg', height = 1000, width = 750)
pheatmap(plot.mat, main = "KO + KOF4 vs WT MRs", fontsize = 12,
         annotation_row = annot.df, annotation_colors = annot.colors,
         cluster_cols = FALSE, show_colnames = TRUE,
         cluster_rows = TRUE, show_rownames = TRUE, fontsize_row = 12,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'vip'))
dev.off()
###############

### rescue v. K0 analysis
###############
k0.seur <- readRDS('k0_seur.rds')
k0f4.seur <- readRDS('k0f4_seur.rds')
## load networks
k0.net <- readRDS('sc-nets/k0-net_pruned.rds')
k0.meta <- readRDS('sc-nets/k0_k5-net_pruned.rds')
k0f4.net <- readRDS('sc-nets/k0f4-net_pruned.rds')
k0f4.meta <- readRDS('sc-nets/k0f4_k5-net_pruned.rds')
## separate counts
k0.counts <- as.matrix(k0.seur@assays$SCT@counts)
k0f4.counts <- as.matrix(k0f4.seur@assays$SCT@counts)
## filter for shared genes
shared.genes <- intersect(rownames(k0.counts), rownames(k0f4.counts))
k0.counts <- k0.counts[shared.genes,]
k0f4.counts <- k0f4.counts[shared.genes,]
## make GES
k0f4.v.k0.ges <- DESeqGES(k0f4.counts, k0.counts)
saveRDS(k0f4.v.k0.ges, file = 'k0f4-v-k0_deseq2-ges.rds')
## run VIPER
ges.mat <- cbind(as.matrix(k0f4.v.k0.ges), as.matrix(k0f4.v.k0.ges))
colnames(ges.mat) <- c('K0F4', 'K0F4')
vip.mat <- viper(ges.mat, k0.meta)
saveRDS(vip.mat[,1], 'k0f4-v-k0_vip.rds')
###############

### versus K0 plots
###############
vip.mat <- readRDS('k0f4-v-k0_vip.rds')
## MR heatmap
num.mrs <- 50
mr.vec <- sort(vip.mat)
mr.set <- c(mr.vec[1:num.mrs], tail(mr.vec, num.mrs))
plot.mat <- as.matrix(vip.mat[names(mr.set)])
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/k0f4-v-k0_mrs.jpg', height = 750, width = 250)
pheatmap(plot.mat, main = "KOF4 vs. KO MRs", fontsize = 12,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = FALSE, show_rownames = TRUE, fontsize_row = 10,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'vip'))
dev.off()
###############

### WT MWKMeans Analysis
###############
wt.seur <- readRDS('wt_seur.rds')
wt.meta <- readRDS('sc-nets/wt-k5-net_pruned.rds')
## separate counts
wt.counts <- as.matrix(wt.seur@assays$SCT@counts)
wt.cpm <- CPMTransform(wt.counts)
## internal GES
wt.ges <- GESTransform(wt.cpm)
## viper
wt.vip <- viper(wt.ges, wt.meta)
saveRDS(wt.vip, file = 'wt-internal_sc-vip.rds')
wt.vip.dist <- as.dist(viperSimilarity(wt.vip))
saveRDS(wt.vip.dist, file = 'wt-internal_sc-vip-dist.rds')
wt.umap <- uwot::umap(t(wt.vip), metric = 'correlation')
rownames(wt.umap) <- colnames(wt.vip)
saveRDS(wt.umap, file = 'wt-internal_sc-vip-umap.rds')
## k2 cluster
wt.pam <- PamKRange(wt.vip.dist, kmin = 2, kmax = 5)
saveRDS(wt.pam, file = 'wt-internal_sc-vip-pam.rds')
## mwkmeans
k2.centers <- cbind(rowMeans(wt.vip[, which(wt.pam$clusterings$k2$clustering == 1)]),
                    rowMeans(wt.vip[, which(wt.pam$clusterings$k2$clustering == 2)]))
wt.k2.mwk <- MWKMeans(wt.vip, k2.centers)
saveRDS(wt.k2.mwk, file = 'wt-internal_sc-vip-mwk.rds')
## identify candidate MRs
phi.vec <- wt.k2.mwk[[1]][2,]
phi.prime.vec <- -4*phi.vec^2 + 4*phi.vec
phi.corr <- apply(wt.vip, 1, function(x) {cor(x, phi.vec)} )
phi.corr <- phi.corr[order(abs(phi.corr), decreasing = TRUE)]
phi.prime.corr <- apply(wt.vip, 1, function(x) {cor(x, phi.prime.vec)} )
phi.prime.corr <- phi.prime.corr[order(abs(phi.prime.corr), decreasing = TRUE)]
wt.mwk.corrs <- list('phi' = phi.corr, 'phi.prime' = phi.prime.corr)
saveRDS(wt.mwk.corrs, file = 'wt-internal_sc-vip-mwk-corrs.rds')
###############

### k0 MWKMeans Analysis
###############
k0.seur <- readRDS('k0_seur.rds')
k0.meta <- readRDS('sc-nets/k0_k5-net_pruned.rds')
## separate counts
k0.counts <- as.matrix(k0.seur@assays$SCT@counts)
k0.cpm <- CPMTransform(k0.counts)
## internal GES
k0.ges <- GESTransform(k0.cpm)
## viper
k0.vip <- viper(k0.ges, k0.meta)
saveRDS(k0.vip, file = 'k0-internal_sc-vip.rds')
k0.vip.dist <- as.dist(viperSimilarity(k0.vip))
saveRDS(k0.vip.dist, file = 'k0-internal_sc-vip-dist.rds')
k0.umap <- uwot::umap(t(k0.vip), metric = 'correlation')
rownames(k0.umap) <- colnames(k0.vip)
saveRDS(k0.umap, file = 'k0-internal_sc-vip-umap.rds')
## k2 cluster
k0.pam <- PamKRange(k0.vip.dist, kmin = 2, kmax = 5)
saveRDS(k0.pam, file = 'k0-internal_sc-vip-pam.rds')
## mwkmeans
k2.centers <- cbind(rowMeans(k0.vip[, which(k0.pam$clusterings$k2$clustering == 1)]),
                    rowMeans(k0.vip[, which(k0.pam$clusterings$k2$clustering == 2)]))
k0.k2.mwk <- MWKMeans(k0.vip, k2.centers)
saveRDS(k0.k2.mwk, file = 'k0-internal_sc-vip-mwk.rds')
## identify candidate MRs
phi.vec <- k0.k2.mwk[[1]][2,]
phi.prime.vec <- -4*phi.vec^2 + 4*phi.vec
phi.corr <- apply(k0.vip, 1, function(x) {cor(x, phi.vec)} )
phi.corr <- phi.corr[order(abs(phi.corr), decreasing = TRUE)]
phi.prime.corr <- apply(k0.vip, 1, function(x) {cor(x, phi.prime.vec)} )
phi.prime.corr <- phi.prime.corr[order(abs(phi.prime.corr), decreasing = TRUE)]
k0.mwk.corrs <- list('phi' = phi.corr, 'phi.prime' = phi.prime.corr)
saveRDS(k0.mwk.corrs, file = 'k0-internal_sc-vip-mwk-corrs.rds')
###############

### WT MWK Plots
###############
wt.seur <- readRDS('wt_seur.rds')
wt.vip <- readRDS('wt-internal_sc-vip.rds')
wt.umap <- readRDS('wt-internal_sc-vip-umap.rds')
wt.pam <- readRDS('wt-internal_sc-vip-pam.rds')
wt.k2.mwk <- readRDS('wt-internal_sc-vip-mwk.rds')
wt.mwk.corrs <- readRDS('wt-internal_sc-vip-mwk-corrs.rds')
marker.list <- c('Plaur', 'Angpt2', 'Lcp2', 'Cxcr4', 'Apln', 'Kcne3', 'Dll4', 'Esm1',
                 'Tie2', 'Tek', 'Jag1', 'Apj')
## re-generate GES
wt.ges <- GESTransform(CPMTransform(as.matrix(wt.seur@assays$SCT@counts)))
## pull out cluster / phi / phi' vectors
clust.vec <- wt.pam$clusterings$k2$clustering
phi.vec <- wt.k2.mwk[[1]][2,]
phi.prime.vec <- -4*phi.vec^2 + 4*phi.vec
cell.order <- names(sort(phi.vec))
## make annotation data frame
annot.df <- data.frame('Cluster' = clust.vec[cell.order],
                       'phi' = phi.vec[cell.order],
                       'phi.prime' = phi.prime.vec[cell.order])
clust.colors <- ClusterColors(2); names(clust.colors) <- unique(clust.vec)
annot.color <- list('Cluster' = clust.colors,
                    'phi' = colorRampPalette(c('grey', 'forestgreen'))(100),
                    'phi.prime' = colorRampPalette(c('grey', 'purple'))(100))
## plot marker activity
plot.mat <- wt.vip[intersect(marker.list, rownames(wt.vip)), cell.order]
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/wt-internal_mwk-marker-vip.jpg', height = 750, width = 750)
pheatmap(plot.mat, main = "MWKMeans: Marker Activity", fontsize = 12,
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = FALSE, show_rownames = TRUE, fontsize_row = 12,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'vip'))
dev.off()
## plot marker expression (GES)
plot.mat <- wt.ges[intersect(marker.list, rownames(wt.ges)), cell.order]
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/wt-internal_mwk-marker-exp.jpg', height = 750, width = 750)
pheatmap(plot.mat, main = "MWKMeans: Marker Expression", fontsize = 12, scale = 'column',
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = FALSE, show_rownames = TRUE, fontsize_row = 12,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'gexp'))
dev.off()
## plot candidate MRs (phi)
num.mrs <- 50
mr.set <- names(wt.mwk.corrs$phi[1:num.mrs])
plot.mat <- wt.vip[intersect(mr.set, rownames(wt.vip)), cell.order]
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/wt-internal_mwk-phi-mrs.jpg', height = 1000, width = 1000)
pheatmap(plot.mat, main = "MWKMeans: Phi MRs", fontsize = 12,
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = TRUE, show_rownames = TRUE, fontsize_row = 10,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'vip'))
dev.off()
## plot candidate MRs (phi-prime)
num.mrs <- 50
mr.set <- names(wt.mwk.corrs$phi.prime[1:num.mrs])
plot.mat <- wt.vip[intersect(mr.set, rownames(wt.vip)), cell.order]
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/wt-internal_mwk-phi-prime-mrs.jpg', height = 1000, width = 1000)
pheatmap(plot.mat, main = "MWKMeans: Phi Prime MRs", fontsize = 12,
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = TRUE, show_rownames = TRUE, fontsize_row = 10,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'vip'))
dev.off()
## phi UMAP
plot.dat <- data.frame('UMAP1' = wt.umap[,1], 'UMAP2' = wt.umap[,2],
                       'phi' = phi.vec)
jpeg('plots/wt-internal_mwk-vip-umap.jpg', height = 750, width = 1250)
ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = phi)) +
  ggtitle('WT: Phi') + 
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  theme(text = element_text(size = 20))
dev.off()
## marker umap
marker.set <- intersect(marker.list, rownames(wt.vip))
plot.dat <- data.frame('UMAP1' = wt.umap[,1], 'UMAP2' = wt.umap[,2])
plot.dat <- cbind(plot.dat, t(wt.vip[marker.set,]))
plot.list <- list()
for (m in marker.set) {
  plot.list[[m]] <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes_string(color = m)) +
    ggtitle(m) + theme(text = element_text(size = 12)) +
    scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0, high = 'red')
}
marker.plot <- ggarrange(plotlist = plot.list, nrow = 2, ncol = length(marker.set) %/% 2 + 1)
jpeg('plots/wt-internal_mwk-marker-vip-umap.jpg', height = 1000, 
     width = 500 * (length(marker.set) %/% 2 + 1))
print(marker.plot)
dev.off()
###############

### k0 MWK Plots
###############
k0.seur <- readRDS('k0_seur.rds')
k0.vip <- readRDS('k0-internal_sc-vip.rds')
k0.umap <- readRDS('k0-internal_sc-vip-umap.rds')
k0.pam <- readRDS('k0-internal_sc-vip-pam.rds')
k0.k2.mwk <- readRDS('k0-internal_sc-vip-mwk.rds')
k0.mwk.corrs <- readRDS('k0-internal_sc-vip-mwk-corrs.rds')
marker.list <- c('Plaur', 'Angpt2', 'Lcp2', 'Cxcr4', 'Apln', 'Kcne3', 'Dll4', 'Esm1',
                 'Tie2', 'Tek', 'Jag1', 'Apj')
## re-generate GES
k0.ges <- GESTransform(CPMTransform(as.matrix(k0.seur@assays$SCT@counts)))
## pull out cluster / phi / phi' vectors
clust.vec <- k0.pam$clusterings$k2$clustering
phi.vec <- k0.k2.mwk[[1]][2,]
phi.prime.vec <- -4*phi.vec^2 + 4*phi.vec
cell.order <- names(sort(phi.vec))
## make annotation data frame
annot.df <- data.frame('Cluster' = clust.vec[cell.order],
                       'phi' = phi.vec[cell.order],
                       'phi.prime' = phi.prime.vec[cell.order])
clust.colors <- ClusterColors(2); names(clust.colors) <- unique(clust.vec)
annot.color <- list('Cluster' = clust.colors,
                    'phi' = colorRampPalette(c('grey', 'forestgreen'))(100),
                    'phi.prime' = colorRampPalette(c('grey', 'purple'))(100))
## plot marker activity
plot.mat <- k0.vip[intersect(marker.list, rownames(k0.vip)), cell.order]
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/k0-internal_mwk-marker-vip.jpg', height = 750, width = 750)
pheatmap(plot.mat, main = "MWKMeans: Marker Activity", fontsize = 12,
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = FALSE, show_rownames = TRUE, fontsize_row = 12,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'vip'))
dev.off()
## plot marker expression (GES)
plot.mat <- k0.ges[intersect(marker.list, rownames(k0.ges)), cell.order]
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/k0-internal_mwk-marker-exp.jpg', height = 750, width = 750)
pheatmap(plot.mat, main = "MWKMeans: Marker Expression", fontsize = 12, scale = 'column',
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = FALSE, show_rownames = TRUE, fontsize_row = 12,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'gexp'))
dev.off()
## plot candidate MRs (phi)
num.mrs <- 50
mr.set <- names(k0.mwk.corrs$phi[1:num.mrs])
plot.mat <- k0.vip[intersect(mr.set, rownames(k0.vip)), cell.order]
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/k0-internal_mwk-phi-mrs.jpg', height = 1000, width = 1000)
pheatmap(plot.mat, main = "MWKMeans: Phi MRs", fontsize = 12,
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = TRUE, show_rownames = TRUE, fontsize_row = 10,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'vip'))
dev.off()
## plot candidate MRs (phi-prime)
num.mrs <- 50
mr.set <- names(k0.mwk.corrs$phi.prime[1:num.mrs])
plot.mat <- k0.vip[intersect(mr.set, rownames(k0.vip)), cell.order]
mat.breaks <- QuantileBreaks(plot.mat, 100)
jpeg('plots/k0-internal_mwk-phi-prime-mrs.jpg', height = 1000, width = 1000)
pheatmap(plot.mat, main = "MWKMeans: Phi Prime MRs", fontsize = 12,
         annotation_col = annot.df, annotation_colors = annot.color,
         cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = TRUE, show_rownames = TRUE, fontsize_row = 10,
         breaks = mat.breaks, color = ColorLevels(length(mat.breaks) - 1, 'vip'))
dev.off()
## phi UMAP
plot.dat <- data.frame('UMAP1' = k0.umap[,1], 'UMAP2' = k0.umap[,2],
                       'phi' = phi.vec)
jpeg('plots/k0-internal_mwk-vip-umap.jpg', height = 750, width = 1250)
ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = phi)) +
  ggtitle('KO: Phi') + 
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  theme(text = element_text(size = 20))
dev.off()
## marker umap
marker.set <- intersect(marker.list, rownames(k0.vip))
plot.dat <- data.frame('UMAP1' = k0.umap[,1], 'UMAP2' = k0.umap[,2])
plot.dat <- cbind(plot.dat, t(k0.vip[marker.set,]))
plot.list <- list()
for (m in marker.set) {
  plot.list[[m]] <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes_string(color = m)) +
    ggtitle(m) + theme(text = element_text(size = 12)) +
    scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0, high = 'red')
}
marker.plot <- ggarrange(plotlist = plot.list, nrow = 2, ncol = length(marker.set) %/% 2 + 1)
jpeg('plots/k0-internal_mwk-marker-vip-umap.jpg', height = 1000, 
     width = 500 * (length(marker.set) %/% 2 + 1))
print(marker.plot)
dev.off()
###############
