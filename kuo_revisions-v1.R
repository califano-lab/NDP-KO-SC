setwd('C://Users/Lukas Vlahos/OneDrive/ac_lab/kuo-colab/')
library(Seurat)
library(PISCES)
library(NaRnEA)
library(reshape2)
library(ggplot2)
library(mclust)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
phi_adjust <- function(x) {
  return(x**2 / (x**2 + (1 - x)**2))
}

### marker sets
###############
zarkada.markers <- list('Venous' = c('Hmgn2', 'Ptgis'),
                        'Proliferative' = c('Mki67', 'Birc5'),
                        'Capillary' = c('Mfsd2a', 'Slc22a8'),
                        'Arterial' = c('Unc5b', 'Bmx'),
                        'Tip' = c('Kcne3', 'Angpt2'),
                        'S-Tip' = c('Angpt2', 'Plvap'),
                        'D-Tip' = c('Cldn5', 'Mfsd2a'))

kuo.markers <- list('tip' = c('Plaur', 'Angpt2', 'Lcp2', 'Cxcr4', 'Apln', 
                              'Kcne3', 'Dll4', 'Esm1', 'Vegfr2', 'Vegfr3', 
                              'Kdr', 'Flk1', 'Flt4'),
                    'stalk' = c('Tie2', 'Tek', 'Jag1', 'Apj', 'Vegfr1', 'Flt1'))
## zarkada gsea
zarkada.gsea.df <- read.csv(file = 'revisions/zarkada_gsea-sets.csv')
zarkada.regulons <- apply(zarkada.gsea.df, 2, function(x) {
  gene.set <- x[which(x != '')]
  gene.set <- unique(gene.set)
  am.vec <- rep(1, length(gene.set))
  names(am.vec) <- gene.set
  aw.vec <- rep(1, length(gene.set))
  names(aw.vec) <- gene.set
  return(list('am' = am.vec, 'aw' = aw.vec))
})
names(zarkada.regulons) <- c('BRB', 'WNT', 'WNT.Targets', 'Glycolysis', 'TCA.Cycle',
                             'TGFb.Targets', 'ECM', 'ECM.Regulators', 'D.Tip', 'S.Tip')
use.regs <- c('BRB', 'WNT', 'WNT.Targets', 'Glycolysis', 'TCA.Cycle',
              'TGFb.Targets', 'D.Tip', 'S.Tip')
zarkada.regulons <- zarkada.regulons[use.regs]
## create marker dfs
zarkada.marker.df <- data.frame('Gene' = unlist(zarkada.markers),
                                'Group' = rep(names(zarkada.markers), sapply(zarkada.markers, length)))
kuo.marker.df <- data.frame('Gene' = unlist(kuo.markers),
                            'Group' = rep(names(kuo.markers), sapply(kuo.markers, length)))
## phi adjust function
phi_adjust <- function(x) {
  return(x**2 / (x**2 + (1 - x)**2))
}
## set annotation colors
cell.group.cols <- group_colors(3); names(cell.group.cols) <- c('Stalk', 'Intermediate', 'Tip')
phi.col.func <-  colorRamp2(c(0, 0.5, 1), c('white', 'grey', 'black'))
zarkada.group.col <- group_colors(length(unique(zarkada.marker.df$Group))); names(zarkada.group.col) <- unique(zarkada.marker.df$Group)
kuo.group.col <- group_colors(length(unique(kuo.marker.df$Group))); names(kuo.group.col) <- unique(kuo.marker.df$Group)
###############

### zarkada differentail expression process
###############
p6.dif.exp <- read.csv('revisions/zarkada_p6-genes.csv', header = TRUE)
colnames(p6.dif.exp) <- c('Tip', 'Arterial', 'Capillary', 'Proliferative', 'Venous')
p10.dif.exp <- read.csv('revisions/zarkada_p10-genes.csv', header = TRUE)
colnames(p10.dif.exp) <- c('Tip', 'Arterial', 'Capillary', 'Proliferative', 'Venous')
## get overlap for each group
dif.gene.sets <- sapply(colnames(p6.dif.exp), function(x) {
  dif.genes <- intersect(p6.dif.exp[, x], p10.dif.exp[, x])
  return(dif.genes)
})
dif.gene.df <- Reduce(rbind, lapply(names(dif.gene.sets), function(x) {
  return(data.frame('Gene' = dif.gene.sets[[x]][1:5], 'Group' = x))
}))
## get s-tip and d-tip markers
ds.dif.exp <- read.csv('revisions/zarkada_d-s-tip-genes.csv', header = TRUE)
colnames(ds.dif.exp) <- c('D-Tip', 'S-Tip')
ds.dif.gene.df <- Reduce(rbind, lapply(colnames(ds.dif.exp), function(x) {
  return(data.frame('Gene' = ds.dif.exp[[x]][1:5], 'Group' = x))
}))
###############

###########################################################################
# Functions
###########################################################################

#' Zarkada GSEA function
zarkada_gsea <- function(counts.mat, clust.vec, plot.title, x.title) {
  
  ## create signature and run narnea
  sig.vecs <- list()
  narnea.vecs <- list()
  for (cn in sort(unique(clust.vec))) {
    print(cn)
    # calculate signature
    test.cells <- names(clust.vec)[which(clust.vec == cn)]
    ref.cells <- names(clust.vec)[which(clust.vec != cn)]
    cn.ges <- apply(counts.mat, 1, function(x) {
      w.test <- wilcox.test(x[test.cells], x[ref.cells], alternative = "two.sided")
      rbs.cor <- 2 * w.test$statistic / (length(test.cells) * length(ref.cells)) - 1
      return(qnorm(1 - w.test$p.val) * sign(rbs.cor))
    })
    cn.ges <- cn.ges[!is.na(cn.ges)]
    sig.vecs[[as.character(cn)]] <- cn.ges
    # run narnea for each regulon
    cn.narnea <- sapply(zarkada.regulons, function(reg.obj) {
      # filter regulon
      use.genes <- intersect(names(cn.ges), names(reg.obj$aw))
      if (length(use.genes) < 25) {
        return(list('nes' = NA, 'pes' = NA))
      }
      aw.vec <- reg.obj$aw[use.genes]
      am.vec <- reg.obj$am[use.genes]
      # run narnea
      x.narnea <- NaRnEA(cn.ges, association.weight = aw.vec, association.mode = am.vec,
                         minimum.size = 20, ledge = FALSE)
    })
    narnea.vecs[[as.character(cn)]] <- cn.narnea
  }
  
  ## generate dot plot
  pes.mat <- sapply(narnea.vecs, function(x) {
    sapply(x, function(y) {y$pes})
  })
  row.has.na <- apply(pes.mat, 1, function(x) {length(which(is.na(x))) > 0})
  pes.mat <- pes.mat[!row.has.na,]
  nes.mat <- sapply(narnea.vecs, function(x) {
    sapply(x, function(y) {y$nes})
  })
  row.has.na <- apply(nes.mat, 1, function(x) {length(which(is.na(x))) > 0})
  nes.mat <- nes.mat[!row.has.na,]
  # plot
  plot.df <- melt(pes.mat)
  colnames(plot.df) <- c('GS', 'Cluster', 'PES')
  plot.df$NES <- melt(nes.mat)[,3]
  plot.df$Cluster <- as.factor(plot.df$Cluster)
  plot.obj <- ggplot(plot.df, aes(x = Cluster, y = GS)) + 
    geom_point(aes(color = PES, size = abs(NES) + 1)) +
    scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) + 
    labs(x = x.title, y = "Gene Set", title = plot.title) +
    guides(size = guide_legend(title = "|NES|"))
  
  ## return
  return(list('narnea' = narnea.vecs, 'plot' = plot.obj))
}

###########################################################################
###########################################################################



###########################################################################
# Gene Expression Analysis - Clustering, Markers, GSEA
###########################################################################

### WT - initial analysis + network prep
###############
wt.raw <- Read10X('sc/WT/filtered_feature_bc_matrix/')
wt.raw <- as.matrix(wt.raw)
## qc
jpeg('revisions/WT/wt_qc-plot.jpg', height = 6, width = 8, units = 'in', res = 200)
qc_plot(wt.raw, species = 'mur', genes = 'symb')
dev.off()
filt.counts <- qc_filt(wt.raw, min.depth = 2000, max.depth = 50000,
                       min.genes = 1000,
                       max.mt = 0.1)
saveRDS(filt.counts, file = 'revisions/WT/wt_filt-counts.rds')
## get cpm counts
cpm.counts <- cpm_norm(filt.counts, l2 = TRUE)
saveRDS(cpm.counts, file = 'revisions/WT/wt_cpm-counts.rds')
## create seurat object w/ discrete clustering
seurat.obj <- CreateSeuratObject(counts = filt.counts, assay = 'RNA', project = 'wt.seur')
seurat.obj <- PercentageFeatureSet(seurat.obj, pattern = "^mt-", col.name = "percent.mt")
seurat.obj <- SCTransform(seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindClusters(seurat.obj, verbose = FALSE)
saveRDS(seurat.obj, file = 'revisions/WT/wt_seur.rds')
## marker heatmap w/ zarkada markers
zarkada.marker.df <- data.frame('Gene' = unlist(zarkada.markers),
                                'Group' = rep(names(zarkada.markers), sapply(zarkada.markers, length)))
jpeg('revisions/WT/wt_seurat-clust_zarkada-markers_heatmap.jpg',
     height = 10, width = 10, res = 300, units = 'in')
cluster_mr_heatmap(seurat.obj@assays$SCT@scale.data, dat.type = 'gexp',
                   clust.vec = seurat.obj$seurat_clusters,
                   marker.set = zarkada.marker.df, 
                   plot.title = 'Cluster Marker Expression - Zarkada Markers')
dev.off()
## marker heatmap w/ Kuo markers
kuo.marker.df <- data.frame('Gene' = unlist(kuo.markers),
                            'Group' = rep(names(kuo.markers), sapply(kuo.markers, length)))
jpeg('revisions/WT/wt_seurat-clust_kuo-markers_heatmap.jpg',
     height = 10, width = 10, res = 300, units = 'in')
cluster_mr_heatmap(seurat.obj@assays$SCT@scale.data, dat.type = 'gexp',
                   clust.vec = seurat.obj$seurat_clusters,
                   marker.set = kuo.marker.df, 
                   plot.title = 'Cluster Marker Expression - Kuo Markers')
dev.off()
## gsea analysis for each cluster versus the rest
clust.vec <- seurat.obj$seurat_clusters
signature.vecs <- list()
narnea.vecs <- list()
for (cn in sort(unique(clust.vec))) {
  print(cn)
  # calculate signature
  test.cells <- names(clust.vec)[which(clust.vec == cn)]
  ref.cells <- names(clust.vec)[which(clust.vec != cn)]
  cn.ges <- apply(cpm.counts, 1, function(x) {
    w.test <- wilcox.test(x[test.cells], x[ref.cells], alternative = "two.sided")
    rbs.cor <- 2 * w.test$statistic / (length(test.cells) * length(ref.cells)) - 1
    return(qnorm(1 - w.test$p.val) * sign(rbs.cor))
  })
  cn.ges <- cn.ges[!is.na(cn.ges)]
  signature.vecs[[as.character(cn)]] <- cn.ges
  # run narnea for each regulon
  cn.narnea <- sapply(zarkada.regulons, function(reg.obj) {
    # filter regulon
    use.genes <- intersect(names(cn.ges), names(reg.obj$aw))
    if (length(use.genes) < 25) {
      return(list('nes' = NA, 'pes' = NA))
    }
    aw.vec <- reg.obj$aw[use.genes]
    am.vec <- reg.obj$am[use.genes]
    # run narnea
    x.narnea <- NaRnEA(cn.ges, association.weight = aw.vec, association.mode = am.vec,
                       minimum.size = 20, ledge = FALSE)
  })
  narnea.vecs[[as.character(cn)]] <- cn.narnea
}
saveRDS(signature.vecs, file = 'revisions/WT/wt_seurat-clust_signatures.rds')
saveRDS(narnea.vecs, file = 'revisions/WT/wt_seurat-clust_zarkada-gs-narnea.rds')
## make gene set plot
pes.mat <- sapply(narnea.vecs, function(x) {
  sapply(x, function(y) {y$pes})
})
row.has.na <- apply(pes.mat, 1, function(x) {length(which(is.na(x))) > 0})
pes.mat <- pes.mat[!row.has.na,]
nes.mat <- sapply(narnea.vecs, function(x) {
  sapply(x, function(y) {y$nes})
})
row.has.na <- apply(nes.mat, 1, function(x) {length(which(is.na(x))) > 0})
nes.mat <- nes.mat[!row.has.na,]
# make dot plot
plot.df <- melt(pes.mat)
colnames(plot.df) <- c('GS', 'Cluster', 'PES')
plot.df$NES <- melt(nes.mat)[,3]
plot.df$Cluster <- as.factor(plot.df$Cluster)
ggplot(plot.df, aes(x = Cluster, y = GS)) + 
  geom_point(aes(color = PES, size = abs(NES))) +
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) + 
  labs(x = "Seurat Cluster", y = "Gene Set", title = "Zarkada Gene Set Enrichment") +
  guides(size = guide_legend(title = "|NES|"))
ggsave('revisions/WT/wt_seurat-clust_zarkada-gs-narnea_dot-plot.jpg',
       width = 10, height = 6, units = 'in', dpi = 250)
## network data prep
wt.dist <- dist(seurat.obj@reductions$pca@cell.embeddings)
wt.metacell <- make_metacells(filt.counts, wt.dist, subset = 500)
saveRDS(wt.metacell[[1]], file = 'revisions/WT/a3_network/wt_k5-meta.rds')
###############

### KO - initial analysis
###############
ko.raw <- Read10X('sc/KO/filtered_feature_bc_matrix/')
ko.raw <- as.matrix(ko.raw)
## qc
jpeg('revisions/KO/ko_qc-plot.jpg', height = 6, width = 8, units = 'in', res = 200)
qc_plot(ko.raw, species = 'mur', genes = 'symb')
dev.off()
filt.counts <- qc_filt(ko.raw, min.depth = 2000, max.depth = 50000,
                       min.genes = 1000,
                       max.mt = 0.1)
saveRDS(filt.counts, file = 'revisions/KO/ko_filt-counts.rds')
## get cpm counts
cpm.counts <- cpm_norm(filt.counts, l2 = TRUE)
saveRDS(cpm.counts, file = 'revisions/KO/ko_cpm-counts.rds')
## create seurat object w/ discrete clustering
seurat.obj <- CreateSeuratObject(counts = filt.counts, assay = 'RNA', project = 'ko.seur')
seurat.obj <- PercentageFeatureSet(seurat.obj, pattern = "^mt-", col.name = "percent.mt")
seurat.obj <- SCTransform(seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindClusters(seurat.obj, verbose = FALSE)
saveRDS(seurat.obj, file = 'revisions/KO/ko_seur.rds')
## marker heatmap w/ zarkada markers
zarkada.marker.df <- data.frame('Gene' = unlist(zarkada.markers),
                                'Group' = rep(names(zarkada.markers), sapply(zarkada.markers, length)))
jpeg('revisions/KO/ko_seurat-clust_zarkada-markers_heatmap.jpg',
     height = 10, width = 10, res = 300, units = 'in')
cluster_mr_heatmap(seurat.obj@assays$SCT@scale.data, dat.type = 'gexp',
                   clust.vec = seurat.obj$seurat_clusters,
                   marker.set = zarkada.marker.df, 
                   plot.title = 'Cluster Marker Expression - Zarkada Markers')
dev.off()
## marker heatmap w/ Kuo markers
kuo.marker.df <- data.frame('Gene' = unlist(kuo.markers),
                            'Group' = rep(names(kuo.markers), sapply(kuo.markers, length)))
jpeg('revisions/KO/ko_seurat-clust_kuo-markers_heatmap.jpg',
     height = 10, width = 10, res = 300, units = 'in')
cluster_mr_heatmap(seurat.obj@assays$SCT@scale.data, dat.type = 'gexp',
                   clust.vec = seurat.obj$seurat_clusters,
                   marker.set = kuo.marker.df, 
                   plot.title = 'Cluster Marker Expression - Kuo Markers')
dev.off()
## gsea analysis for each cluster versus the rest
clust.vec <- seurat.obj$seurat_clusters
signature.vecs <- list()
narnea.vecs <- list()
for (cn in sort(unique(clust.vec))) {
  print(cn)
  # calculate signature
  test.cells <- names(clust.vec)[which(clust.vec == cn)]
  ref.cells <- names(clust.vec)[which(clust.vec != cn)]
  cn.ges <- apply(cpm.counts, 1, function(x) {
    w.test <- wilcox.test(x[test.cells], x[ref.cells], alternative = "two.sided")
    rbs.cor <- 2 * w.test$statistic / (length(test.cells) * length(ref.cells)) - 1
    return(qnorm(1 - w.test$p.val) * sign(rbs.cor))
  })
  cn.ges <- cn.ges[!is.na(cn.ges)]
  signature.vecs[[as.character(cn)]] <- cn.ges
  # run narnea for each regulon
  cn.narnea <- sapply(zarkada.regulons, function(reg.obj) {
    # filter regulon
    use.genes <- intersect(names(cn.ges), names(reg.obj$aw))
    if (length(use.genes) < 25) {
      return(list('nes' = NA, 'pes' = NA))
    }
    aw.vec <- reg.obj$aw[use.genes]
    am.vec <- reg.obj$am[use.genes]
    # run narnea
    x.narnea <- NaRnEA(cn.ges, association.weight = aw.vec, association.mode = am.vec,
                       minimum.size = 20, ledge = FALSE)
  })
  narnea.vecs[[as.character(cn)]] <- cn.narnea
}
saveRDS(signature.vecs, file = 'revisions/KO/ko_seurat-clust_signatures.rds')
saveRDS(narnea.vecs, file = 'revisions/KO/ko_seurat-clust_zarkada-gs-narnea.rds')
## make gene set plot
pes.mat <- sapply(narnea.vecs, function(x) {
  sapply(x, function(y) {y$pes})
})
row.has.na <- apply(pes.mat, 1, function(x) {length(which(is.na(x))) > 0})
pes.mat <- pes.mat[!row.has.na,]
nes.mat <- sapply(narnea.vecs, function(x) {
  sapply(x, function(y) {y$nes})
})
row.has.na <- apply(nes.mat, 1, function(x) {length(which(is.na(x))) > 0})
nes.mat <- nes.mat[!row.has.na,]
# make dot plot
plot.df <- melt(pes.mat)
colnames(plot.df) <- c('GS', 'Cluster', 'PES')
plot.df$NES <- melt(nes.mat)[,3]
plot.df$Cluster <- as.factor(plot.df$Cluster)
ggplot(plot.df, aes(x = Cluster, y = GS)) + 
  geom_point(aes(color = PES, size = abs(NES))) +
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) + 
  labs(x = "Seurat Cluster", y = "Gene Set", title = "KO - Zarkada Gene Set Enrichment") +
  guides(size = guide_legend(title = "|NES|"))
ggsave('revisions/KO/ko_seurat-clust_zarkada-gs-narnea_dot-plot.jpg',
       width = 10, height = 6, units = 'in', dpi = 250)
###############

### KOF4 - initial analysis
###############
kof4.raw <- Read10X('sc/KOF4/filtered_feature_bc_matrix/')
kof4.raw <- as.matrix(kof4.raw)
## qc
jpeg('revisions/KOF4/kof4_qc-plot.jpg', height = 6, width = 8, units = 'in', res = 200)
qc_plot(kof4.raw, species = 'mur', genes = 'symb')
dev.off()
filt.counts <- qc_filt(kof4.raw, min.depth = 2000, max.depth = 40000,
                       min.genes = 1000,
                       max.mt = 0.1)
saveRDS(filt.counts, file = 'revisions/KOF4/kof4_filt-counts.rds')
## get cpm counts
cpm.counts <- cpm_norm(filt.counts, l2 = TRUE)
saveRDS(cpm.counts, file = 'revisions/KOF4/kof4_cpm-counts.rds')
## create seurat object w/ discrete clustering
seurat.obj <- CreateSeuratObject(counts = filt.counts, assay = 'RNA', project = 'ko.seur')
seurat.obj <- PercentageFeatureSet(seurat.obj, pattern = "^mt-", col.name = "percent.mt")
seurat.obj <- SCTransform(seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:30, verbose = FALSE)
seurat.obj <- FindClusters(seurat.obj, verbose = FALSE)
saveRDS(seurat.obj, file = 'revisions/KOF4/kof4_seur.rds')
## marker heatmap w/ zarkada markers
zarkada.marker.df <- data.frame('Gene' = unlist(zarkada.markers),
                                'Group' = rep(names(zarkada.markers), sapply(zarkada.markers, length)))
jpeg('revisions/KOF4/kof4_seurat-clust_zarkada-markers_heatmap.jpg',
     height = 10, width = 10, res = 300, units = 'in')
cluster_mr_heatmap(seurat.obj@assays$SCT@scale.data, dat.type = 'gexp',
                   clust.vec = seurat.obj$seurat_clusters,
                   marker.set = zarkada.marker.df, 
                   plot.title = 'Cluster Marker Expression - Zarkada Markers')
dev.off()
## marker heatmap w/ Kuo markers
kuo.marker.df <- data.frame('Gene' = unlist(kuo.markers),
                            'Group' = rep(names(kuo.markers), sapply(kuo.markers, length)))
jpeg('revisions/KOF4/kof4_seurat-clust_kuo-markers_heatmap.jpg',
     height = 10, width = 10, res = 300, units = 'in')
cluster_mr_heatmap(seurat.obj@assays$SCT@scale.data, dat.type = 'gexp',
                   clust.vec = seurat.obj$seurat_clusters,
                   marker.set = kuo.marker.df, 
                   plot.title = 'Cluster Marker Expression - Kuo Markers')
dev.off()
## gsea analysis for each cluster versus the rest
clust.vec <- seurat.obj$seurat_clusters
signature.vecs <- list()
narnea.vecs <- list()
for (cn in sort(unique(clust.vec))) {
  print(cn)
  # calculate signature
  test.cells <- names(clust.vec)[which(clust.vec == cn)]
  ref.cells <- names(clust.vec)[which(clust.vec != cn)]
  cn.ges <- apply(cpm.counts, 1, function(x) {
    w.test <- wilcox.test(x[test.cells], x[ref.cells], alternative = "two.sided")
    rbs.cor <- 2 * w.test$statistic / (length(test.cells) * length(ref.cells)) - 1
    return(qnorm(1 - w.test$p.val) * sign(rbs.cor))
  })
  cn.ges <- cn.ges[!is.na(cn.ges)]
  signature.vecs[[as.character(cn)]] <- cn.ges
  # run narnea for each regulon
  cn.narnea <- sapply(zarkada.regulons, function(reg.obj) {
    # filter regulon
    use.genes <- intersect(names(cn.ges), names(reg.obj$aw))
    if (length(use.genes) < 25) {
      return(list('nes' = NA, 'pes' = NA))
    }
    aw.vec <- reg.obj$aw[use.genes]
    am.vec <- reg.obj$am[use.genes]
    # run narnea
    x.narnea <- NaRnEA(cn.ges, association.weight = aw.vec, association.mode = am.vec,
                       minimum.size = 20, ledge = FALSE)
  })
  narnea.vecs[[as.character(cn)]] <- cn.narnea
}
saveRDS(signature.vecs, file = 'revisions/KOF4/kof4_seurat-clust_signatures.rds')
saveRDS(narnea.vecs, file = 'revisions/KOF4/kof4_seurat-clust_zarkada-gs-narnea.rds')
## make gene set plot
pes.mat <- sapply(narnea.vecs, function(x) {
  sapply(x, function(y) {y$pes})
})
row.has.na <- apply(pes.mat, 1, function(x) {length(which(is.na(x))) > 0})
pes.mat <- pes.mat[!row.has.na,]
nes.mat <- sapply(narnea.vecs, function(x) {
  sapply(x, function(y) {y$nes})
})
row.has.na <- apply(nes.mat, 1, function(x) {length(which(is.na(x))) > 0})
nes.mat <- nes.mat[!row.has.na,]
# make dot plot
plot.df <- melt(pes.mat)
colnames(plot.df) <- c('GS', 'Cluster', 'PES')
plot.df$NES <- melt(nes.mat)[,3]
plot.df$Cluster <- as.factor(plot.df$Cluster)
ggplot(plot.df, aes(x = Cluster, y = GS)) + 
  geom_point(aes(color = PES, size = abs(NES))) +
  scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0) + 
  labs(x = "Seurat Cluster", y = "Gene Set", title = "KOF4 - Zarkada Gene Set Enrichment") +
  guides(size = guide_legend(title = "|NES|"))
ggsave('revisions/KOF4/kof4_seurat-clust_zarkada-gs-narnea_dot-plot.jpg',
       width = 10, height = 6, units = 'in', dpi = 250)
###############

###########################################################################
###########################################################################



###########################################################################
# Anchored Analysis
###########################################################################

### normalization, clustering, gsea
###############
wt.counts <- readRDS('revisions/WT/wt_filt-counts.rds')
ko.counts <- readRDS('revisions/KO/ko_filt-counts.rds')
kof4.counts <- readRDS('revisions/KOF4/kof4_filt-counts.rds')
## get shared genes
shared.genes <- Reduce(intersect, list('wt' = rownames(wt.counts), 'ko' = rownames(ko.counts), 'kof4' = rownames(kof4.counts)))
## modify cell names
colnames(wt.counts) <- paste('WT', colnames(wt.counts), sep = '_')
colnames(ko.counts) <- paste('KO', colnames(ko.counts), sep = '_')
colnames(kof4.counts) <- paste('KOF4', colnames(kof4.counts), sep = '_')
## create anchored object
wt.obj <- CreateSeuratObject(counts = wt.counts, project = 'WT')
wt.obj <- PercentageFeatureSet(wt.obj, pattern = "^mt-", col.name = "percent.mt")
wt.obj <- SCTransform(wt.obj, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = nrow(wt.counts))
ko.obj <- CreateSeuratObject(counts = ko.counts, project = 'KO')
ko.obj <- PercentageFeatureSet(ko.obj, pattern = "^mt-", col.name = "percent.mt")
ko.obj <- SCTransform(ko.obj, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = nrow(ko.counts))
kof4.obj <- CreateSeuratObject(counts = kof4.counts, project = 'KOF4')
kof4.obj <- PercentageFeatureSet(kof4.obj, pattern = "^mt-", col.name = "percent.mt")
kof4.obj <- SCTransform(kof4.obj, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = nrow(kof4.counts))
obj.list <- list('wt' = wt.obj, 'ko' = ko.obj, 'kof4' = kof4.obj)
int.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = length(shared.genes))
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = int.features)
int.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", 
                                      anchor.features = int.features, verbose = T, reference = 1)
anchored.obj  <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = T)
## clustering analysis
anchored.obj <- ScaleData(anchored.obj, verbose = FALSE)
anchored.obj <- RunPCA(anchored.obj, npcs = 30, verbose = FALSE)
anchored.obj <- RunUMAP(anchored.obj, reduction = "pca", dims = 1:30)
anchored.obj <- FindNeighbors(anchored.obj, reduction = "pca", dims = 1:30)
anchored.obj <- FindClusters(anchored.obj, resolution = 0.5)
saveRDS(anchored.obj, file = 'revisions/anchored/anchored_seur.rds')
## marker analysis
cluster.markers <- list()
for (clust.name in levels(anchored.obj$seurat_clusters)) {
  print(clust.name)
  cn.markers <- FindMarkers(anchored.obj, ident.1 = clust.name)
  cn.markers <- cn.markers[order(cn.markers$p_val_adj, -abs(cn.markers$avg_log2FC)),]
  cluster.markers[[as.character(clust.name)]] <- cn.markers
}
saveRDS(cluster.markers, file = 'revisions/anchored/anchored_seurat-clust_markers.rds')
## gsea
z.gsea <- zarkada_gsea(as.matrix(anchored.obj@assays$SCT@scale.data), clust.vec = anchored.obj$seurat_clusters,
                       plot.title = 'Anchored Analysis - Zarkada GSEA')
saveRDS(z.gsea, file = 'revisions/anchored/anchored_zarkada-gsea.rds')
jpeg('revisions/anchored/anchored_zarkada-gsea.jpg', height = 8, width = 10, units = 'in', res = 250)
print(z.gsea$plot)
dev.off()
## cluster change statistical significance
clust.exp.table <- table(ct.exp, anchored.obj$seurat_clusters)
wt.v.ko <- chisq.test(x = t(clust.exp.table[c('WT', 'KO'),]))
ko.v.kof4 <- chisq.test(x = t(clust.exp.table[c('KO', 'KOF4'),]))
chisq.tests <- list('wt.v.ko' = wt.v.ko, 'ko.v.kof4' = ko.v.kof4)
saveRDS(chisq.tests, file = 'revisions/anchored/anchored_chisq-tests.rds')
# make dot plot for cluster differences
wt.v.ko.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(wt.v.ko$stdres), lower.tail = FALSE)), ncol = 2)
ko.v.kof4.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(ko.v.kof4$stdres), lower.tail = FALSE)), ncol = 2)
pval.vec <- c(wt.v.ko.adjp[,2], ko.v.kof4.adjp[,2])
sign.vec <- c(sign(wt.v.ko$stdres[,2]), sign(ko.v.kof4$stdres[,2]))
dif.vec <- rep('NS', 2*ncol(clust.exp.table))
dif.vec[intersect(which(pval.vec < 0.05), which(sign.vec == 1))] <- 'UP'
dif.vec[intersect(which(pval.vec < 0.05), which(sign.vec == -1))] <- 'DOWN'
plot.dat <- data.frame('dif' = dif.vec,
                       'cluster' = factor(rep(colnames(clust.exp.table), 2), levels = colnames(clust.exp.table)),
                       'test' = rep(c('KO.v.WT', 'KOF4.v.KO'), each = ncol(clust.exp.table)))
ggplot(plot.dat, aes(x = test, y = cluster)) + geom_point(aes(color = dif), size = 3) +
  scale_color_manual(values = c('UP' = 'red', 'NS' = 'darkgrey', 'DOWN' = 'blue')) +
  labs(x = 'Condition Comparison', y = 'Cluster', title = 'Significance in Changes of Cluster Frequency by Experimental Condition') +
  guides(color = guide_legend(title = "Change"))
ggsave('revisions/anchored/anchored_cluster-exp-significant-changes.jpg', height = 10, width = 8, units = 'in', dpi = 300)
###############

### umap + bargraph plots
###############
anchored.obj <- readRDS('revisions/anchored/anchored_seur.rds')
## create plot data frame
ct.exp <- sapply(rownames(anchored.obj@reductions$umap@cell.embeddings), function(x) {strsplit(x, '_')[[1]][1]})
ct.exp <- factor(ct.exp, levels = c('WT', 'KO', 'KOF4'))
plot.df <- data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[,1],
                      'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[,2],
                      'Cluster' = anchored.obj$seurat_clusters,
                      'Exp' = ct.exp)
# cluster umap plot
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.75) +
  ggtitle('Anchored Analysis - Seurat Clusters')
ggsave('revisions/anchored/anchored_seurat-clust_umap.jpg', height = 8, width = 10, units = 'in', dpi = 300)
# exp umap plot
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Exp), size = 0.75) +
  ggtitle('Anchored Analysis - Experimental Condition') +
  guides(color = guide_legend(title="Condition"))
ggsave('revisions/anchored/anchored_seurat-exp_umap.jpg', height = 8, width = 10, units = 'in', dpi = 300)
# seperate exp umap plots
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.5) +
  ggtitle('Anchored Analysis - Experimental Condition') +
  facet_grid(cols = vars(Exp))
ggsave('revisions/anchored/anchored_seurat-clust-by-exp_umap.jpg', height = 6, width = 12, units = 'in', dpi = 300)
## percentage of experimental groups in each cluster
clust.exp.table <- table(ct.exp, anchored.obj$seurat_clusters)
barplot.df <- melt(clust.exp.table)
colnames(barplot.df) <- c('Exp', 'Cluster', 'Count')
barplot.df$Cluster <- as.factor(barplot.df$Cluster)
ggplot(barplot.df, aes(x = Cluster, y = Count, fill = Exp)) + geom_bar(position = 'fill', stat = "identity") +
  guides(fill = guide_legend(title="Condition")) +
  labs(y = "Proportion") + 
  ggtitle('Anchored Analysis - Cluster Condition Percentages')
ggsave('revisions/anchored/anchored_clust-exp-proportion_barplot.jpg', height = 8, width = 10, units = 'in', dpi = 300)
###############

### heatmaps
###############
anchored.obj <- readRDS('revisions/anchored/anchored_seur.rds')
ct.exp <- sapply(rownames(anchored.obj@reductions$umap@cell.embeddings), function(x) {strsplit(x, '_')[[1]][1]})
ct.exp <- factor(ct.exp, levels = c('WT', 'KO', 'KOF4'))
## create column annotations
cell.order <- names(sort(anchored.obj$seurat_clusters))
cluster.cols <- group_colors(length(unique(anchored.obj$seurat_clusters)))
names(cluster.cols) <- sort(unique(anchored.obj$seurat_clusters))
condition.cols <- group_colors(3); names(condition.cols) <- c('WT', 'KO', 'KOF4')
column.annot <- columnAnnotation('Cluster' = anchored.obj$seurat_clusters[cell.order],
                                 'Condition' = ct.exp[cell.order],
                                 col = list('Cluster' = cluster.cols,
                                            'Condition' = condition.cols))
col.gaps <- anchored.obj$seurat_clusters[cell.order]
## zarkada markers
plot.mat <- as.matrix(anchored.obj@assays$SCT@scale.data[, cell.order])
marker.use.inds <- which(zarkada.marker.df$Gene %in% rownames(plot.mat))
plot.mat <- plot.mat[zarkada.marker.df$Gene[marker.use.inds],]
# set row annotation and gaps
row.annot <- rowAnnotation('Marker' = zarkada.marker.df$Group[marker.use.inds],
                           col = list('Marker' = zarkada.group.col),
                           show_annotation_name = FALSE)
row.gaps <- zarkada.marker.df$Group[marker.use.inds]
# scale and set colors
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
# make plot
map.title <- 'Anchored Analysis - Zarkada Marker Expression'
jpeg('revisions/anchored/anchored_zarkada-marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        column_split = col.gaps, row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
## kuo markers
plot.mat <- as.matrix(anchored.obj@assays$SCT@scale.data[, cell.order])
marker.use.inds <- which(kuo.marker.df$Gene %in% rownames(plot.mat))
plot.mat <- plot.mat[kuo.marker.df$Gene[marker.use.inds],]
# set row annotation and gaps
row.annot <- rowAnnotation('Marker' = kuo.marker.df$Group[marker.use.inds],
                           col = list('Marker' = kuo.group.col),
                           show_annotation_name = FALSE)
row.gaps <- kuo.marker.df$Group[marker.use.inds]
# scale and set colors
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
# make plot
map.title <- 'Anchored Analysis - Kuo Marker Expression'
jpeg('revisions/anchored/anchored_kuo-marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        column_split = col.gaps, row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
## top markers
num.markers <- 10
marker.set <- lapply(names(cluster.markers), function(x) {
  cm.markers <- cluster.markers[[x]]
  c.markers <- cm.markers[order(cm.markers$avg_log2FC, decreasing = TRUE),]
  return(data.frame('gene' = rownames(c.markers)[1:num.markers],
                    'group' = rep(x, num.markers)))
})
marker.set <- Reduce(rbind, marker.set)
marker.set <- marker.set[which(marker.set$gene %in% rownames(anchored.obj@assays$SCT@scale.data)),]
marker.set$group <- as.numeric(marker.set$group)
plot.mat <- as.matrix(anchored.obj@assays$SCT@scale.data[marker.set$gene, cell.order])
# set row annotations and gaps
row.annot <- rowAnnotation('Cluster' = marker.set$group,
                           col = list('Cluster' = cluster.cols),
                           show_annotation_name = FALSE)
row.gaps <- marker.set$group
# scale and set colors
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
# make plot
map.title <- 'Anchored Analysis - Seurat Cluster Markers'
jpeg('revisions/anchored/anchored_cluster-markers_heatmap.jpg',
     height = 15, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        column_split = col.gaps, row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL,
        row_names_gp = gpar(fontsize = 8))
dev.off()
###############

### marker umaps
###############
anchored.obj <- readRDS('revisions/anchored/anchored_seur.rds')
sct.mat <- anchored.obj@assays$SCT@scale.data
ct.exp <- sapply(rownames(anchored.obj@reductions$umap@cell.embeddings), function(x) {strsplit(x, '_')[[1]][1]})
ct.exp <- factor(ct.exp, levels = c('WT', 'KO', 'KOF4'))
## zarkada markers
marker.use.inds <- which(zarkada.marker.df$Gene %in% rownames(sct.mat))
for (ind in marker.use.inds) {
  m.gene <- zarkada.marker.df$Gene[ind]
  m.group <- zarkada.marker.df$Group[ind]
  # create plot df
  plot.df <- data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[,1],
                        'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[,2],
                        'marker' = sct.mat[m.gene,],
                        'Exp' = ct.exp)
  # make plot
  ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = marker), size = 0.5) + 
    ggtitle(paste(m.gene, ' (', m.group, ')', sep = '')) +
    scale_color_gradient2(low = 'green', mid = 'darkgrey', high = 'purple', name = m.gene) +
    facet_grid(cols = vars(Exp))
  ggsave(paste('revisions/anchored/zarkada-marker_umap/', m.gene, '_', m.group, '_sct-umap.jpg', sep = ''),
         height = 6, width = 12, units = 'in', dpi = 300)
}
## kuo markers
marker.use.inds <- which(kuo.marker.df$Gene %in% rownames(sct.mat))
for (ind in marker.use.inds) {
  m.gene <- kuo.marker.df$Gene[ind]
  m.group <- kuo.marker.df$Group[ind]
  # create plot df
  plot.df <- data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[,1],
                        'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[,2],
                        'marker' = sct.mat[m.gene,])
  # make plot
  ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = marker)) + 
    ggtitle(paste(m.gene, ' (', m.group, ')', sep = '')) +
    scale_color_gradient2(low = 'green', mid = 'darkgrey', high = 'purple', name = m.gene)
  ggsave(paste('revisions/anchored/kuo-marker_umap/', m.gene, '_', m.group, '_sct-umap.jpg', sep = ''),
         height = 8, width = 10, units = 'in', dpi = 300)
}
## additional marker set
gene.set <- c('Angpt2', 'Mfsd2a', 'Bmx', 'Mki67', 'Spock2', 
              'Slc22a8 ', 'Selenop', 'Itm2a', 'Bsg', 'Cdh5', 
              'PECAM', 'Cdh5', 'Vwf', 'CD34')
gene.set <- intersect(gene.set, rownames(sct.mat))
plot.df.list <- lapply(gene.set, function(x) {
  return(data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[,1],
                    'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[,2],
                    'Expression' = sct.mat[x,],
                    'Marker' = x))
})
plot.df <- Reduce(rbind, plot.df.list)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Expression), size = 0.5) +
  facet_wrap(vars(Marker), nrow = 3) +
  scale_color_gradient2(low = 'green', mid = 'darkgrey', high = 'purple') + 
  ggtitle('Anchored Marker Expression') + 
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20))
ggsave('revisions/anchored/ct_plots/anchored_marker-expression-v2_umap.jpg', height = 12, width = 12, units = 'in', dpi = 300)
###############

### marker umaps 2.0
###############
anchored.obj <- readRDS('revisions/anchored/anchored_seur.rds')
sct.mat <- anchored.obj@assays$SCT@scale.data
ct.exp <- sapply(rownames(anchored.obj@reductions$umap@cell.embeddings), function(x) {strsplit(x, '_')[[1]][1]})
ct.exp <- factor(ct.exp, levels = c('WT', 'KO', 'KOF4'))
## zarkada differnetially expression genes - p6 and p10 groups
marker.use.inds <- which(dif.gene.df$Gene %in% rownames(sct.mat))
for (ind in marker.use.inds) {
  m.gene <- dif.gene.df$Gene[ind]
  print(m.gene)
  m.group <- dif.gene.df$Group[ind]
  # create plot df
  plot.df <- data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[,1],
                        'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[,2],
                        'marker' = sct.mat[m.gene,])
  # make plot
  ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = marker)) + 
    ggtitle(paste(m.gene, ' (', m.group, ')', sep = '')) +
    scale_color_gradient2(low = 'green', mid = 'darkgrey', high = 'purple', name = m.gene)
  ggsave(paste('revisions/anchored/zarkada-dif-exp_umap/', m.gene, '_', m.group, '_sct-umap.jpg', sep = ''),
         height = 8, width = 10, units = 'in', dpi = 300)
}
## d and s tip
marker.use.inds <- which(ds.dif.gene.df$Gene %in% rownames(sct.mat))
for (ind in marker.use.inds) {
  m.gene <- ds.dif.gene.df$Gene[ind]
  print(m.gene)
  m.group <- ds.dif.gene.df$Group[ind]
  # create plot df
  plot.df <- data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[,1],
                        'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[,2],
                        'marker' = sct.mat[m.gene,])
  # make plot
  ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = marker)) + 
    ggtitle(paste(m.gene, ' (', m.group, ')', sep = '')) +
    scale_color_gradient2(low = 'green', mid = 'darkgrey', high = 'purple', name = m.gene)
  ggsave(paste('revisions/anchored/zarkada-dif-exp_umap/', m.gene, '_', m.group, '_sct-umap.jpg', sep = ''),
         height = 8, width = 10, units = 'in', dpi = 300)
}
###############

### exp condition differential expression
###############
anchored.obj <- readRDS('revisions/anchored/anchored_seur.rds')
sct.mat <- anchored.obj@assays$SCT@scale.data
ct.exp <- sapply(rownames(anchored.obj@reductions$umap@cell.embeddings), function(x) {strsplit(x, '_')[[1]][1]})
ct.exp <- factor(ct.exp, levels = c('WT', 'KO', 'KOF4'))
## get groups
wt.cells <- names(ct.exp)[which(ct.exp == 'WT')]
ko.cells <- names(ct.exp)[which(ct.exp == 'KO')]
kof4.cells <- names(ct.exp)[which(ct.exp == 'KOF4')]
## set genes
gene.list <- c('Bmx', 'Mki67', 'Ptgis', 'Angpt2', 'Mfsd2a', 'Itm2a')
## ko v wt differential expression analysis
ko.v.wt.difexp <- sapply(gene.list, function(x) {
  test.obj <- wilcox.test(sct.mat[x, ko.cells], sct.mat[x, wt.cells], alternative = 'two.sided')
  p.val <- test.obj$p.val
  rbs.cor <- 2 * test.obj$statistic / (length(ko.cells) * length(wt.cells)) - 1
  return(data.frame('p.val' = p.val, 'rbs.cor' = rbs.cor,
                    'wt.mean' = mean(sct.mat[x, wt.cells]),
                    'ko.mean' = mean(sct.mat[x, ko.cells])))
})
ko.v.wt.difexp <- t(ko.v.wt.difexp)
write.csv(ko.v.wt.difexp, file = 'revisions/anchored/dif-exp_analysis/ko-v-wt_9-26-22-dif-exp.csv', row.names = TRUE, quote = FALSE)
## kof4 v ko differential expression analysis
kof4.v.ko.difexp <- sapply(gene.list, function(x) {
  test.obj <- wilcox.test(sct.mat[x, kof4.cells], sct.mat[x, ko.cells], alternative = 'two.sided')
  p.val <- test.obj$p.val
  rbs.cor <- 2 * test.obj$statistic / (length(kof4.cells) * length(ko.cells)) - 1
  return(data.frame('p.val' = p.val, 'rbs.cor' = rbs.cor,
                    'ko.mean' = mean(sct.mat[x, ko.cells]),
                    'kof4.mean' = mean(sct.mat[x, kof4.cells])))
})
kof4.v.ko.difexp <- t(kof4.v.ko.difexp)
write.csv(kof4.v.ko.difexp, file = 'revisions/anchored/dif-exp_analysis/kof4-v-ko_9-26-22-dif-exp.csv', row.names = TRUE, quote = FALSE)
###############

### cluster-cell type grouping
###############
anchored.obj <- readRDS('revisions/anchored/anchored_seur.rds')
sct.mat <- as.matrix(anchored.obj@assays$SCT@scale.data)
ct.exp <- sapply(rownames(anchored.obj@reductions$umap@cell.embeddings), function(x) {strsplit(x, '_')[[1]][1]})
ct.exp <- factor(ct.exp, levels = c('WT', 'KO', 'KOF4'))
## map from clusters to cell types
ct.vec <- anchored.obj$seurat_clusters
ct.vec <- plyr::mapvalues(ct.vec, from = 0:12,
                          to = c('Capillary', 'Capillary', 'Arterial', 'D-Tip', 'S-Tip',
                                 'Capillary', 'Proliferative', 'Proliferative', 'Venous', 'Venous',
                                 'Venous', 'Venous', 'Venous'))
## sample umaps
plot.df <- data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[,1],
                      'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[,2],
                      'Cell.Type' = ct.vec,
                      'Condition' = ct.exp)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cell.Type), size = 0.5) +
  ggtitle('Anchored Analysis - Experimental Condition') +
  guides(color = guide_legend(title = "Cell Type",
                              override.aes = list(size = 5))) +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20)) +
  facet_grid(cols = vars(Condition))
ggsave('revisions/anchored/ct_plots/anchored_ct-by-exp_umap.jpg', height = 6, width = 12, units = 'in', dpi = 300)
# merged umap
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cell.Type), size = 0.5) +
  ggtitle('Anchored Analysis - Cell Types') +
  guides(color = guide_legend(title = "Cell Type",
                              override.aes = list(size = 5))) +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20))
ggsave('revisions/anchored/ct_plots/anchored_ct_umap.jpg', height = 6, width = 8, units = 'in', dpi = 300)
## condition - cell type bar graph
ct.exp.table <- table(ct.exp, ct.vec)
barplot.df <- melt(ct.exp.table)
colnames(barplot.df) <- c('Exp', 'Cell.Type', 'Count')
barplot.df$Cell.Type <- as.factor(barplot.df$Cell.Type)
# change to percentage
total.vec <- rep(rowSums(ct.exp.table), ncol(ct.exp.table))
barplot.df$Percent <- barplot.df$Count / total.vec
# make plot
ggplot(barplot.df, aes(x = Cell.Type, y = Percent, fill = Exp)) + geom_bar(position = position_dodge(), stat = "identity") +
  guides(fill = guide_legend(title="Condition")) +
  labs(y = "Percentage", x = "Cell Type") + 
  scale_y_continuous(labels=scales::percent) +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20)) +
  ggtitle('Anchored Analysis - Cell Type Percentage by Condition')
ggsave('revisions/anchored/ct_plots/anchored_ct-exp-proportion_barplot.jpg', 
       height = 8, width = 10, units = 'in', dpi = 300)
## condition - cell type statistical analysis + dot plot
wt.v.ko <- chisq.test(x = t(ct.exp.table[c('WT', 'KO'),]))
ko.v.kof4 <- chisq.test(x = t(ct.exp.table[c('KO', 'KOF4'),]))
chisq.tests <- list('wt.v.ko' = wt.v.ko, 'ko.v.kof4' = ko.v.kof4)
saveRDS(chisq.tests, file = 'revisions/anchored/anchored_ct-chisq-tests.rds')
# make dot plot for cluster differences
wt.v.ko.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(wt.v.ko$stdres), lower.tail = FALSE)), ncol = 2)
ko.v.kof4.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(ko.v.kof4$stdres), lower.tail = FALSE)), ncol = 2)
pval.vec <- c(wt.v.ko.adjp[,2], ko.v.kof4.adjp[,2])
sign.vec <- c(sign(wt.v.ko$stdres[,2]), sign(ko.v.kof4$stdres[,2]))
dif.vec <- rep('NS', 2*ncol(ct.exp.table))
dif.vec[intersect(which(pval.vec < 0.01), which(sign.vec == 1))] <- 'UP'
dif.vec[intersect(which(pval.vec < 0.01), which(sign.vec == -1))] <- 'DOWN'
plot.dat <- data.frame('dif' = dif.vec,
                       'cluster' = factor(rep(colnames(ct.exp.table), 2), levels = colnames(ct.exp.table)),
                       'test' = rep(c('KO.v.WT', 'KOF4.v.KO'), each = ncol(ct.exp.table)))
ggplot(plot.dat, aes(x = test, y = cluster)) + geom_point(aes(color = dif), size = 3) +
  scale_color_manual(values = c('UP' = 'red', 'NS' = 'darkgrey', 'DOWN' = 'blue')) +
  labs(x = 'Condition Comparison', y = 'Cluster', 
       title = 'Significant Changes of Cell Type Frequency by Experimental Condition') +
  guides(color = guide_legend(title = "Change"))
ggsave('revisions/anchored/ct_plots/anchored_ct-exp-significant-changes.jpg', 
       height = 10, width = 8, units = 'in', dpi = 300)
# format to csv
chisq.df <- data.frame('KO.v.WT.pval' = wt.v.ko.adjp[,2],
                       'KO.v.WT.sign' = sign(wt.v.ko$stdres[,2]),
                       'KOF4.v.KO.pval' = ko.v.kof4.adjp[,2],
                       'KOF4.v.KO.sign' = sign(ko.v.kof4$stdres[,2]))
write.csv(chisq.df, file = 'revisions/anchored/anchored_ct-chisq-tests.csv', 
          row.names = TRUE, quote = FALSE)
## zarkada heatmap
cell.order <- names(sort(ct.vec))
ct.cols <- group_colors(length(unique(ct.vec)))
names(ct.cols) <- sort(unique(ct.vec))
condition.cols <- group_colors(3); names(condition.cols) <- c('WT', 'KO', 'KOF4')
column.annot <- columnAnnotation('Cell Type' = ct.vec[cell.order],
                                 'Condition' = ct.exp[cell.order],
                                 col = list('Cell Type' = ct.cols,
                                            'Condition' = condition.cols))
col.gaps <- ct.vec[cell.order]
# set plot matrix
plot.mat <- as.matrix(anchored.obj@assays$SCT@scale.data[, cell.order])
marker.use.inds <- which(zarkada.marker.df$Gene %in% rownames(plot.mat))
plot.mat <- plot.mat[zarkada.marker.df$Gene[marker.use.inds],]
# set row annotation and gaps
row.annot <- rowAnnotation('Marker' = zarkada.marker.df$Group[marker.use.inds],
                           col = list('Marker' = zarkada.group.col),
                           show_annotation_name = FALSE)
row.gaps <- zarkada.marker.df$Group[marker.use.inds]
# scale and set colors
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
# make plot
map.title <- 'Anchored Analysis - Zarkada Marker Expression'
jpeg('revisions/anchored/ct_plots/anchored_zarkada-ct-marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        column_split = col.gaps, row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
## zarkada gsea
z.gsea <- zarkada_gsea(as.matrix(anchored.obj@assays$SCT@scale.data), clust.vec = ct.vec,
                       plot.title = 'Anchored Analysis - Zarkada GSEA',
                       x.title = 'Cell Type Classification')
saveRDS(z.gsea, file = 'revisions/anchored/anchored_ct-zarkada-gsea.rds')
jpeg('revisions/anchored/ct_plots/anchored_ct-zarkada-gsea.jpg', height = 8, width = 10, units = 'in', res = 250)
print(z.gsea$plot + scale_size(range = c(5, 15)) + labs(x = 'Cell Type Classification'))
dev.off()
## maker umaps
gene.list <- c('Angpt2', 'Mfsd2a', 'Bmx', 'Mki67', 'Ptgis', 'Itm2a')
plot.df.list <- lapply(gene.list, function(x) {
  return(data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[,1],
                    'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[,2],
                    'Expression' = sct.mat[x,],
                    'Marker' = x))
})
plot.df <- Reduce(rbind, plot.df.list)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Expression), size = 0.5) +
  facet_wrap(vars(Marker), nrow = 3) +
  scale_color_gradient2(low = 'green', mid = 'darkgrey', high = 'purple') + 
  ggtitle('Anchored Marker Expression') + 
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20))
ggsave('revisions/anchored/ct_plots/anchored_marker-expression_umap.jpg', height = 12, width = 8, units = 'in', dpi = 300)
###############

###########################################################################
###########################################################################



###########################################################################
# Check previous viper trajectory analysis
###########################################################################

### recompute KOF4 viper activity due to corrupted file
###############
wt.seur <- readRDS('revisions/WT/wt_seur.rds')
kof4.seur <- readRDS('revisions/KOF4/KOF4_seur.rds')
wt.meta <- readRDS('sc/sc-nets/wt-k5-net_pruned.rds')
## create cpm matrices
wt.cpm <- cpm_norm(as.matrix(wt.seur@assays$SCT@counts))
kof4.cpm <- cpm_norm(as.matrix(kof4.seur@assays$SCT@counts))
shared.genes <- intersect(rownames(wt.cpm), rownames(kof4.cpm))
## create signature
wt.mean <- apply(wt.cpm[shared.genes,], 1, mean)
wt.sd <- apply(wt.cpm[shared.genes,], 1, sd)
wt.ges <- apply(wt.cpm[shared.genes,], 2, function(x) { (x - wt.mean) / wt.sd })
kof4.ges <- apply(kof4.cpm[shared.genes,], 2, function(x) { (x - wt.mean) / wt.sd })
## run viper
kof4.vip <- viper::viper(kof4.ges, wt.meta, method = 'none')
saveRDS(kof4.vip, file = 'revisions/kof4/kof4_vip.rds')
###############

### unified trajectory trajectory
###############
###############

### wt trained trajectory
###############
wt.mwk <- readRDS('sc/trajectory/wt_vip-mwk.rds')
wt.phi <- wt.mwk[[1]][1,]
wt.phi <- phi_adjust(wt.phi)
wt.model <- readRDS('sc/trajectory/wt_vip-mwk-rf.rds')
ko.phi <- readRDS('sc/trajectory/ko_vip-mwk_wt-pred.rds')
kof4.phi <- readRDS('sc/trajectory/kof4_vip-mwk_wt-pred.rds')
## load expression and viper data
wt.seur <- readRDS('revisions/WT/wt_seur.rds')
ko.seur <- readRDS('revisions/KO/ko_seur.rds')
kof4.seur <- readRDS('revisions/KOF4/kof4_seur.rds')
wt.vip <- readRDS('sc/WT/wt-internal_sc-vip.rds')
ko.vip <- readRDS('sc/KO/k0-v-wt_sc-vip.rds')
kof4.vip <- readRDS('revisions/kof4/kof4_vip.rds')
## fit gmm to wt
wt.fit <- Mclust(wt.phi, G = 3)
print(wt.fit$loglik)
## predict
wt.pred <- wt.fit$classification
ko.pred <- predict(wt.fit, ko.phi)$classification 
names(ko.pred) <- names(ko.phi)
kof4.pred <- predict(wt.fit, kof4.phi)$classification 
names(kof4.pred) <- names(kof4.phi)
wt.table <- table(wt.pred)
ko.table <- table(ko.pred)
kof4.table <- table(kof4.pred)
## plot
phi.vec <- c(ko.phi, kof4.phi, wt.phi)
exp.vec <- factor(c(rep('KO', length(ko.phi)), rep('KOF4', length(kof4.phi)), rep('WT', length(wt.phi))),
                  levels = c('KO', 'KOF4', 'WT'))
phi.dense.df <- data.frame('Phi' = phi.vec, 'Exp' = exp.vec)
ggplot(phi.dense.df, aes(x = Phi)) + geom_density(aes(fill = Exp), alpha = 0.4) +
  stat_function(fun = function(x) {dnorm(x, wt.fit$parameters$mean[1], sqrt(wt.fit$parameters$variance$sigmasq[1])) * 0.2},
                color = 'blue', linetype = 'dashed', size = 1) +
  stat_function(fun = function(x) {dnorm(x, wt.fit$parameters$mean[2], sqrt(wt.fit$parameters$variance$sigmasq[2])) * 0.2},
                color = 'purple', linetype = 'dashed', size = 1) +
  stat_function(fun = function(x) {dnorm(x, wt.fit$parameters$mean[3], sqrt(wt.fit$parameters$variance$sigmasq[3])) * 0.2},
                color = 'red', linetype = 'dashed', size = 1) +
  facet_wrap(vars(Exp), nrow = 3) +
  labs(title = "Phi Densities by Experiment", y = 'Density') +
  theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 20), 
        axis.title = element_text(hjust = 0.5, colour = "black", size = 16), 
        axis.text.x = element_text(hjust = 0.5, colour = "black", size = 12), 
        axis.text.y = element_text(hjust = 1, colour = "black", size = 12), 
        legend.title = element_text(hjust = 0, colour = "black", size = 16), 
        legend.text = element_text(hjust = 0, colour = "black", size = 12), legend.position = "right")
ggsave('revisions/wt_trajectory/phi-dense-plot.jpg',
       height = 6, width = 10, units = 'in', dpi = 250)
## stacked bar graph
plot.df <- data.frame('Count' = c(table(wt.pred), table(ko.pred), table(kof4.pred)),
                      'Exp' = factor(rep(c('WT', 'KO', 'KOF4'), each = 3), levels = c('WT', 'KO', 'KOF4')),
                      'Group' = factor(rep(c('Stalk', 'Intermediate', 'Tip'), 3), levels = c('Stalk', 'Intermediate', 'Tip')))
ggplot(plot.df, aes(fill = Group, y = Count, x = Exp)) +
  geom_bar(position = 'fill', stat = 'identity') +
  labs(x = 'Experimental Group', y = 'Percentage', title = 'WT Trajectory - Cell Groups')
ggsave('revisions/wt_trajectory/cell-group_bar-chart.jpg', height = 6, width = 10, units = 'in', dpi = 250)
## WT trajectory heatmaps
# zarkada marker expression
plot.mat <- as.matrix(wt.seur@assays$SCT@scale.data)
marker.use.inds <- which(zarkada.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('WT_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[zarkada.marker.df$Gene[marker.use.inds], intersect(names(sort(wt.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = wt.phi[colnames(plot.mat)],
                                 'Group' = factor(wt.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = zarkada.marker.df$Group[marker.use.inds],
                           col = list('Marker' = zarkada.group.col),
                           show_annotation_name = FALSE)
row.gaps <- zarkada.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
map.title <- 'WT - Phi Trajectory w/ Zarkada Markers'
jpeg('revisions/wt_trajectory/wt_zarkada-marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
# kuo expression
plot.mat <- as.matrix(wt.seur@assays$SCT@scale.data)
marker.use.inds <- which(kuo.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('WT_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[kuo.marker.df$Gene[marker.use.inds], intersect(names(sort(wt.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = wt.phi[colnames(plot.mat)],
                                 'Group' = factor(wt.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = kuo.marker.df$Group[marker.use.inds],
                           col = list('Marker' = kuo.group.col),
                           show_annotation_name = FALSE)
row.gaps <- kuo.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
map.title <- 'WT - Phi Trajectory w/ Kuo Markers'
jpeg('revisions/wt_trajectory/wt_kuo-marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
# zarkada activity
plot.mat <- wt.vip
marker.use.inds <- which(zarkada.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('WT_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[zarkada.marker.df$Gene[marker.use.inds], intersect(names(sort(wt.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = wt.phi[colnames(plot.mat)],
                                 'Group' = factor(wt.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = zarkada.marker.df$Group[marker.use.inds],
                           col = list('Marker' = zarkada.group.col),
                           show_annotation_name = FALSE)
row.gaps <- zarkada.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
map.title <- 'WT - Phi Trajectory w/ Zarkada Marker Activity'
jpeg('revisions/wt_trajectory/wt_zarkada-marker-activity_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'VIPER',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
# kuo activity
plot.mat <- wt.vip
marker.use.inds <- which(kuo.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('WT_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[kuo.marker.df$Gene[marker.use.inds], intersect(names(sort(wt.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = wt.phi[colnames(plot.mat)],
                                 'Group' = factor(wt.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = kuo.marker.df$Group[marker.use.inds],
                           col = list('Marker' = kuo.group.col),
                           show_annotation_name = FALSE)
row.gaps <- kuo.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
map.title <- 'WT - Phi Trajectory w/ Kuo Marker Activity'
jpeg('revisions/wt_trajectory/wt_kuo-marker-activity_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'VIPER',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
## KO trajectory heatmaps
# zarkada marker expression
plot.mat <- as.matrix(ko.seur@assays$SCT@scale.data)
marker.use.inds <- which(zarkada.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('KO_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[zarkada.marker.df$Gene[marker.use.inds], intersect(names(sort(ko.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = ko.phi[colnames(plot.mat)],
                                 'Group' = factor(ko.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = zarkada.marker.df$Group[marker.use.inds],
                           col = list('Marker' = zarkada.group.col),
                           show_annotation_name = FALSE)
row.gaps <- zarkada.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
map.title <- 'KO - Phi Trajectory w/ Zarkada Markers'
jpeg('revisions/wt_trajectory/ko_zarkada-marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
# kuo expression
plot.mat <- as.matrix(ko.seur@assays$SCT@scale.data)
marker.use.inds <- which(kuo.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('KO_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[kuo.marker.df$Gene[marker.use.inds], intersect(names(sort(ko.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = ko.phi[colnames(plot.mat)],
                                 'Group' = factor(ko.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = kuo.marker.df$Group[marker.use.inds],
                           col = list('Marker' = kuo.group.col),
                           show_annotation_name = FALSE)
row.gaps <- kuo.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
map.title <- 'KO - Phi Trajectory w/ Kuo Markers'
jpeg('revisions/wt_trajectory/ko_kuo-marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
# zarkada activity
plot.mat <- ko.vip
marker.use.inds <- which(zarkada.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('KO_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[zarkada.marker.df$Gene[marker.use.inds], intersect(names(sort(ko.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = ko.phi[colnames(plot.mat)],
                                 'Group' = factor(ko.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = zarkada.marker.df$Group[marker.use.inds],
                           col = list('Marker' = zarkada.group.col),
                           show_annotation_name = FALSE)
row.gaps <- zarkada.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
map.title <- 'KO - Phi Trajectory w/ Zarkada Marker Activity'
jpeg('revisions/wt_trajectory/ko_zarkada-marker-activity_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'VIPER',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
# kuo activity
plot.mat <- ko.vip
marker.use.inds <- which(kuo.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('KO_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[kuo.marker.df$Gene[marker.use.inds], intersect(names(sort(ko.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = ko.phi[colnames(plot.mat)],
                                 'Group' = factor(ko.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = kuo.marker.df$Group[marker.use.inds],
                           col = list('Marker' = kuo.group.col),
                           show_annotation_name = FALSE)
row.gaps <- kuo.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
map.title <- 'ko - Phi Trajectory w/ Kuo Marker Activity'
jpeg('revisions/wt_trajectory/ko_kuo-marker-activity_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'VIPER',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
## KOF4 trajectory heatmaps
# zarkada marker expression
plot.mat <- as.matrix(kof4.seur@assays$SCT@scale.data)
marker.use.inds <- which(zarkada.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('KOF4_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[zarkada.marker.df$Gene[marker.use.inds], intersect(names(sort(kof4.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = kof4.phi[colnames(plot.mat)],
                                 'Group' = factor(kof4.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = zarkada.marker.df$Group[marker.use.inds],
                           col = list('Marker' = zarkada.group.col),
                           show_annotation_name = FALSE)
row.gaps <- zarkada.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
map.title <- 'KOF4 - Phi Trajectory w/ Zarkada Markers'
jpeg('revisions/wt_trajectory/kof4_zarkada-marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
# kuo expression
plot.mat <- as.matrix(kof4.seur@assays$SCT@scale.data)
marker.use.inds <- which(kuo.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('KOF4_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[kuo.marker.df$Gene[marker.use.inds], intersect(names(sort(kof4.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = kof4.phi[colnames(plot.mat)],
                                 'Group' = factor(kof4.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = kuo.marker.df$Group[marker.use.inds],
                           col = list('Marker' = kuo.group.col),
                           show_annotation_name = FALSE)
row.gaps <- kuo.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
map.title <- 'KOF4 - Phi Trajectory w/ Kuo Markers'
jpeg('revisions/wt_trajectory/kof4_kuo-marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
# zarkada activity
plot.mat <- kof4.vip
marker.use.inds <- which(zarkada.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('KOF4_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[zarkada.marker.df$Gene[marker.use.inds], intersect(names(sort(kof4.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = kof4.phi[colnames(plot.mat)],
                                 'Group' = factor(kof4.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = zarkada.marker.df$Group[marker.use.inds],
                           col = list('Marker' = zarkada.group.col),
                           show_annotation_name = FALSE)
row.gaps <- zarkada.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
map.title <- 'KOF4 - Phi Trajectory w/ Zarkada Marker Activity'
jpeg('revisions/wt_trajectory/kof4_zarkada-marker-activity_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'VIPER',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
# kuo activity
plot.mat <- kof4.vip
marker.use.inds <- which(kuo.marker.df$Gene %in% rownames(plot.mat))
colnames(plot.mat) <- paste('KOF4_', colnames(plot.mat), sep = '')
plot.mat <- plot.mat[kuo.marker.df$Gene[marker.use.inds], intersect(names(sort(kof4.phi)), colnames(plot.mat))]
column.annot <- columnAnnotation('Phi' = kof4.phi[colnames(plot.mat)],
                                 'Group' = factor(kof4.pred[colnames(plot.mat)], levels = c(1, 2, 3), labels = c('Stalk', 'Intermediate', 'Tip')),
                                 col = list('Phi' = phi.col.func,
                                            'Group' = cell.group.cols))
row.annot <- rowAnnotation('Marker' = kuo.marker.df$Group[marker.use.inds],
                           col = list('Marker' = kuo.group.col),
                           show_annotation_name = FALSE)
row.gaps <- kuo.marker.df$Group[marker.use.inds]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
map.title <- 'KOF4 - Phi Trajectory w/ Kuo Marker Activity'
jpeg('revisions/wt_trajectory/kof4_kuo-marker-activity_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'VIPER',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
###############

###########################################################################
###########################################################################



###########################################################################
# New NaRnEA Analysis
###########################################################################

### process networks
###############
wt.tf <- readRDS('revisions/narnea_analysis/wt_net/wt_tf_regulon-list.rds')
wt.surf <- readRDS('revisions/narnea_analysis/wt_net/wt_surf_regulon-list.rds')
wt.list <- c(wt.tf, wt.surf)
saveRDS(wt.list, file = 'revisions/narnea_analysis/wt_net/wt_network.rds')
###############

### run NaRnEA w/ WT as reference
###############
wt.net <- readRDS('revisions/narnea_analysis/wt_net/wt_network.rds')
wt.cpm <- readRDS('revisions/WT/wt_cpm-counts.rds')
ko.cpm <- readRDS('revisions/KO/ko_cpm-counts.rds')
kof4.cpm <- readRDS('revisions/KOF4/kof4_cpm-counts.rds')
## get reference data
shared.genes <- intersect(rownames(wt.cpm), intersect(rownames(ko.cpm), rownames(kof4.cpm)))
wt.mean <- rowMeans(wt.cpm[shared.genes,])
wt.sd <- apply(wt.cpm[shared.genes,], 1, sd)
## create signatures
wt.ges <- (wt.cpm[shared.genes,] - wt.mean) / wt.sd
ko.ges <- (ko.cpm[shared.genes,] - wt.mean) / wt.sd
kof4.ges <- (kof4.cpm[shared.genes,] - wt.mean) / wt.sd
sig.list <- list('wt' = wt.ges, 'ko' = ko.ges, 'kof4' = kof4.ges)
saveRDS(sig.list, file = 'revisions/narnea_analysis/v-wt_ges-list.rds')
## run narnea
wt.narnea <- matrix_narnea(wt.ges, wt.net)
ko.narnea <- matrix_narnea(ko.ges, wt.net)
kof4.narnea <- matrix_narnea(kof4.ges, wt.net)
narnea.list <- list('wt' = wt.narnea, 'ko' = ko.narnea, 'kof4' = kof4.narnea)
saveRDS(narnea.list, file = 'revisions/narnea_analysis/v-wt_narnea-list.rds')
###############

### standard clustering analysis
###############
wt.cpm <- readRDS('revisions/WT/wt_cpm-counts.rds')
ko.cpm <- readRDS('revisions/KO/ko_cpm-counts.rds')
kof4.cpm <- readRDS('revisions/KOF4/kof4_cpm-counts.rds')
narnea.list <- readRDS('revisions/narnea_analysis/v-wt_narnea-list.rds')
## wt clustering
narnea.name <- 'wt'
narnea.obj <- narnea.list$wt
narnea.pca <- fast_pca(narnea.obj$PES)
narnea.dist <- dist(narnea.pca$x)
narnea.clust <- louvain_k(narnea.dist)
saveRDS(narnea.dist, file = paste('revisions/narnea_analysis/', narnea.name, '/', 
                                  narnea.name, '_narnea-dist.rds', sep = ''))
saveRDS(narnea.clust, file = paste('revisions/narnea_analysis/', narnea.name, '/', 
                                   narnea.name, '_narnea-clust.rds', sep = ''))
# MR analysis
narnea.mr <- kw_cluster_mrs(narnea.obj, clust.vec = narnea.clust$opt.clust)
saveRDS(narnea.mr, file = paste('revisions/narnea_analysis/', narnea.name, '/',
                                narnea.name, '_narnea-mrs.rds', sep = ''))
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_narnea-mr-heatmap_regs.jpg', sep = ''),
     height = 12, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact',  
                   clust.vec = narnea.clust$opt.clust,
                   mr.list = narnea.mr, reg.class = 'regulator', 
                   plot.title = paste(toupper(narnea.name), ': NaRnEA MRs', sep = ''))
dev.off()
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_narnea-mr-heatmap_markers.jpg', sep = ''),
     height = 12, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact',  
                   clust.vec = narnea.clust$opt.clust,
                   mr.list = narnea.mr, reg.class = 'marker', 
                   plot.title = paste(toupper(narnea.name), ': NaRnEA MRs', sep = ''))
dev.off()
# zarkada and kuo marker heatmaps
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_zarkada-marker-pact.jpg', sep = ''),
     height = 8, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust,
                   marker.set = zarkada.marker.df,
                   plot.title = paste(toupper(narnea.name), ': Zarkada Markers', sep = ''))
dev.off()
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_kuo-marker-pact.jpg', sep = ''),
     height = 8, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust,
                   marker.set = kuo.marker.df,
                   plot.title = paste(toupper(narnea.name), ': Zarkada Markers', sep = ''))
dev.off()
# zarkada gsea
z.gsea <- zarkada_gsea(wt.cpm, narnea.clust$opt.clust,
                       plot.title = paste(toupper(narnea.name), ': Zarkada GSEA', sep = ''))
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', 
           narnea.name, '_zarkada-gsea-dot-plot.jpg', sep = ''),
     height = 8, width = 10, units = 'in', res = 250)
print(z.gsea$plot)
dev.off()
## ko clustering
narnea.name <- 'ko'
narnea.obj <- narnea.list$ko
narnea.pca <- fast_pca(narnea.obj$PES)
narnea.dist <- dist(narnea.pca$x)
narnea.clust <- louvain_k(narnea.dist)
saveRDS(narnea.dist, file = paste('revisions/narnea_analysis/', narnea.name, '/', 
                                  narnea.name, '_narnea-dist.rds', sep = ''))
saveRDS(narnea.clust, file = paste('revisions/narnea_analysis/', narnea.name, '/', 
                                   narnea.name, '_narnea-clust.rds', sep = ''))
# MR analysis
narnea.mr <- kw_cluster_mrs(narnea.obj, clust.vec = narnea.clust$opt.clust)
saveRDS(narnea.mr, file = paste('revisions/narnea_analysis/', narnea.name, '/',
                                narnea.name, '_narnea-mrs.rds', sep = ''))
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_narnea-mr-heatmap_regs.jpg', sep = ''),
     height = 12, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact',  
                   clust.vec = narnea.clust$opt.clust,
                   mr.list = narnea.mr, reg.class = 'regulator', 
                   plot.title = paste(toupper(narnea.name), ': NaRnEA MRs', sep = ''))
dev.off()
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_narnea-mr-heatmap_markers.jpg', sep = ''),
     height = 12, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact',  
                   clust.vec = narnea.clust$opt.clust,
                   mr.list = narnea.mr, reg.class = 'marker', 
                   plot.title = paste(toupper(narnea.name), ': NaRnEA MRs', sep = ''))
dev.off()
# zarkada and kuo marker heatmaps
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_zarkada-marker-pact.jpg', sep = ''),
     height = 8, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust,
                   marker.set = zarkada.marker.df,
                   plot.title = paste(toupper(narnea.name), ': Zarkada Markers', sep = ''))
dev.off()
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_kuo-marker-pact.jpg', sep = ''),
     height = 8, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust,
                   marker.set = kuo.marker.df,
                   plot.title = paste(toupper(narnea.name), ': Zarkada Markers', sep = ''))
dev.off()
# zarkada gsea
z.gsea <- zarkada_gsea(ko.cpm, narnea.clust$opt.clust,
                       plot.title = paste(toupper(narnea.name), ': Zarkada GSEA', sep = ''))
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', 
           narnea.name, '_zarkada-gsea-dot-plot.jpg', sep = ''),
     height = 8, width = 10, units = 'in', res = 250)
print(z.gsea$plot)
dev.off()
## kof4 clustering
narnea.name <- 'kof4'
narnea.obj <- narnea.list$kof4
narnea.pca <- fast_pca(narnea.obj$PES)
narnea.dist <- dist(narnea.pca$x)
narnea.clust <- louvain_k(narnea.dist)
saveRDS(narnea.dist, file = paste('revisions/narnea_analysis/', narnea.name, '/', 
                                  narnea.name, '_narnea-dist.rds', sep = ''))
saveRDS(narnea.clust, file = paste('revisions/narnea_analysis/', narnea.name, '/', 
                                   narnea.name, '_narnea-clust.rds', sep = ''))
# MR analysis
narnea.mr <- kw_cluster_mrs(narnea.obj, clust.vec = narnea.clust$opt.clust)
saveRDS(narnea.mr, file = paste('revisions/narnea_analysis/', narnea.name, '/',
                                narnea.name, '_narnea-mrs.rds', sep = ''))
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_narnea-mr-heatmap_regs.jpg', sep = ''),
     height = 12, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact',  
                   clust.vec = narnea.clust$opt.clust,
                   mr.list = narnea.mr, reg.class = 'regulator', 
                   plot.title = paste(toupper(narnea.name), ': NaRnEA MRs', sep = ''))
dev.off()
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_narnea-mr-heatmap_markers.jpg', sep = ''),
     height = 12, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact',  
                   clust.vec = narnea.clust$opt.clust,
                   mr.list = narnea.mr, reg.class = 'marker', 
                   plot.title = paste(toupper(narnea.name), ': NaRnEA MRs', sep = ''))
dev.off()
# zarkada and kuo marker heatmaps
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_zarkada-marker-pact.jpg', sep = ''),
     height = 8, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust,
                   marker.set = zarkada.marker.df,
                   plot.title = paste(toupper(narnea.name), ': Zarkada Markers', sep = ''))
dev.off()
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', narnea.name, '_kuo-marker-pact.jpg', sep = ''),
     height = 8, width = 10, units = 'in', res = 250)
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust,
                   marker.set = kuo.marker.df,
                   plot.title = paste(toupper(narnea.name), ': Zarkada Markers', sep = ''))
dev.off()
# zarkada gsea
z.gsea <- zarkada_gsea(kof4.cpm, narnea.clust$opt.clust,
                       plot.title = paste(toupper(narnea.name), ': Zarkada GSEA', sep = ''))
jpeg(paste('revisions/narnea_analysis/', narnea.name, '/', 
           narnea.name, '_zarkada-gsea-dot-plot.jpg', sep = ''),
     height = 8, width = 10, units = 'in', res = 250)
print(z.gsea$plot)
dev.off()
###############

### generate WT trajectory and train results for KO and KOF4
###############
narnea.list <- readRDS('revisions/narnea_analysis/v-wt_narnea-list.rds')
## generate clustering
wt.dist <- readRDS('revisions/narnea_analysis/wt/wt_narnea-dist.rds')
wt.k2 <- pam_k(wt.dist, kmin = 2, kmax = 3)
clust.vec <- wt.k2$opt.clust
k2.centers <- cbind(rowMeans(narnea.list$wt$PES[, which(clust.vec == 2)]),
                    rowMeans(narnea.list$wt$PES[, which(clust.vec == 1)]))
wt.k2.mwk <- MWKMeans(narnea.list$wt$PES, k2.centers)
saveRDS(clust.vec, file = 'revisions/narnea_analysis/wt/wt_narnea-k2.rds')
saveRDS(wt.k2.mwk, file = 'revisions/narnea_analysis/wt/wt_narnea-mwk.rds')
## wt plots
phi.vec <- wt.k2.mwk[[1]][2,]
phi.vec <- phi_adjust(phi.vec)
cell.order <- names(sort(phi.vec))
## column annotation
column.annot <- columnAnnotation('Phi' = phi.vec[cell.order],
                                 col = list('Phi' = phi.col.func))
## create kuo plot
use.inds <- which(kuo.marker.df$Gene %in% rownames(narnea.list$wt$PES))
row.annot <- rowAnnotation('Marker' = kuo.marker.df$Group[use.inds],
                           col = list('Marker' = kuo.group.col),
                           show_annotation_name = FALSE)
row.gaps <- kuo.marker.df$Group[use.inds]
plot.mat <- narnea.list$wt$PES[, cell.order]
plot.mat <- plot.mat[kuo.marker.df$Gene[use.inds],]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
map.title <- 'WT: NaRnEA MWK Analysis w/ Kuo Markers'
jpeg(paste('revisions/narnea_analysis/wt/wt_mwk_kuo-marker-heatmap.jpg'),
     height = 8, width = 10, units = 'in', res = 250)
Heatmap(plot.mat, name = 'PES',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
## create zarkada plot
use.inds <- which(zarkada.marker.df$Gene %in% rownames(narnea.list$wt$PES))
row.annot <- rowAnnotation('Marker' = zarkada.marker.df$Group[use.inds],
                           col = list('Marker' = zarkada.group.col),
                           show_annotation_name = FALSE)
row.gaps <- zarkada.marker.df$Group[use.inds]
plot.mat <- narnea.list$wt$PES[, cell.order]
plot.mat <- plot.mat[zarkada.marker.df$Gene[use.inds],]
plot.mat <- t(apply(plot.mat, 1, scale))
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('pact', length(col.breaks)))
map.title <- 'WT: NaRnEA MWK Analysis w/ Zarkada Markers'
jpeg(paste('revisions/narnea_analysis/wt/wt_mwk_zarkada-marker-heatmap.jpg'),
     height = 8, width = 10, units = 'in', res = 250)
Heatmap(plot.mat, name = 'PES',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = map.title, row_title = NULL)
dev.off()
###############

