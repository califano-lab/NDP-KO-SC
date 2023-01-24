setwd('C://Users/Lukas Vlahos/OneDrive/ac_lab/kuo-colab/')
library(Seurat)
library(PISCES)
library(NaRnEA)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

### marker sets and color palettes
###############
marker.set <- list('D-tip' = c('Cldn5', 'Mfsd2a', 'Spock2'),
                   'S-tip' = c('Angpt2', 'Kcn3', 'Esm1'),
                   'Arterial' = c('Bmx', 'Unc5b'),
                   'Venous' = c('Hmgn2', 'Ptgis'),
                   'Capillary' = c('Slc22a8'),
                   'Proliferative' = c('Mki67', 'Birc5'),
                   'Pan-Endothelial' = c('Cdh5', 'Pecam1', 'Vwf', 'CD34'))
marker.df <- data.frame('Gene' = unlist(marker.set),
                        'Group' = rep(names(marker.set), sapply(marker.set, length)))
## wnt and proliferative set
wnt.proliferative.markers <- c('Mki67', 'Bicr5', 'Axin2', 'Plvap', 'Apcdd1', 'Cldn5')
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
## zarkada dge sets
zarkada.p6.dge <- read.csv(file = 'revisions/zarkada_p6-genes.csv')
colnames(zarkada.p6.dge) <- c('Tip', 'Arterial', 'Capillary', 'Proliferative', 'Venous')
zarkada.p6.regulons <- apply(zarkada.p6.dge, 2, function(x) {
  gene.set <- x[which(x != '')]
  gene.set <- unique(gene.set)
  am.vec <- rep(1, length(gene.set))
  names(am.vec) <- gene.set
  aw.vec <- rep(1, length(gene.set))
  names(aw.vec) <- gene.set
  return(list('am' = am.vec, 'aw' = aw.vec))
})
names(zarkada.p6.regulons) <- colnames(zarkada.p6.dge)
zarkada.p10.dge <- read.csv(file = 'revisions/zarkada_p10-genes.csv')
colnames(zarkada.p10.dge) <- c('Tip', 'Arterial', 'Capillary', 'Proliferative', 'Venous')
zarkada.p10.regulons <- apply(zarkada.p10.dge, 2, function(x) {
  gene.set <- x[which(x != '')]
  gene.set <- unique(gene.set)
  am.vec <- rep(1, length(gene.set))
  names(am.vec) <- gene.set
  aw.vec <- rep(1, length(gene.set))
  names(aw.vec) <- gene.set
  return(list('am' = am.vec, 'aw' = aw.vec))
})
names(zarkada.p10.regulons) <- colnames(zarkada.p10.dge)
## colors
clust.colors <- group_colors(13); names(clust.colors) <- 0:12
condition.colors <- group_colors(3); names(condition.colors) <- c('WT', 'NdpKO', 'NdpKO + L6-F4-2')
marker.colors <- group_colors(length(unique(marker.df$Group))); names(marker.colors) <- sort(unique(marker.df$Group))
ct.names <- sort(unique(marker.df$Group))
ct.colors <- group_colors(length(ct.names)); names(ct.colors) <- ct.names
## plot theme
plot.theme <- theme(text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20))
###############

#' Zarkada GSEA function
zarkada_gsea <- function(counts.mat, clust.vec, plot.title, x.title, regulon.list) {
  
  ## create signature and run narnea
  sig.vecs <- list()
  narnea.vecs <- list()
  for (cn in sort(unique(clust.vec))) {
    print(cn)
    # calculate signature
    cat('Calculating signature...\n')
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
    cat('Running NaRnEA...\n')
    cn.narnea <- sapply(regulon.list, function(reg.obj) {
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
    return(unlist(x['pes',]))
  })
  row.has.na <- apply(pes.mat, 1, function(x) {length(which(is.na(x))) > 0})
  pes.mat <- pes.mat[!row.has.na,]
  nes.mat <- sapply(narnea.vecs, function(x) {
    return(unlist(x['nes',]))
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

### cluster analysis
###############
anchored.obj <- readRDS('revisions/anchored/anchored_seur.rds')
sample.treatment <- sapply(colnames(anchored.obj@assays$SCT@scale.data), function(x) {
  strsplit(x, '_')[[1]][1]
})
sample.treatment <- factor(sample.treatment, levels = c('WT', 'KO', 'KOF4'))
sample.treatment.name <- plyr::mapvalues(sample.treatment, from = c('WT', 'KO', 'KOF4'),
                                         to = c('WT', 'NdpKO', 'NdpKO + L6-F4-2'))
## cluster umaps
plot.df <- data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[,1],
                      'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[,2],
                      'Condition' = sample.treatment.name,
                      'Cluster' = anchored.obj$seurat_clusters)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Condition), size = 0.5) +
  ggtitle('Anchored Analysis - Experimental Condition') +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  plot.theme
ggsave('revisions/anchored/v2_results/cluster_experimental-condition-umap.jpg',
       height = 6, width = 8, units = 'in', dpi = 300)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.5) +
  ggtitle('Anchored Analysis - Seurat Clusters') +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  plot.theme
ggsave('revisions/anchored/v2_results/cluster_seurat-cluster-umap.jpg',
       height = 6, width = 8, units = 'in', dpi = 300)
# condition specific UMAP
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster), size = 0.5) +
  ggtitle('Anchored Analysis - Clusters by Experimental Condition') +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  plot.theme + 
  facet_grid(cols = vars(Condition))
ggsave('revisions/anchored/v2_results/cluster_seurat-clust-by-condition.jpg', 
       height = 6, width = 12, units = 'in', dpi = 300)
## visualization of CDH5 for cluster exclusion
plot.df <- data.frame('Expression' = anchored.obj@assays$SCT@scale.data['Cdh5',],
                      'Gene' = 'Cdh5',
                      'Cluster' = anchored.obj$seurat_clusters)
ggplot(plot.df, aes(x = Cluster, y = Expression)) + 
  geom_violin(aes(group = Cluster, color = Cluster, fill = Cluster)) + 
  labs(title = 'SCT Normalized Cdh5 expression across clusters',
       x = 'Seurat Cluster', y = 'Cdh5')
ggsave('revisions/anchored/v2_results/cdh5_cluster-expression.jpg',
       height = 6, width = 8, units = 'in', dpi = 300)
## one-sided Cdh5 p-value for each cluster
cdh5.vec <- anchored.obj@assays$SCT@scale.data['Cdh5',]
cdh5.df <- sapply(sort(unique(anchored.obj$seurat_clusters)), function(x) {
  test.samps <- which(anchored.obj$seurat_clusters == x)
  ref.samps <- which(anchored.obj$seurat_clusters != x)
  test.obj <- wilcox.test(cdh5.vec[test.samps], cdh5.vec[ref.samps],
                          alternative = 'less')
  rbs.cor <- 2 * test.obj$statistic / (length(test.samps) * length(ref.samps)) - 1
  return(c(mean(cdh5.vec[test.samps]), test.obj$p.val, rbs.cor))
})
cdh5.df <- t(cdh5.df)
colnames(cdh5.df) <- c('Mean', 'p.val', 'RBSC')
rownames(cdh5.df) <- sort(unique(anchored.obj$seurat_clusters))
write.csv(cdh5.df, 'revisions/anchored/v2_results/cdh5_cluster-expression_wilcox-test.csv', 
          row.names = TRUE, quote = FALSE)
## visualization of Pecam1 for cluster exclusion
plot.df <- data.frame('Expression' = anchored.obj@assays$SCT@scale.data['Pecam1',],
                      'Gene' = 'Pecam1',
                      'Cluster' = anchored.obj$seurat_clusters)
ggplot(plot.df, aes(x = Cluster, y = Expression)) + 
  geom_violin(aes(group = Cluster, color = Cluster, fill = Cluster)) + 
  labs(title = 'SCT Normalized Pecam1 expression across clusters',
       x = 'Seurat Cluster', y = 'Pecam1')
ggsave('revisions/anchored/v2_results/pecam1_cluster-expression.jpg',
       height = 6, width = 8, units = 'in', dpi = 300)
## one-sided Pecam1 p-value for each cluster
pecam1.vec <- anchored.obj@assays$SCT@scale.data['Pecam1',]
pecam1.df <- sapply(sort(unique(anchored.obj$seurat_clusters)), function(x) {
  test.samps <- which(anchored.obj$seurat_clusters == x)
  ref.samps <- which(anchored.obj$seurat_clusters != x)
  test.obj <- wilcox.test(pecam1.vec[test.samps], pecam1.vec[ref.samps],
                          alternative = 'less')
  rbs.cor <- 2 * test.obj$statistic / (length(test.samps) * length(ref.samps)) - 1
  return(c(mean(pecam1.vec[test.samps]), test.obj$p.val, rbs.cor))
})
pecam1.df <- t(pecam1.df)
colnames(pecam1.df) <- c('Mean', 'p.val', 'RBSC')
rownames(pecam1.df) <- sort(unique(anchored.obj$seurat_clusters))
write.csv(pecam1.df, 'revisions/anchored/v2_results/pecam1_cluster-expression_wilcox-test.csv', 
          row.names = TRUE, quote = FALSE)
## set clusters and samples for use
use.clusters <- c(0, 1, 2, 3, 4, 5, 6, 7, 9)
use.samps <- names(anchored.obj$seurat_clusters)[which(anchored.obj$seurat_clusters %in% use.clusters)]
## marker heatmap
clust.vec <- anchored.obj$seurat_clusters[use.samps]
cell.order <- names(sort(clust.vec))
# column annotation
column.annot <- columnAnnotation('Cluster' = clust.vec[cell.order],
                                 'Condition' = sample.treatment.name[cell.order],
                                 col = list('Condition' = condition.colors,
                                            'Cluster' = clust.colors[sort(unique(clust.vec))]))
col.gaps <- clust.vec[cell.order]
# plot matrix
plot.mat <- as.matrix(anchored.obj@assays$SCT@scale.data[, cell.order])
marker.use.inds <- which(marker.df$Gene %in% rownames(plot.mat))
plot.mat <- plot.mat[marker.df$Gene[marker.use.inds],]
# row annotation
row.annot <- rowAnnotation('Marker' = marker.df$Group[marker.use.inds],
                           col = list('Marker' = marker.colors),
                           show_annotation_name = FALSE)
row.gaps <- marker.df$Group[marker.use.inds]
# plot colors and title
plot.mat <- t(apply(plot.mat, 1, scale))
colnames(plot.mat) <- cell.order
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
plot.title <- 'Seurat Cluster - Marker Expression'
# generate plot
jpeg('revisions/anchored/v2_results/cluster_marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        column_split = col.gaps, row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_dend_reorder = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = plot.title, row_title = NULL)
dev.off()
## zarkada gsea
z.gsea <- zarkada_gsea(as.matrix(anchored.obj@assays$SCT@scale.data)[, use.samps], 
                       clust.vec = anchored.obj$seurat_clusters[use.samps],
                       plot.title = 'Seurat Clusters - Zarkada GSEA',
                       x.title = 'Seurat Cluster')
saveRDS(z.gsea, file = 'revisions/anchored/v2_results/cluster_zarkada-gsea.rds')
jpeg('revisions/anchored/v2_results/cluster_zarkada-gsea.jpg', 
     height = 8, width = 10, units = 'in', res = 250)
print(z.gsea$plot + scale_size(range = c(5, 15)))
dev.off()
## p6 cluster dge
z.p6.gsea <- zarkada_gsea(as.matrix(anchored.obj@assays$SCT@scale.data)[, use.samps],
                          clust.vec = anchored.obj$seurat_clusters[use.samps],
                          plot.title = 'Seurat Clusters - Zarkada P6 DGE Enrichment',
                          x.title = 'Seruat Clusters',
                          regulon.list = zarkada.p6.regulons)
saveRDS(z.p6.gsea, file = 'revisions/anchored/v2_results/cluster_zarkada-p6-dge-gsea.rds')
jpeg('revisions/anchored/v2_results/cluster_zarkada-p6-dge-gsea.jpg', 
     height = 8, width = 10, units = 'in', res = 250)
print(z.p6.gsea$plot + scale_size(range = c(5, 15)))
dev.off()
## p6 cluster dge
z.p10.gsea <- zarkada_gsea(as.matrix(anchored.obj@assays$SCT@scale.data)[, use.samps],
                           clust.vec = anchored.obj$seurat_clusters[use.samps],
                           plot.title = 'Seurat Clusters - Zarkada P10 DGE Enrichment',
                           x.title = 'Seruat Clusters',
                           regulon.list = zarkada.p10.regulons)
saveRDS(z.p10.gsea, file = 'revisions/anchored/v2_results/cluster_zarkada-p10-dge-gsea.rds')
jpeg('revisions/anchored/v2_results/cluster_zarkada-p10-dge-gsea.jpg', 
     height = 8, width = 10, units = 'in', res = 250)
print(z.p10.gsea$plot + scale_size(range = c(5, 15)))
dev.off()
## proliferative and wnt marker violin plot
sct.mat <- as.matrix(anchored.obj@assays$SCT@scale.data[, use.samps])
gene.list <- intersect(rownames(sct.mat), wnt.proliferative.markers)
plot.df.list <- lapply(gene.list, function(x) {
  return(data.frame('Expression' = sct.mat[x, use.samps],
                    'Gene' = x,
                    'Condition' = sample.treatment.name[use.samps],
                    'Cluster' = anchored.obj$seurat_clusters[use.samps]))
})
plot.df <- Reduce(rbind, plot.df.list)
ggplot(plot.df, aes(x = Condition, y = Expression)) +
  geom_violin(aes(color = Cluster)) +
  facet_wrap(vars(Gene), nrow = 3) +
  ggtitle('Marker Expression across Treatments and Clusters') + 
  plot.theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('revisions/anchored/v2_results/cluster_wnt-proliferative-marker-expression_violin.jpg', 
       height = 8, width = 10, units = 'in', dpi = 300)
## find cluster markers
sub.obj <- anchored.obj[, use.samps]
clust.vec <- sub.obj$seurat_clusters
clust.markers <- list()
for (cv in sort(unique(clust.vec))) {
  print(cv)
  cn.markers <- FindMarkers(sub.obj, ident.1 = cv)
  cn.markers <- cn.markers[order(cn.markers$p_val_adj, -abs(cn.markers$avg_log2FC)),]
  clust.markers[[as.character(cv)]] <- cn.markers
}
saveRDS(clust.markers, file = 'revisions/anchored/v2_results/cluster_top-genes.rds')
## make csv with top 50 genes
top.markers <- lapply(clust.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = TRUE),]
  return(data.frame('gene' = rownames(c.markers)[1:50],
                    'log2FC' = c.markers$avg_log2FC[1:50]))
})
top.marker.df <- Reduce(cbind, top.markers)
colnames(top.marker.df) <- paste(rep(paste('c', names(top.markers), sep = ''), each = 2),
                                 colnames(top.marker.df), sep = '.')
write.table(top.marker.df, file = 'revisions/anchored/v2_results/cluster_top-genes.csv',
            sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
## create heatmap with top 10 genes
clust.vec <- anchored.obj$seurat_clusters[use.samps]
cell.order <- names(sort(clust.vec))
marker.df <- sapply(clust.markers, function(x) {
  c.markers <- x[order(x$avg_log2FC, decreasing = TRUE),]
  return(rownames(c.markers)[1:10])
})
marker.df <- as.data.frame(Reduce(rbind, marker.df))
marker.df$Cluster <- rep(names(clust.markers), each = 10)
colnames(marker.df) <- c('Gene', 'Cluster')
# row annotation
row.annot <- rowAnnotation('Cluster' = marker.df$Cluster,
                           col = list('Cluster' = clust.colors), 
                           show_annotation_name = FALSE,
                           show_legend = TRUE)
row.gaps <- marker.df$Cluster
# create plot data
plot.mat <- as.matrix(anchored.obj@assays$SCT@scale.data[, cell.order])
plot.mat <- plot.mat[marker.df$Gene,]
#plot.mat <- t(apply(plot.mat, 1, scale))
colnames(plot.mat) <- cell.order
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
plot.title <- 'Seurat Cluster - Top Genes'
# make plot
jpeg('revisions/anchored/v2_results/cluster_top-genes.jpg',
     height = 10, width = 15, units = 'in', res = 300)
print(Heatmap(plot.mat, name = 'GEXP',
              col = col.fun,
              top_annotation = column.annot, column_split = col.gaps,
              left_annotation = row.annot, row_split = row.gaps,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = FALSE,
              column_title = plot.title, row_title = NULL,
              row_names_gp = gpar(fontsize = 8)))
dev.off()
###############

### cell type analysis
###############
anchored.obj <- readRDS('revisions/anchored/anchored_seur.rds')
sample.treatment <- sapply(colnames(anchored.obj@assays$SCT@scale.data), function(x) {
  strsplit(x, '_')[[1]][1]
})
sample.treatment <- factor(sample.treatment, levels = c('WT', 'KO', 'KOF4'))
sample.treatment.name <- plyr::mapvalues(sample.treatment, from = c('WT', 'KO', 'KOF4'),
                                         to = c('WT', 'NdpKO', 'NdpKO + L6-F4-2'))
use.clusters <- c(0, 1, 2, 3, 4, 5, 6, 7, 9)
use.samps <- names(anchored.obj$seurat_clusters)[which(anchored.obj$seurat_clusters %in% use.clusters)]
sct.mat <- as.matrix(anchored.obj@assays$SCT@scale.data[, use.samps])
## map clusters to cell types
sample.ct <- anchored.obj$seurat_clusters[use.samps]
sample.ct <- plyr::mapvalues(sample.ct, from = 0:12,
                             to = c('Pan-Endothelial', 'Pan-Endothelial', 'Arterial', 'D-tip', 'S-tip',
                                    'Pan-Endothelial', 'Proliferative', 'Proliferative', 'Exclude', 'Venous',
                                    'Exclude', 'Exclude', 'Exclude'))
sample.ct <- factor(sample.ct, levels = c('Arterial', 'D-tip', 'Pan-Endothelial',
                                          'Proliferative', 'S-tip', 'Venous'))
## ct umaps
plot.df <- data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[use.samps,1],
                      'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[use.samps,2],
                      'Condition' = sample.treatment.name[use.samps],
                      'Cell.Type' = sample.ct[use.samps])
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cell.Type), size = 0.5) +
  ggtitle('Anchored Analysis - Cell Type') +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  plot.theme
ggsave('revisions/anchored/v2_results/ct_cell-types-umap.jpg',
       height = 6, width = 8, units = 'in', dpi = 300)
# condition specific UMAP
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cell.Type), size = 0.5) +
  ggtitle('Anchored Analysis - Cell Type by Experimental Condition') +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  plot.theme + 
  facet_grid(cols = vars(Condition))
ggsave('revisions/anchored/v2_results/ct_cell-type-by-condition.jpg', 
       height = 6, width = 18, units = 'in', dpi = 300)
## heatmap of markers
cell.order <- names(sort(sample.ct))
# column annotation
column.annot <- columnAnnotation('Cell Type' = sample.ct[cell.order],
                                 'Condition' = sample.treatment.name[cell.order],
                                 col = list('Condition' = condition.colors,
                                            'Cell Type' = ct.colors))
col.gaps <- sample.ct[cell.order]
# plot matrix
plot.mat <- as.matrix(anchored.obj@assays$SCT@scale.data[, cell.order])
marker.use.inds <- which(marker.df$Gene %in% rownames(plot.mat))
plot.mat <- plot.mat[marker.df$Gene[marker.use.inds],]
# row annotation
row.annot <- rowAnnotation('Marker' = marker.df$Group[marker.use.inds],
                           col = list('Marker' = marker.colors),
                           show_annotation_name = FALSE)
row.gaps <- marker.df$Group[marker.use.inds]
# plot colors and title
plot.mat <- t(apply(plot.mat, 1, scale))
colnames(plot.mat) <- cell.order
col.breaks <- quantile_breaks(plot.mat)
col.fun <- colorRamp2(col.breaks, color_levels('gexp', length(col.breaks)))
plot.title <- 'Cell Type Classification - Marker Expression'
# generate plot
jpeg('revisions/anchored/v2_results/ct_marker-expression_heatmap.jpg',
     height = 8, width = 8, units = 'in', res = 300)
Heatmap(plot.mat, name = 'Expression',
        col = col.fun,
        top_annotation = column.annot, left_annotation = row.annot, 
        column_split = col.gaps, row_split = row.gaps,
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_dend_reorder = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_title = plot.title, row_title = NULL)
dev.off()
## zarkada gsea
ct.z.gsea <- zarkada_gsea(as.matrix(anchored.obj@assays$SCT@scale.data)[, use.samps],
                          clust.vec = sample.ct,
                          plot.title = 'Cell Type Clasification - Zarkada GSEA',
                          x.title = 'Cell Type',
                          regulon.list = zarkada.regulons)
saveRDS(ct.z.gsea, file = 'revisions/anchored/v2_results/ct_zarkada-gsea.rds')
jpeg('revisions/anchored/v2_results/ct_zarkada-gsea.jpg', 
     height = 8, width = 10, units = 'in', res = 250)
print(ct.z.gsea$plot + scale_size(range = c(5, 15)))
dev.off()
## marker umaps
gene.list <- intersect(rownames(sct.mat), marker.df$Gene)
plot.df.list <- lapply(gene.list, function(x) {
  return(data.frame('UMAP1' = anchored.obj@reductions$umap@cell.embeddings[use.samps,1],
                    'UMAP2' = anchored.obj@reductions$umap@cell.embeddings[use.samps,2],
                    'Expression' = sct.mat[x,use.samps],
                    'Marker' = paste(x, ' (', marker.df$Group[which(marker.df$Gene == x)], ')', sep = '')))
})
plot.df <- Reduce(rbind, plot.df.list)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Expression), size = 0.5) +
  facet_wrap(vars(Marker), ncol = 3) +
  scale_color_gradient2(low = 'green', mid = 'darkgrey', high = 'purple') + 
  ggtitle('Anchored Marker Expression') + 
  plot.theme
ggsave('revisions/anchored/v2_results/ct_marker-expression_umap.jpg', 
       height = 12, width = 10, units = 'in', dpi = 300)
## marker violin plots
gene.list <- intersect(rownames(sct.mat), marker.df$Gene)
plot.df.list <- lapply(gene.list, function(x) {
  return(data.frame('Expression' = sct.mat[x,use.samps],
                    'Gene' = x,
                    'Condition' = sample.treatment.name[use.samps],
                    'Cell.Type' = sample.ct[use.samps]))
})
plot.df <- Reduce(rbind, plot.df.list)
ggplot(plot.df, aes(x = Condition, y = Expression)) +
  geom_violin(aes(color = Cell.Type)) +
  facet_wrap(vars(Gene), nrow = 3) +
  ggtitle('Marker Expression across Treatments and Cell Types') + 
  plot.theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('revisions/anchored/v2_results/ct_marker-expression_violin.jpg', 
       height = 8, width = 10, units = 'in', dpi = 300)
## wnt and proliferative marker violin plot
gene.list <- intersect(rownames(sct.mat), wnt.proliferative.markers)
plot.df.list <- lapply(gene.list, function(x) {
  return(data.frame('Expression' = sct.mat[x, use.samps],
                    'Gene' = x,
                    'Condition' = sample.treatment.name[use.samps],
                    'Cell.Type' = sample.ct[use.samps]))
})
plot.df <- Reduce(rbind, plot.df.list)
ggplot(plot.df, aes(x = Condition, y = Expression)) +
  geom_violin(aes(color = Cell.Type)) +
  facet_wrap(vars(Gene), nrow = 3) +
  ggtitle('Marker Expression across Treatments and Cell Types') + 
  plot.theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('revisions/anchored/v2_results/ct_wnt-proliferative-marker-expression_violin.jpg', 
       height = 8, width = 10, units = 'in', dpi = 300)
## bar graph of ct frequency
con.ct.table <- table(sample.treatment[use.samps], sample.ct)
con.ct.table <- con.ct.table[,c(1:5, 7)]
barplot.df <- melt(con.ct.table)
colnames(barplot.df) <- c('Exp', 'Cell.Type', 'Count')
barplot.df$Cell.Type <- as.factor(barplot.df$Cell.Type)
total.vec <- rep(rowSums(con.ct.table), ncol(con.ct.table))
barplot.df$Percent <- barplot.df$Count / total.vec
barplot.df$Exp <- plyr::mapvalues(barplot.df$Exp,
                                  from = c('WT', 'KO', 'KOF4'),
                                  to = c('WT', 'NdpKO', 'NdpKO + L6-F4-2'))
# make plot
ggplot(barplot.df, aes(x = Cell.Type, y = Percent, fill = Exp)) + 
  geom_bar(position = position_dodge(), stat = "identity") +
  guides(fill = guide_legend(title="Condition")) +
  labs(y = "Percentage", x = "Cell Type") + 
  scale_y_continuous(labels=scales::percent) +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20)) +
  ggtitle('Anchored Analysis - Cell Type Percentage by Condition')
ggsave('revisions/anchored/v2_results/ct_clust-frequency-bargraph.jpg', 
       height = 8, width = 10, units = 'in', dpi = 300)
## p-values + dot plot for ct frequency comparisons
wt.v.ko <- chisq.test(x = t(con.ct.table[c('WT', 'KO'),]))
ko.v.kof4 <- chisq.test(x = t(con.ct.table[c('KO', 'KOF4'),]))
chisq.tests <- list('wt.v.ko' = wt.v.ko, 'ko.v.kof4' = ko.v.kof4)
saveRDS(chisq.tests, file = 'revisions/anchored/v2_results/ct_clust-frequency-chisq.rds')
# create plot
wt.v.ko.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(wt.v.ko$stdres), lower.tail = FALSE)), ncol = 2)
ko.v.kof4.adjp <- matrix(p.adjust(method = 'BH', p = pnorm(abs(ko.v.kof4$stdres), lower.tail = FALSE)), ncol = 2)
pval.vec <- c(wt.v.ko.adjp[,2], ko.v.kof4.adjp[,2])
sign.vec <- c(sign(wt.v.ko$stdres[,2]), sign(ko.v.kof4$stdres[,2]))
dif.vec <- rep('NS', 2*ncol(con.ct.table))
dif.vec[intersect(which(pval.vec < 0.01), which(sign.vec == 1))] <- 'UP'
dif.vec[intersect(which(pval.vec < 0.01), which(sign.vec == -1))] <- 'DOWN'
plot.dat <- data.frame('dif' = dif.vec,
                       'cluster' = factor(rep(colnames(con.ct.table), 2), levels = colnames(con.ct.table)),
                       'test' = rep(c('KO.v.WT', 'KOF4.v.KO'), each = ncol(con.ct.table)))
ggplot(plot.dat, aes(x = test, y = cluster)) + geom_point(aes(color = dif), size = 5) +
  scale_color_manual(values = c('UP' = 'red', 'NS' = 'darkgrey', 'DOWN' = 'blue')) +
  labs(x = 'Condition Comparison', y = 'Cluster', 
       title = 'Significant Changes of Cell Type Frequency') +
  guides(color = guide_legend(title = "Change")) +
  plot.theme
ggsave('revisions/anchored/v2_results/ct_clust-frequency-pval-dot-plot.jpg', 
       height = 10, width = 8, units = 'in', dpi = 300)
# save as csv
chisq.df <- data.frame('KO.v.WT.pval' = wt.v.ko.adjp[,2],
                       'KO.v.WT.sign' = sign(wt.v.ko$stdres[,2]),
                       'KOF4.v.KO.pval' = ko.v.kof4.adjp[,2],
                       'KOF4.v.KO.sign' = sign(ko.v.kof4$stdres[,2]))
write.csv(chisq.df, file = 'revisions/anchored/v2_results/ct_clust-frequency-chisq.csv', 
          row.names = TRUE, quote = FALSE)
###############



clust.name <- 0
test.samps <- which(clust.vec == clust.name)
ref.samps <- which(clust.vec != clust.name)
clust.dge <- apply(sct.mat, 1, function(x) {
  w.test <- wilcox.test(x[test.samps], x[ref.samps], alternative = "two.sided")
  rbs.cor <- 2 * w.test$statistic / (length(test.samps) * length(ref.samps)) - 1
  return(list('p.val' = w.test$p.val, 'rbsc' = rbs.cor))
})
clust.dge.df <- as.data.frame(Reduce(rbind, clust.dge))
rownames(clust.dge.df) <- names(clust.dge)
clust.dge.df$p.val <- p.adjust(clust.dge.df$p.val, method = 'BH')
clust.dge.df <- clust.dge.df[which(clust.dge.df$p.val < 0.05),]
clust.dge.df <- clust.dge.df[order(unlist(clust.dge.df$rbsc), decreasing = TRUE),]

