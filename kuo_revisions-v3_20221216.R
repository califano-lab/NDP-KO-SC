setwd('C://Users/Lukas Vlahos/OneDrive/ac_lab/kuo-colab/')
library(Seurat)
library(ggplot2)
library(NaRnEA)


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
    return(sapply(x, function(y) {y$pes}))
  })
  row.has.na <- apply(pes.mat, 1, function(x) {length(which(is.na(x))) > 0})
  pes.mat <- pes.mat[!row.has.na,]
  nes.mat <- sapply(narnea.vecs, function(x) {
    return(sapply(x, function(y) {y$nes}))
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
    guides(size = guide_legend(title = "|NES|", order = 2),
           color = guide_legend(order = 1)) 
  
  ## return
  return(list('narnea' = narnea.vecs, 'plot' = plot.obj))
}

### figure 3E for each condition
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
## gsea for cts within each condition
for (cond.name in names(table(sample.treatment))) {
  print(cond.name)
  
  cond.cells <- intersect(use.samps, names(sample.treatment)[which(sample.treatment == cond.name)])
  cond.mat <- sct.mat[, cond.cells]
  print(dim(cond.mat))
  
  cond.gsea <- zarkada_gsea(cond.mat, clust.vec = sample.ct[cond.cells],
                            plot.title = paste(cond.name, ': Cell Type Classification - Zarkada GSEA', sep = ''),
                            x.title = 'Cell Type', regulon.list = zarkada.regulons)
  saveRDS(cond.gsea, file = paste('revisions-2/', cond.name, '_ct-zarkada-sea.rds', sep = ''))
  
  jpeg(paste('revisions-2/', cond.name , '_ct-zarkada-gsea.jpg', sep = ''), 
       height = 8, width = 10, units = 'in', res = 250)
  print(cond.gsea$plot + scale_size(range = c(5, 15)))
  dev.off()
}
## gsea between conditions
cond.gsea <- zarkada_gsea(sct.mat, clust.vec = sample.treatment[use.samps],
                          plot.title = 'Experimental Ccondition - Zarkada GSEA',
                          x.title = 'Experimental Condition',
                          regulon.list = zarkada.regulons)
saveRDS(cond.gsea, file = 'revisions-2/condition_zarkada-gsea.rds')

jpeg('revisions-2/condition_zarkada-gsea.jpg', 
     height = 8, width = 10, units = 'in', res = 250)
print(cond.gsea$plot + scale_size(range = c(5, 15)))
dev.off()
###############



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