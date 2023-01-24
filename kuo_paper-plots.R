setwd('C://Users/lvlah/linux/ac_lab/kuo-colab/')
library(PISCES)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(plot3D)
library(RColorBrewer)
library(mclust)
library(circlize)
library(ComplexHeatmap)
marker.list <- list('tip' = c('Plaur', 'Angpt2', 'Lcp2', 'Cxcr4', 'Apln', 
                              'Kcne3', 'Dll4', 'Esm1', 'Vegfr2', 'Vegfr3', 'Kdr', 'Flk1', 'Flt4'),
                    'stalk' = c('Tie2', 'Tek', 'Jag1', 'Apj', 'Vegfr1', 'Flt1'),
                    'general' = c('Fzd4', 'Plvap'))
phi_adjust <- function(x) {
  return(x**2 / (x**2 + (1 - x)**2))
}

### load data and set plot params
###############
anchor.gexp <- readRDS('sc/anchored/anchor_gexp.rds')
vwt.vip <- readRDS('sc/KO/k0-k0f4-v-wt_vip.rds')
anch.vip <- readRDS('sc/anchored/anchor_vip.rds')
anch.umap <- readRDS('sc/anchored/anchor_vip-umap.rds')
anch.pca <- readRDS('sc/anchored/anchor_vip-pca.rds')
anchor.mwk <- readRDS('sc/anchored/anchor_vip-k2-mwk.rds')
mr.df <- readRDS('sc/anchored/anchor_ko-kof4-mr-table.rds')
phi.vec <- phi_adjust(anchor.mwk[[1]][1,])
exp.vec <- as.factor(sapply(colnames(anch.vip), function(x) { strsplit(x, '_')[[1]][1] }))
marker.set <- intersect(unlist(marker.list), rownames(anch.vip))
phi.corr <- readRDS('sc/anchored/anchor_vip-k2-mwk-phi-corr.rds')
wt.fit <- readRDS('sc/anchored/anchor_vip-k2-mwk_wt-gmm.rds')
## ggplot theme
plot.theme <- theme(text = element_text(family = "sans"),
                    plot.title = element_text(hjust = 0.5, colour = "black", size = 40), 
                    axis.title = element_text(hjust = 0.5, colour = "black", size = 32), 
                    axis.text.x = element_text(hjust = 0.5, colour = "black", size = 24), 
                    axis.text.y = element_text(hjust = 1, colour = "black", size = 24), 
                    legend.title = element_text(hjust = 0, colour = "black", size = 32), 
                    legend.text = element_text(hjust = 0, colour = "black", size = 24),
                    strip.text.x = element_text(hjust = 0.5, colour = "black", size = 32),
                    strip.text.y = element_text(hjust = 0.5, colour = "black", size = 32),
                    legend.position = "right"
                    )
y.umap.lim <- c(-5, 5)
x.umap.lim <- c(-8, 8)
exp.colors <- ClusterColors(3); names(exp.colors) <- c('WT', 'KO', 'KOF4')
legend_params <- list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 12))
###############

### fig 1b
#################
cell.order <- names(sort(phi.vec))
num.breaks <- 100
plot.mat <- anch.vip[marker.set, cell.order]
col.func <- colorRamp2(QuantileBreaks(plot.mat, num.breaks), ColorLevels(num.breaks, 'vip'))
# make annotations
legend_params <- list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 16))
m.source <- rep(c('Tip', 'Stalk'), c(length(which(marker.set %in% marker.list$tip)),
                                     length(which(marker.set %in% marker.list$stalk))))
names(m.source) <- marker.set
row_color <- c('darkgreen', 'purple'); names(row_color) <- c('Tip', 'Stalk')
row_annot <- rowAnnotation('Marker' = m.source, col = list('Marker' = row_color),
                           show_annotation_name = FALSE,
                           annotation_legend_param = legend_params)
col_annot <- columnAnnotation('Phi' = phi.vec[cell.order],
                              'Experiment' = factor(exp.vec[cell.order], levels = names(exp.colors)),
                              col = list('Phi' = colorRamp2(c(0, 0.5, 1), c('white', 'grey', 'black')),
                                         'Experiment' = exp.colors),
                              annotation_legend_param = legend_params,
                              annotation_name_gp = gpar(fontsize =))
# make plot
jpeg('paper_components/fig-1_b.jpg', height = 1000, width = 750)
Heatmap(plot.mat, name = 'Activity',
        column_title = 'Marker Activity', 
        column_title_gp = grid::gpar(fontsize = 36),
        row_names_gp = grid::gpar(fontsize = 24),
        col = col.func, 
        left_annotation = row_annot, top_annotation = col_annot,
        heatmap_legend_param = legend_params,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE)
dev.off()
#################

### fig 1c
###############
plot.df <- data.frame('UMAP1' = anch.umap[,1], 'UMAP2' = anch.umap[,2],
                      'Phi' = phi.vec)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Phi)) + 
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  ggtitle('MWKMeans: Phi Values') +
  plot.theme + xlim(x.umap.lim) + ylim(y.umap.lim)
ggsave(filename = 'paper_components/fig-1_c.jpg', dpi = 500, height = 10, width = 12)
###############

### fig 1d
###############
cell.order <- names(sort(exp.vec))
exp.label <- factor(sort(exp.vec), levels = c('WT', 'KO', 'KOF4'))
## umap version
plot.df <- data.frame('UMAP1' = anch.umap[cell.order, 1], 'UMAP2' = anch.umap[cell.order, 2],
                      'Experiment' = exp.label)
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Experiment)) + 
  geom_density_2d(color = 'black') +
  facet_wrap(vars(Experiment), ncol = 3) +
  ggtitle('Experimental Condition-Specific Densities') +
  plot.theme + xlim(x.umap.lim) + ylim(y.umap.lim)
ggsave(filename = 'paper_components/fig-1_d.jpg', dpi = 500, height = 10, width = 32)
## pca version
plot.df <- data.frame('PC1' = anch.pca$x[cell.order, 1], 'PC2' = anch.pca$x[cell.order, 2],
                      'Experiment' = exp.label)
ggplot(plot.df, aes(PC1, PC2)) + geom_point(aes(color = Experiment)) + 
  geom_density_2d(color = 'black') +
  facet_wrap(vars(Experiment), ncol = 3) +
  ggtitle('Experimental Condition-Specific Densities') +
  plot.theme
ggsave(filename = 'paper_components/fig-1_d_pca.jpg', dpi = 500, height = 10, width = 32)
###############

### fig 1e
###############
cell.order <- names(sort(exp.vec))
exp.label <- factor(sort(exp.vec), levels = c('WT', 'KO', 'KOF4'))
plot.df <- data.frame('Phi' = phi.vec[cell.order], 'Experiment' = exp.label)
ggplot(plot.df, aes(x = Phi)) + geom_density(aes(fill = Experiment), alpha = 0.6) +
  facet_wrap(vars(Experiment), nrow = 3) +  
  stat_function(fun = function(x) {dnorm(x, wt.fit$parameters$mean[1], sqrt(wt.fit$parameters$variance$sigmasq[1])) * 0.2},
                color = 'blue', linetype = 'dashed', size = 1) +
  stat_function(fun = function(x) {dnorm(x, wt.fit$parameters$mean[2], sqrt(wt.fit$parameters$variance$sigmasq[2])) * 0.2},
                color = 'purple', linetype = 'dashed', size = 1) +
  stat_function(fun = function(x) {dnorm(x, wt.fit$parameters$mean[3], sqrt(wt.fit$parameters$variance$sigmasq[3])) * 0.2},
                color = 'red', linetype = 'dashed', size = 1) +
  ggtitle('Phi Densities by Experiment') +
  plot.theme
ggsave(filename = 'paper_components/fig-1_e.jpg', dpi = 500, height = 10, width = 12)
###############

### fig 1A
###############
num.mrs <- 50
exp.mr.sets <- apply(vwt.vip, 2, function(x) {
  x <- sort(x)
  mr.vec <- c(names(head(x, num.mrs / 2)), names(tail(x, num.mrs / 2)))
  return(mr.vec)
})
mr.set <- unique(unlist(as.list(exp.mr.sets)))
mr.label <- rep('KO', length(mr.set)); names(mr.label) <- mr.set
mr.label[which(mr.set %in% exp.mr.sets[,'K0F4'])] <- 'KOF4'
mr.label[which(mr.set %in% intersect(exp.mr.sets[, 'K0'], exp.mr.sets[, 'K0F4']))] <- 'Shared'
mr.label <- sort(mr.label, decreasing = FALSE)
mr.set <- names(mr.label)
## make plot
plot.mat <- vwt.vip[mr.set,]; colnames(plot.mat) <- c('KO', 'KOF4')
num.breaks <- 25
col.func <- colorRamp2(QuantileBreaks(plot.mat, num.breaks), ColorLevels(num.breaks, 'vip'))
# legend paramaters
legend_params <- list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 12))
# make annotation
source.col <- c(exp.colors[2:3], 'grey'); names(source.col) <- c('KO', 'KOF4', 'Shared')
row_annot <- rowAnnotation('Source' = mr.label, col = list('Source' = source.col), 
                           show_annotation_name = FALSE,
                           annotation_legend_param = legend_params)
col_annot <- columnAnnotation(foo = anno_text(c('KO', 'KOF4'), rot = 0, just = 'center',
                                              gp = gpar(fontsize = 16)))
# make plot
jpeg('paper_components/fig-1_a.jpg', height = 1000, width = 500)
Heatmap(plot.mat, name = 'Activity',
        column_title = 'Master Regulator Overlap', 
        column_title_gp = grid::gpar(fontsize = 28),
        row_names_gp = grid::gpar(fontsize = 12),
        col = col.func, 
        row_split = mr.label,
        row_title_gp = grid::gpar(fontsize = 16), row_title_rot = 0,
        left_annotation = row_annot,
        top_annotation = col_annot,
        heatmap_legend_param = legend_params,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE)
dev.off()
## transposed version
source.col <- c(exp.colors[2:3], 'grey'); names(source.col) <- c('KO', 'KOF4', 'Shared')
t_col_annot <- columnAnnotation('Source' = mr.label, col = list('Source' = source.col), 
                           show_annotation_name = FALSE,
                           annotation_legend_param = legend_params)
t_row_annot <- rowAnnotation(foo = anno_text(c('KO', 'KOF4'), rot = 0, just = 'left',
                                              gp = gpar(fontsize = 16)))
t.plot.mat <- t(plot.mat)
jpeg('paper_components/fig-1_a_transpose.jpg', height = 500, width = 1200)
Heatmap(t.plot.mat, name = 'Activity',
        #column_title = 'Master Regulator Overlap', 
        column_title_gp = grid::gpar(fontsize = 16),
        row_title_gp = grid::gpar(fontsize = 16), row_title_rot = 0,
        column_names_gp = grid::gpar(fontsize = 12),
        col = col.func, 
        column_split = mr.label,
        heatmap_legend_param = legend_params,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = TRUE,
        left_annotation = t_row_annot,
        top_annotation = t_col_annot)
dev.off()
###############

### fig s2A/B/C
###############
exp.label <- factor(sort(exp.vec), levels = c('WT', 'KO', 'KOF4'))
marker.df <- data.frame('UMAP1' = anch.umap[,1], 'UMAP2' = anch.umap[,2],
                        'Experiment' = exp.label)
# add protein activity
marker.set <- intersect(unlist(marker.list), rownames(anch.vip))
marker.activity <- unlist(as.list(t(anch.vip[marker.set,])))
marker.label <- factor(rep(marker.set, each = nrow(anch.umap)), levels = rev(unique(marker.set)))
marker.df <- do.call('rbind', replicate(length(marker.set), marker.df, simplify = FALSE))
marker.df$Activity <- marker.activity
marker.df$Protein <- marker.label
## stalk marker plot
stalk.marker.df <- marker.df[which(marker.df$Protein %in% marker.list$stalk),]
ggplot(stalk.marker.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Activity)) +  
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  facet_grid(vars(Experiment), vars(Protein)) +
  ggtitle('Stalk Marker Activity') +
  plot.theme + xlim(x.umap.lim) + ylim(y.umap.lim)
ggsave(filename = 'paper_components/fig-s2_a.jpg', dpi = 500, height = 17, width = 20)
## stalk marker plot
tip.marker.df <- marker.df[which(marker.df$Protein %in% marker.list$tip),]
ggplot(tip.marker.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Activity)) +  
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  facet_grid(vars(Experiment), vars(Protein)) +
  ggtitle('Tip Marker Activity') +
  plot.theme + xlim(x.umap.lim) + ylim(y.umap.lim)
ggsave(filename = 'paper_components/fig-s2_b.jpg', dpi = 500, height = 17, width = 30)
## stalk marker plot
general.marker.df <- marker.df[which(marker.df$Protein %in% marker.list$general),]
ggplot(general.marker.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Activity)) +  
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  facet_grid(vars(Experiment), vars(Protein)) +
  ggtitle('Marker Activity') +
  plot.theme + xlim(x.umap.lim) + ylim(y.umap.lim)
ggsave(filename = 'paper_components/fig-s2_c.jpg', dpi = 500, height = 17, width = 15)
###############

### fig s3 - gene expression UMAP plots
###############
plot.dat <- data.frame('UMAP1' = anchor.gexp@reductions$umap@cell.embeddings[,1],
                       'UMAP2' = anchor.gexp@reductions$umap@cell.embeddings[,2],
                       'Seurat' = anchor.gexp@active.ident)
exp.vec <- sapply(rownames(plot.dat), function(x) {strsplit(x, '_')[[1]][1]} )
exp.vec <- factor(exp.vec, levels = c('WT', 'KO', 'KOF4'))
plot.dat$Experiment <- exp.vec
plot.dat$Phi <- phi.vec[rownames(plot.dat)]
## add markers
sct.mat <- anchor.gexp@assays$integrated@scale.data
marker.set <- intersect(unlist(marker.list), rownames(sct.mat))
marker.expression <- unlist(as.list(t(sct.mat[marker.set,])))
marker.label <- factor(rep(marker.set, each = nrow(plot.dat)), levels = rev(unique(marker.set)))
marker.df <- do.call('rbind', replicate(length(marker.set), plot.dat[, c('UMAP1', 'UMAP2', 'Experiment')],
                                       simplify = FALSE))
marker.df$Expression <- marker.expression
marker.df$Gene <- marker.label
## get fzd4 expression
fzd4.vec <- anchor.gexp@assays$SCT@counts['Fzd4',]
fzd4.vec <- (fzd4.vec / sum(fzd4.vec)) * 1e6
fzd4.vec <- (fzd4.vec - mean(fzd4.vec)) / sd(fzd4.vec)
fzd4.df <- plot.dat[, c('UMAP1', 'UMAP2', 'Experiment')]
fzd4.df$Expression <- fzd4.vec; fzd4.df$Gene <- rep('Fzd4', length(fzd4.vec))
marker.df <- rbind(marker.df, fzd4.df)
## experimental label plot
ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = Experiment)) +
  facet_wrap(vars(Experiment), ncol = 3) +
  geom_density_2d(color = 'black') +
  ggtitle('Experimental Conditions') +
  plot.theme
ggsave(filename = 'paper_components/fig-s3_a.jpg', dpi = 500, height = 10, width = 32)
## phi values
ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = Phi), size = 2) +
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  facet_wrap(vars(Experiment), ncol = 3) +
  ggtitle('Phi Values') +
  plot.theme
ggsave(filename = 'paper_components/fig-s3_b.jpg', dpi = 500, height = 10, width = 32)
## tip marker plot
tip.marker.df <- marker.df[which(marker.df$Gene %in% marker.list$tip),]
ggplot(tip.marker.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Expression)) +  
  scale_color_gradient2(low = 'green', mid = 'grey', midpoint = 0.5, high = 'purple') +
  facet_grid(vars(Experiment), vars(Gene)) +
  ggtitle('Tip Marker Expression') +
  plot.theme
ggsave(filename = 'paper_components/fig-s3_c.jpg', dpi = 500, height = 17, width = 40)
## stalk marker plot
stalk.marker.df <- marker.df[which(marker.df$Gene %in% marker.list$stalk),]
ggplot(stalk.marker.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Expression)) +  
  scale_color_gradient2(low = 'green', mid = 'grey', midpoint = 0.5, high = 'purple') +
  facet_grid(vars(Experiment), vars(Gene)) +
  ggtitle('Stalk Marker Expression') +
  plot.theme
ggsave(filename = 'paper_components/fig-s3_d.jpg', dpi = 500, height = 17, width = 20)
## stalk marker plot
general.marker.df <- marker.df[which(marker.df$Gene %in% marker.list$general),]
ggplot(general.marker.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Expression)) +  
  scale_color_gradient2(low = 'green', mid = 'grey', midpoint = 0.5, high = 'purple') +
  facet_grid(vars(Experiment), vars(Gene)) +
  ggtitle('Marker Expression') +
  plot.theme
ggsave(filename = 'paper_components/fig-s3_e.jpg', dpi = 500, height = 17, width = 15)
## seurat clusters
ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(color = Seurat), size = 2) +
  facet_wrap(vars(Experiment), ncol = 3) +
  ggtitle('Seurat Clusters') +
  plot.theme + 
  theme(legend.key.size = unit(1.5, 'cm')) +
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave(filename = 'paper_components/fig-s3_f.jpg', dpi = 500, height = 10, width = 32)
## 3D plot
anchor.gexp <- RunUMAP(anchor.gexp, reduction = "pca", dims = 1:30, n.components = 3)
anchor.3d.umap <- anchor.gexp@reductions$umap@cell.embeddings
jpeg('paper_components/fig-s3_g.jpg', height = 750, width = 750)
scatter3D(anchor.3d.umap[,1], anchor.3d.umap[,2], anchor.3d.umap[,3],
          colvar = as.integer(plot.dat$Seurat),
          col = ClusterColors(length(levels(plot.dat$Seurat))),
          xlab = 'UMAP1', ylab = 'UMAP2', zlab = 'UMAP3', main = "Seurat Clusters")
dev.off()
###############

########## OLD FIGURE VERSIONS ##########

### fig 2 OLD
###############
exp.label <- factor(sort(exp.vec), levels = c('WT', 'KO', 'KOF4'))
marker.df <- data.frame('UMAP1' = anch.umap[,1], 'UMAP2' = anch.umap[,2],
                        'Experiment' = exp.label)
# add protein activity
marker.set <- intersect(unlist(marker.list), rownames(anch.vip))
marker.activity <- unlist(as.list(t(anch.vip[marker.set,])))
marker.label <- factor(rep(marker.set, each = nrow(anch.umap)), levels = rev(unique(marker.set)))
marker.df <- do.call('rbind', replicate(length(marker.set), marker.df, simplify = FALSE))
marker.df$Activity <- marker.activity
marker.df$Protein <- marker.label
## WT plot
wt.cells <- which(marker.df$Experiment == 'WT')
wt.markers <- ggplot(marker.df[wt.cells,], aes(UMAP1, UMAP2)) + geom_point(aes(color = Activity)) +
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  facet_wrap(~Protein, nrow = 3, ncol = 3) +
  ggtitle('WT Marker Activity') +
  plot.theme + xlim(x.umap.lim) + ylim(y.umap.lim)
print(wt.markers)
ggsave(filename = 'paper_components/fig-2_a.jpg', dpi = 500, height = 17, width = 17)
## KO plot
ko.cells <- which(marker.df$Experiment == 'KO')
ko.markers <- ggplot(marker.df[ko.cells,], aes(UMAP1, UMAP2)) + geom_point(aes(color = Activity)) +
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  facet_wrap(~Protein, nrow = 3, ncol = 3) +
  ggtitle('KO Marker Activity') +
  plot.theme + xlim(x.umap.lim) + ylim(y.umap.lim)
print(ko.markers)
ggsave(filename = 'paper_components/fig-2_b.jpg', dpi = 500, height = 17, width = 17)
## KOF4 plot
kof4.cells <- which(marker.df$Experiment == 'KOF4')
kof4.markers <- ggplot(marker.df[kof4.cells,], aes(UMAP1, UMAP2)) + geom_point(aes(color = Activity)) +
  scale_color_gradient2(low = 'blue', mid = 'grey', midpoint = 0.5, high = 'red') +
  facet_wrap(~Protein, nrow = 3, ncol = 3) +
  ggtitle('KOF4 Marker Activity') +
  plot.theme + xlim(x.umap.lim) + ylim(y.umap.lim)
print(kof4.markers)
ggsave(filename = 'paper_components/fig-2_c.jpg', dpi = 500, height = 17, width = 17)
###############
