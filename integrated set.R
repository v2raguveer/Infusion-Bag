library(Seurat)
library(cowplot)
ctrl.data <- Read10X(data.dir = "/Users/vanitharaguveer/Documents/RData/cellrangerFiltered/s284_filterMatrix")
acti.data <- Read10X(data.dir = "/Users/vanitharaguveer/Documents/RData/cellrangerFiltered/s285_filterMatrix")

ctrl <- CreateSeuratObject(counts = ctrl.data, project = "Immune_ctrl", min.cells = 5)
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA>500 & nFeature_RNA<3500)
ctrl <- PercentageFeatureSet(ctrl, pattern = "^MT-", col.name = "percent.mt")
ctrl <- SCTransform(ctrl, vars.to.regress = "percent.mt", verbose = FALSE)

stim <- CreateSeuratObject(counts = acti.data, project = "Immune_stim", min.cells = 5)
stim$stim <- "STIM"
stim <- subset(stim, subset = nFeature_RNA>500 & nFeature_RNA<3500)
stim <- PercentageFeatureSet(stim, pattern = "^MT-", col.name = "percent.mt")
stim <- SCTransform(stim, vars.to.regress = "percent.mt", verbose = FALSE)

immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl,stim), dims=1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose=FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = 'pca', dims = 1:20)
immune.combined <- RunTSNE(immune.combined, reduction = 'pca', dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = 'pca', dims=1:20)
immune.combined <- FindClusters(immune.combined, resolution = .5)
p1 <- DimPlot(immune.combined, reduction = 'umap', group.by = 'stim')
p2 <- DimPlot(immune.combined, reduction = 'umap', label = TRUE)
p3 <- DimPlot(immune.combined, reduction = 'tsne', group.by = 'stim')
p4 <- DimPlot(immune.combined, reduction = 'tsne', label = TRUE)
plot(p1)
plot(p2)
plot(p3)
plot(p4)
DimPlot(immune.combined, reduction = 'tsne', split.by = "stim")

DefaultAssay(immune.combined)