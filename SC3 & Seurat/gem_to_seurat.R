library(Seurat)

# Reading GEM table and setting rownames + colnames
gem <- read.table('gem.tbl', sep='\t', header = FALSE)
rownames(gem) <- gem[, 1]
gem <- gem[, -1]
colnames_gem <- read.csv('colnames.txt', sep='\n', header = FALSE)
colnames(gem) <- t(colnames_gem)

# init Seurat object
seu <- CreateSeuratObject(raw.data = gem, min.cells = 3, min.genes = 200)

# plot nGene and nUMI 
VlnPlot(object = seu, features.plot = c("nGene", "nUMI"), nCol = 2)

# gene-gene relationship plot
GenePlot(object = seu, gene1 = "nUMI", gene2 = "nGene")

# gene filter / NowmalizeData only works if this step is skipped
seu <- FilterCells(object = seu, subset.names = c("nGene"), 
                   low.thresholds = c(200), high.thresholds = c(2500))

# normalizing the data
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

# detecting variable genes
seu <- FindVariableGenes(object = seu, mean.function = ExpMean, 
                         dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, 
                         y.cutoff = 0.5)
length(x = seu@var.genes) # [1] 8042

# scaling the data
seu <- ScaleData(object = seu, vars.to.regress = c("nUMI"))

# performing pca
seu <- RunPCA(object = seu, pc.genes = seu@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PrintPCA(object = seu, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = seu, pcs.use = 1:2)
PCAPlot(object = seu, dim.1 = 1, dim.2 = 2)
seu <- ProjectPCA(object = seu, do.print = FALSE)
PCHeatmap(object = seu, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = seu, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

# jackstraw and elbow plot
seu <- JackStraw(object = seu, num.replicate = 100)
JackStrawPlot(object = seu, PCs = 1:12)
PCElbowPlot(object = seu)

# clustering cell types
seu <- FindClusters(object = seu, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = seu)

# t-SNE
seu <- RunTSNE(object = seu, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = seu)

# GEM
# sum(seu@ident == realIdent)
# [1] 194

# GTM
# sum(seu@ident == realIdent)
# [1] 201