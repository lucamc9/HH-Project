# Reading GEM table and setting rownames + colnames
gtm <- read.table('GTM_filtered.tbl', sep='\t', header = FALSE)
rownames(gtm) <- gtm[, 1]
gtm <- gtm[, -1]
colnames_gem <- read.csv('colnames.txt', sep='\n', header = FALSE)
colnames(gtm) <- t(colnames_gem)

# Creating sce
sct <- SingleCellExperiment(assays = list(counts = as.matrix(gtm, mode="numeric")))
rowData(sct)$feature_symbol <- rownames(sct)

# QC steps
sct <- calculateQCMetrics(sct, nmads = 3, pct_feature_controls_threshold = 80)
# Plots
# par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sct$total_counts/1e6, xlab="Library sizes (millions)", main="",
     breaks=30, col="grey80", ylab="Number of cells")
hist(sct$total_features, xlab="Number of expressed genes", main="",
     breaks=30, col="grey80", ylab="Number of cells")
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotQC(sct, type = "highest-expression", n=50) + fontsize

# Removing outlier cells with QC
libsize.drop <- isOutlier(sct$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sct$total_features, nmads=3, type="lower", log=TRUE)
sct <- sct[,!(libsize.drop | feature.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           Remaining=ncol(sct))

# Normalization steps
sct <- computeSumFactors(sct)
summary(sizeFactors(sct)) # Summary info
sct <- normalize(sct)

# SC3
sct <- sc3(sct, ks = 2:4, gene_filter = FALSE, pct_dropout_min = 10, pct_dropout_max = 90, d_region_min = 0.04, d_region_max = 0.07, svm_max = 5000, kmeans_nstart = 50, kmeans_iter_max = 50, biology = TRUE)
# sc3_interactive(sct)
# Optional
sc3_plot_consensus(sct, k = 3)
# sc3_plot_cluster_stability(sct, k = 3)
# sc3_plot_silhouette(sct, k = 3)
# getting counts > counts <- assays(sct)$counts > counts['lyc62_e', ]
