# Reading GEM table and setting rownames + colnames
gem <- read.table('gem.tbl', sep='\t', header = FALSE)
rownames(gem) <- gem[, 1]
gem <- gem[, -1]
colnames_gem <- read.csv('colnames.txt', sep='\n', header = FALSE)
colnames(gem) <- t(colnames_gem)

# separating mice
# gem_t <- data.frame(gem[1:64], gem[65:128], gem[257:320])
# gem_n <- data.frame(gem[129:192], gem[193:256], gem[321:384])

# Creating sce
sce <- SingleCellExperiment(assays = list(counts = as.matrix(gem, mode="numeric")))
rowData(sce)$feature_symbol <- rownames(sce)

# QC steps
sce <- calculateQCMetrics(sce, nmads = 3, pct_feature_controls_threshold = 80)
# Plots
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="",
     breaks=30, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="",
     breaks=30, col="grey80", ylab="Number of cells")
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotQC(sce, type = "highest-expression", n=30) + fontsize

# Removing outlier cells with QC
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
sce <- sce[,!(libsize.drop | feature.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
Remaining=ncol(sce))

# Normalization steps
sce <- computeSumFactors(sce)
summary(sizeFactors(sce)) # Summary info
sce <- normalize(sce)

# SC3
sce <- sc3(sce, ks = 2:4, pct_dropout_min = 10, pct_dropout_max = 90, d_region_min = 0.04, d_region_max = 0.07, svm_max = 5000, kmeans_nstart = 50, kmeans_iter_max = 50, biology = TRUE)
sc3_interactive(sce)
# Optional
sc3_plot_consensus(sce, k = 3)
sc3_plot_cluster_stability(sce, k = 3)
sc3_plot_silhouette(sce, k = 3)

# N stats GEM
#   ByLibSize ByFeature Remaining
# 1        14         4       177
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.04532 0.45323 0.82502 1.00000 1.34833 5.12583 

# T stats GEM
# ByLibSize ByFeature Remaining
# 1         9         0       183
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09064 0.45549 0.83086 1.00000 1.34214 3.26550

# N stats GTM
# ByLibSize ByFeature Remaining
# 1        10         4       181
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09114 0.45584 0.83909 1.00000 1.37377 4.54675 

# T stats GTM
# ByLibSize ByFeature Remaining
# 1         9         1       182
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1255  0.4887  0.8815  1.0000  1.3397  2.7860 

# scatter plot
logcounts <- assays(sce)$logcounts
lyc_norm <- logcounts['Ly6c2_e',]
nra_norm <- logcounts['Nr4a1_e',]
cols <- colnames(sce)
clusters <- colData(sce)$sc3_3_clusters
plt <- data.frame(cols, lyc_norm, nra_norm, clusters)
rownames(plt) <- t(plt["cols"])
plt <- plt[, -1]
ggscatter(plt, x = "lyc_norm", y = "nra_norm", color = "clusters",
          palette = c("#00AFBB", "#E7B800", "#FC4E07") )

# comparing clusters GEM
# [1] 191

# comparing clusters GTM
# [1] 193


real_clusters_sce_N[1:64] = 1
real_clusters_sce_N[65:128] = 1
real_clusters_sce_N[257:320] = 1
real_clusters_sce_N[129:192] = 2
real_clusters_sce_N[193:256] = 2
real_clusters_sce_N[321:384] = 2
real_clusters_sce_T[1:64] = 2
real_clusters_sce_T[65:128] = 2
real_clusters_sce_T[257:320] = 2
real_clusters_sce_T[129:192] = 1
real_clusters_sce_T[193:256] = 1
real_clusters_sce_T[321:384] = 1
