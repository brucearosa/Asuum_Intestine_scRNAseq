# Load necessary libraries
library(Seurat)
library(future)
library(reshape2)
library(clustree)
library(scales)
library(ggplot2)

# Function to integrate Seurat objects from different samples
# Allows specifying a reference sample
integrate_samples <- function(samples, ref = 1) {
  for(i in samples) {
    # Read previously prepared and processed sample objects from RDS files
    assign(i, readRDS(paste(i, ".rds", sep = "")))
    # Find most variable and informative genes for each sample
    assign(i, FindVariableFeatures(get(i)))
  }
  # Select integration features and prepare for SCT integration
  int.features <- SelectIntegrationFeatures(object.list = mget(samples), nFeatures = 2000)
  int.list <- PrepSCTIntegration(object.list = mget(samples), anchor.features = int.features)
  
  # Increase memory limit and request multiple processors for FindIntegrationAnchors
  options(future.globals.maxSize = 4800000000)
  plan("multiprocess", workers = 8)
  
  # Find integration anchors
  int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", anchor.features = int.features, reference = ref)
  
  # Finding the overlapping gene set for integration
  subgenes <- Reduce(intersect, lapply(mget(samples), row.names))
  
  # Increase memory limit for IntegrateData
  options(future.globals.maxSize = 6400000000)
  
  # Integrate data and perform dimensionality reduction
  all <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", features.to.integrate = subgenes)
  all <- RunPCA(all) # Run PCA
  all <- RunTSNE(all, dims = 1:30) # Run t-SNE using 30 PCs
  all <- RunUMAP(all, dims = 1:30) # Run UMAP using 30 PCs
  return(all)
}

# Function to find clusters in an integrated Seurat object
# Uses resolution range from 0.1 to 1, removes results with >15 clusters
find_clusters <- function(all) {
  all <- FindNeighbors(all, dims = 1:30)
  options(future.globals.maxSize = 2400000000)
  plan("multiprocess", workers = 8)
  all <- FindClusters(all, resolution = (1:10)/10, save.SNN = TRUE)
  all@meta.data[,names(which(apply(all@meta.data[,grep("integrated_snn_res", colnames(all@meta.data))], 2, function(x) {length(levels(as.factor(x))) > 15})))] <- NULL
  for(i in colnames(all@meta.data)[grep("integrated_snn_res", colnames(all@meta.data))]){
    levels(all@meta.data[,i]) = as.character(1:length(levels(all@meta.data[,i])))
  }
  return(all)
}

makeplots <- function(all, name = "all", reps = 2, clusts = 9) {
  # UMAP plots colored by treatment
  pdf(paste(name, "umap.pdf", sep = "."))
  DimPlot(all, reduction = "umap", group.by = "treat") %>% print()
  dev.off()
  
  # UMAP colored by sample, split by treatment
  pdf(paste(name, "umap.splitbytreat.pdf", sep = "."), width = length(unique(all@meta.data$orig.ident)) * 2.5, height = 5)
  DimPlot(all, reduction = "umap", group.by = "orig.ident", split.by = "treat") %>% print()
  dev.off()
  
  # Clustree plot
  pdf(paste(name, "clustree.pdf", sep = "."))
  clustree(all, prefix = "integrated_snn_res.") %>% print()
  dev.off()
  
  # UMAP plot colored by cluster identity with labels
  pdf(paste(name, ".umap.", clusts, "clust.pdf", sep = ""))
  DimPlot(all, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.1, label.size = 5) + NoLegend() %>% print()
  dev.off()
  
  # Preparing data for bar plots
  dat <- melt(table(all@meta.data$seurat_clusters, all@meta.data$treat))
  colnames(dat) <- c("Cluster", "Treatment", "Cell_Count")
  cols <- hue_pal()(length(unique(dat$Treatment)))

  # Barplot for cell count for each cluster
  pdf(paste(name, "treatment.barplot.pdf", sep = "."))
  ggplot(dat, aes(x = Cluster, y = Cell_Count, fill = Treatment)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    xlab("\nCluster") +
    ylab("Cell Count\n") %>% print()
  dev.off()

  # Barplot for cell fractions for each cluster, colored by treatment
  pdf(paste(name, "treatment.fractionbarplot.pdf", sep = "."))
  ggplot(dat, aes(x = Cluster, y = Cell_Count, fill = Treatment)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    xlab("\nCluster") +
    ylab("Fraction cells\n") %>% print()
  dev.off()

  # Generate and format data for normalized cell counts and fractions
  tbl <- table(all@meta.data$seurat_clusters, all@meta.data$treat) 
  normtbl <- melt(1e4 * prop.table(tbl, 2))
  colnames(normtbl) <- c("Cluster", "Treatment", "Cell_prop")

  # Barplot for normalized cell counts for each cluster, colored by treatment
  pdf(paste(name, "treatment.normalizedbarplot.pdf", sep = "."))
  ggplot(normtbl, aes(x = Cluster, y = Cell_prop, fill = Treatment)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    xlab("\nCluster") +
    ylab("Normalized count\n") %>% print()
  dev.off()

  # Barplot for normalized cell fractions for each cluster, colored by treatment
  pdf(paste(name, "treatment.normalizedfractionbarplot.pdf", sep = "."))
  ggplot(normtbl, aes(x = Cluster, y = Cell_prop, fill = Treatment)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    xlab("\nCluster") +
    ylab("Fraction cells\n") %>% print()
  dev.off()

  # If replicates are present, additional plots can be generated to compare them.
  # This part of the function could be expanded with specific logic for handling replicates.
}

# Function for differential expression analysis using MAST
DEfiles <- function(all, name, treat, ctrl, super1, super2) {
  options(future.globals.maxSize = 2400000000)
  plan("multiprocess", workers = 8)
  
  # Treatment vs. all others comparison using MAST
  trvall <- FindMarkers(all, ident.1 = treat, test.use = "MAST", group.by = "seurat_clusters")
  
  # Treatment vs. control comparison using MAST
  trvctrl <- FindMarkers(all, ident.1 = treat, ident.2 = ctrl, test.use = "MAST", group.by = "seurat_clusters")
  
  # Supercluster comparison using MAST
  supers <- FindMarkers(all, ident.1 = super1, ident.2 = super2, test.use = "MAST", group.by = "seurat_clusters")
  
  # Write DE analysis results to files
  write.table(trvall, file = paste(name, ".MAST.", treat, "vall.tsv", sep = ""), quote = FALSE, sep = "\t")
  write.table(trvctrl, file = paste(name, ".MAST.", treat, "v", ctrl, ".tsv", sep = ""), quote = FALSE, sep = "\t")
  write.table(supers, file = paste(name, ".MAST.1v2superclusters.tsv", sep = ""), quote = FALSE, sep = "\t")
}

# Function to generate a table of normalized fractions for analysis
fractable <- function(all) {
  tbl <- table(all@meta.data$seurat_clusters, all@meta.data$orig.ident)
  normtbl <- prop.table(tbl, 1) * 10000
  normtbl <- as.data.frame(normtbl)
  colnames(normtbl) <- levels(all@meta.data$orig.ident)
  rownames(normtbl) <- levels(all@meta.data$seurat_clusters)
  
  # Convert to fractions
  fractable <- sweep(normtbl, 1, rowSums(normtbl), FUN = "/")
  return(fractable)
}

# Function to aggregate p-values from multiple analyses
aggregatePvals <- function(pvalMatrix, method = "fishers", pAdjustMethod = "BH", order = TRUE) {
  if(!is.matrix(pvalMatrix)) stop("pvalMatrix must be a matrix.")
  if(is.null(rownames(pvalMatrix))) stop("pvalMatrix must have row names.")
  if(!(method %in% c("stouffers", "fishers"))) stop('method must be "stouffers" or "fishers".')
  
  aggr.pval <- numeric(nrow(pvalMatrix))
  names(aggr.pval) <- rownames(pvalMatrix)
  
  if (method == "fishers") {
    pvalMatrixLogged <- -2 * log(pvalMatrix)
    aggr.pval <- apply(pvalMatrixLogged, 1, sum)
    aggr.pval <- pchisq(aggr.pval, df = 2 * ncol(pvalMatrix), lower.tail = FALSE)
  } else if (method == "stouffers") {
    z.scores <- apply(pvalMatrix, 1, function(x) sum(qnorm(x)) / sqrt(length(x)))
    aggr.pval <- pnorm(z.scores, lower.tail = FALSE)
  }
  
  adj.aggr.pval <- p.adjust(aggr.pval, method = pAdjustMethod)
  aggr.pval.table <- cbind(aggr.pval, adj.aggr.pval)
  rownames(aggr.pval.table) <- names(aggr.pval)
  colnames(aggr.pval.table) <- c("Aggregated p-value", "Adjusted aggregated p-value")
  
  if(order) aggr.pval.table <- aggr.pval.table[order(aggr.pval.table[, "Aggregated p-value"]), ]
  return(aggr.pval.table)
}

# Function to find DE genes, calculate average expression, and write to files
DEexp <- function(obj, name1 = 1, name2 = 2, nameobj) {
  # Check for existing files to avoid overwriting
  fileDE <- paste(nameobj, ".MAST.", name1, "v", name2, ".tsv", sep = "")
  if(file.exists(fileDE)) stop(paste("File", fileDE, "already exists! Exiting."))
  
  fileExp <- paste(nameobj, name1, name2, "avgexpression.tsv", sep = ".")
  if(file.exists(fileExp)) stop(paste("File", fileExp, "already exists! Exiting."))
  
  # Differential expression analysis
  plan("multiprocess", workers = 8)
  DE12 <- FindMarkers(obj, ident.1 = name1, ident.2 = name2, test.use = "MAST", group.by = "DE.Idents")
  write.table(DE12, file = fileDE, quote = FALSE, sep = "\t")
  
  # Average expression calculation
  avgExp <- AverageExpression(obj, return.seurat = TRUE)
  write.table(avgExp$RNA, file = fileExp, quote = FALSE, sep = "\t", row.names = TRUE)
}
