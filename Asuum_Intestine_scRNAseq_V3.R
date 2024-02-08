# Load necessary libraries
library(Seurat)
library(future)
library(reshape2)
library(clustree)
library(scales)
library(ggplot2)



# ____ PRIMARY ANALYSIS WORKFLOW ____
# 1. read the sample filtered data separately and make separate Seurat objects using standard filters (min.cells,min.features. I used min.features 100 rather than 200 to avoid filtering too many cells from "lower mode"). normalize using SCTransform
# 2. find anchors and integrate them into s134 intergated Seurat object
# 3. run PCA, identify clusters using multiple resolutions. Used resolution 0.4 and find doublets in s134
# 4. Make new S1,S3,S4 objects (as in step 1 above). filter the doublets found in step 3 above from the corresponding sample object.
# 5. calculate cell cycle scores and regress them out
# 6. Integrate the samples to cluster and plot them together
# 7. run PCA,tSNE,UMAP. identify clusters using multiple resolutions. 
# 8. plot tSNE, violinplot. plot clustreee. used resolution 0.1 and 0.4 finally.
# 9. transfer meta data from s134 to individual sample objects. Then merge objects
# 10. find differentially expressed genes. using MAST. Also output cluster averages.

# s1, s3 and s4 refer to the F1_UT, F2_UT and F3_UT samples, respectively. 

# Step 1. read the sample filtered data separately and make separate Seurat objects using standard filters (min.cells,min.features. I used min.features 100 rather than 200 to avoid filtering too many cells from "lower mode"). normalize using SCTransform
s1d<-Read10X(data.dir="<path_to_S1_filtered_feature_bc_matrix>")
s3d<-Read10X(data.dir="<path_to_S3_filtered_feature_bc_matrix>")
s4d<-Read10X(data.dir="<path_to_S4_filtered_feature_bc_matrix>")
s1<-CreateSeuratObject(counts = s1d,project="S1",min.cells = 3, min.features = 100)
s3<-CreateSeuratObject(counts = s3d,project="S3",min.cells = 3, min.features = 100)
s4<-CreateSeuratObject(counts = s4d,project="S4",min.cells = 3, min.features = 100)
# Remove cells with high occurence of mitochondrial genes
mtgenes=c("AgB02-g498","AgE51-g006","AgE51-g007","AgR009X-g325","AgR009X-g327","AgR009X-g326","AgR022-g083")
mtgenes=paste("gene:",mtgenes,sep="")
s1[["percent.mt"]] <- PercentageFeatureSet(s1, features=mtgenes)
s1<-subset(s1,subset = percent.mt < 1)
s3[["percent.mt"]] <- PercentageFeatureSet(s3, features=mtgenes)
s3<-subset(s3,subset = percent.mt < 1)
s4[["percent.mt"]] <- PercentageFeatureSet(s4, features=mtgenes)
s4<-subset(s4,subset = percent.mt < 1)
# Normalize using SCTransform
for(i in c("s1","s3","s4")){assign(i,SCTransform(get(i),verbose=FALSE))}

# Step 2. find anchors and integrate them into s134 intergated Seurat object
int.features<-SelectIntegrationFeatures(object.list = c(s1,s3,s4), nFeatures=2000)
int.list <- PrepSCTIntegration(object.list = c(s1,s3,s4),anchor.features = int.features)
int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", anchor.features=int.features)
s134<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT")

# Step 3. run PCA, identify clusters using multiple resolutions. Used resolution 0.4 and find doublets in s134
s134<-RunPCA(s134,verbose=FALSE)
ElbowPlot(s134) # Point of maximum curvature gives rough estimate of no. clusters
# Prepare data to find Doublets
sweep.res.list<-paramSweep_v3(s134,PCs=1:20,sct=TRUE,num.cores=8)
sweep.stats<-summarizeSweep(sweep.res.list,GT=FALSE)
bcmvn<-find.pK(sweep.stats)
bcmvn
s134<-FindNeighbors(s134,dims=1:30)
s134<-FindClusters(s134,resolution= c(0.4, 0.8, 1.2),save.SNN=TRUE)
s134<-RunTSNE(s134,dims=1:30)
s134<-RunUMAP(s134,dims=1:30)
# export tSNE plot and UMAP plot. the res.0.4 is selected after looking at 
# the resolutions described above. The idea is to use the resolution that 
# gives a reasonable number of clusters that can be interpreted. I have 
# usually found ~8-10 for a high resolution look and ~3-4 for a low resolution 
# look. although in this case it was 17 clusters.
DimPlot(s134,reduction="tsne",group.by="integrated_snn_res.0.4")
dev.copy(pdf,"S134.withdoublets.tsne.17clustres0.4.pdf")
dev.off()
DimPlot(s134,reduction="umap",group.by="integrated_snn_res.0.4")
dev.copy(pdf,"S134.withdoublets.umap.17clustres0.4.pdf")
dev.off()
# Find Doublets
annotations<-s134@meta.data$integrated_snn_res.0.4
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*length(s134@meta.data$integrated_snn_res.0.4))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# using 7.5% default
s134<-doubletFinder_v3(s134,PCs=1:30, pN=0.25, pK=0.001, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
head(s134@meta.data,n=2) # finding the number (1895 in this case)
# using adjustment for homotypic proportion
s134<-doubletFinder_v3(s134,PCs=1:30, pN=0.25, pK=0.001, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.001_1895", sct = TRUE)
DimPlot(s134,reduction="tsne",group.by="DF.classifications_0.25_0.001_1695")
dev.copy(pdf,"S134.withdoublets.tsne.doublets1695.pdf")
dev.off()
VlnPlot(s134,features="nFeature_RNA",pt.size=0,group.by="DF.classifications_0.25_0.001_1695")
dev.copy(pdf,"S134.withdoublets.genes.doublets1695.violin.nopoints.pdf")
dev.off()

# Step 4. filter the doublets found in step 3 above from the corresponding sample object.
s1.sing<-gsub("_1","",row.names(s134@meta.data[s134@meta.data$orig.ident == "S1" & s134@meta.data$DF.classifications_0.25_0.001_1695 == "Singlet",]))
s3.sing<-gsub("_2","",row.names(s134@meta.data[s134@meta.data$orig.ident == "S3" & s134@meta.data$DF.classifications_0.25_0.001_1695 == "Singlet",]))
s4.sing<-gsub("_3","",row.names(s134@meta.data[s134@meta.data$orig.ident == "S4" & s134@meta.data$DF.classifications_0.25_0.001_1695 == "Singlet",]))
s1 <- subset(s1,cells=s1.sing)
s3 <- subset(s3,cells=s3.sing)
s4 <- subset(s4,cells=s4.sing)

# Step 5. calculate CC score and regress it out. For this, the gene names are read from the files As.s.genes.manual (for S phase genes) and As.g2m.genes.manual for G2M genes.
s1@active.assay="RNA"
s3@active.assay="RNA"
s4@active.assay="RNA"
s1[["SCT"]] <-NULL
s3[["SCT"]] <-NULL
s4[["SCT"]] <-NULL
for(i in c("s1","s3","s4")){assign(i,SCTransform(get(i),variable.features.n=13000,assay="RNA",new.assay.name="SCT",verbose=FALSE))}
Ass<-read.table("As.s.genes.manual",header=F,sep="\t")
Ass<-as.character(Ass$V1[grep("[0-9]",Ass$V1)])
Ass<-gsub("^","gene:",gsub("_","-",Ass))
Asg2m<-read.table("As.g2m.genes.manual",header=F,sep="\t")
Asg2m<-as.character(Asg2m$V1[grep("[0-9]",Asg2m$V1)])
Asg2m<-gsub("^","gene:",gsub("_","-",Asg2m))
Asg2m<-gsub("-t[0-9]+$","",Asg2m)
s1<-CellCycleScoring(s1,s.features = Ass,g2m.features=Asg2m)
s3<-CellCycleScoring(s3,s.features = Ass,g2m.features=Asg2m)
s4<-CellCycleScoring(s4,s.features = Ass,g2m.features=Asg2m)
DimPlot(s1,group.by="Phase")
s1$CC.Difference<-s1$S.Score-s1$G2M.Score
s3$CC.Difference<-s3$S.Score-s3$G2M.Score
s4$CC.Difference<-s4$S.Score-s4$G2M.Score
s1<-ScaleData(s1,features=rownames(s1),vars.to.regress="CC.Difference")
s3<-ScaleData(s3,features=rownames(s3),vars.to.regress="CC.Difference")
s4<-ScaleData(s4,features=rownames(s4),vars.to.regress="CC.Difference")

# Step 6. Integrate the samples to cluster and plot them together
int.features<-SelectIntegrationFeatures(object.list = c(s1,s3,s4), nFeatures=2000)
int.list <- PrepSCTIntegration(object.list = c(s1,s3,s4),anchor.features = int.features)
int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", anchor.features=int.features)
library(dplyr)
all_features<-lapply(int.list,row.names) %>% Reduce(intersect, .)
rm(s134)
s134<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT",features.to.integrate = all_features)

# Step 7. run PCA,tSNE,UMAP. identify clusters using multiple resolutions.
s134<-RunPCA(s134,verbose=FALSE)
ElbowPlot(s134)
s134<-RunTSNE(s134,dims=1:30)
s134<-RunUMAP(s134,dims=1:30)
s134<-FindNeighbors(s134,dims=1:30)
s134<-FindClusters(s134,resolution= c(0.1,0.4, 0.8),save.SNN=TRUE)
levels(s134@meta.data$integrated_snn_res.0.1) <-as.character(1:3)
levels(s134@meta.data$integrated_snn_res.0.4) <-as.character(1:7)
levels(s134@meta.data$integrated_snn_res.0.8) <-as.character(1:11)
Idents(s134)=s134@meta.data$integrated_snn_res.0.4

# Step 8. plot tSNE, violinplot. plot clustreee. used resolution 0.1 and 0.4 finally.
DimPlot(s134,reduction="tsne",group.by="integrated_snn_res.0.4")
dev.copy(pdf,"S134.tsne.7clustres0.4.pdf")
dev.off()
DimPlot(s134,reduction="tsne",split.by="orig.ident",group.by="integrated_snn_res.0.4")
dev.copy(pdf,"S134.tsne.7clustres0.4.splitbysample.pdf")
dev.off()
library(clustree)
clustree(s134,prefix="integrated_snn_res.",node_colour = "sc3_stability")
dev.copy(pdf,"S134.clustree.pdf")
dev.off()
DimPlot(s134,reduction="tsne",group.by="integrated_snn_res.0.1")
dev.copy(pdf,"S134.tsne.3clustres0.1.pdf")
dev.off()
VlnPlot(s134,features="nFeature_RNA",pt.size=0,group.by="integrated_snn_res.0.4")
dev.copy(pdf,"S134.genes.7clustres0.4.violin.nopoints.pdf")
dev.off()
DimPlot(s134,reduction="umap",group.by="integrated_snn_res.0.4")
dev.copy(pdf,"S134.umap.7clustres0.4.pdf")
dev.off()
DimPlot(s134,reduction="umap",split.by="orig.ident",group.by="integrated_snn_res.0.4")
dev.copy(pdf,"S134.umap.7clustres0.4.splitbysample.pdf")
dev.off()
DimPlot(s134,reduction="umap",group.by="integrated_snn_res.0.1")
dev.copy(pdf,"S134.umap.3clustres0.1.pdf")
dev.off()
# plot violin plot with sample content coloring
S1clusts<-summary(s134@meta.data[s134@meta.data$orig.ident=="S1",]$integrated_snn_res.0.4)
S3clusts<-summary(s134@meta.data[s134@meta.data$orig.ident=="S3",]$integrated_snn_res.0.4)
S4clusts<-summary(s134@meta.data[s134@meta.data$orig.ident=="S4",]$integrated_snn_res.0.4)
r=S1clusts/(S1clusts+S3clusts+S4clusts)
g=S3clusts/(S1clusts+S3clusts+S4clusts)
b=S4clusts/(S1clusts+S3clusts+S4clusts)
p<-VlnPlot(s134,features="nFeature_RNA",pt.size=0,group.by="integrated_snn_res.0.4",col=rgb(r,g,b),y.max=6000)+
  annotate(geom="rect",xmin=(1:7)-0.2,xmax=(1:7)+0.2,ymin=4500-((S1clusts+S3clusts+S4clusts)*500/7700),ymax=4500+((S1clusts+S3clusts+S4clusts)*500/7700),alpha=0.2)+
  annotate(geom="text",x=1:7,y=4500,label=(S1clusts+S3clusts+S4clusts),angle=90,size=3)
p<-p+annotate(geom="text",x=1-0.2,y=5500,label=paste("S1=",round(r[1],2)),angle=90,size=2.5,col="red")+annotate(geom="text",x=(2:7)-0.2,y=5500,label=round(r[2:7],2),angle=90,size=2.5,col="red")+annotate(geom="text",x=1,y=5500,label=paste("S3=",round(g[1],2)),angle=90,size=2.5,col="green")+annotate(geom="text",x=2:7,y=5500,label=round(g[2:7],2),angle=90,size=2.5,col="green")+annotate(geom="text",x=1+0.2,y=5500,label=paste("S4=",round(b[1],2)),angle=90,size=2.5,col="blue")+annotate(geom="text",x=(2:7)+0.2,y=5500,label=round(b[2:7],2),angle=90,size=2.5,col="blue")
p
dev.copy(pdf,"S134.genes.7clustres0.4.colorbysample.violin.nopoints.pdf")
dev.off()

# Step 9. transfer meta data from s134 to individual sample objects. Then merge objects.
s1meta<-s134@meta.data[s134@meta.data$orig.ident=="S1",]
s3meta<-s134@meta.data[s134@meta.data$orig.ident=="S3",]
s4meta<-s134@meta.data[s134@meta.data$orig.ident=="S4",]
row.names(s1meta)<-gsub("_1$","",row.names(s1meta))
row.names(s3meta)<-gsub("_2$","",row.names(s3meta))
row.names(s4meta)<-gsub("_3$","",row.names(s4meta))
s1@meta.data$seurat_clusters=s1meta$integrated_snn_res.0.1
s3@meta.data$seurat_clusters=s3meta$integrated_snn_res.0.1
s4@meta.data$seurat_clusters=s4meta$integrated_snn_res.0.1
s134<-merge(s1,c(s3,s4),add.cell.ids=c("S1","S3","S4"))

# Step 10. find differentially expressed genes. I tried wilcox, MAST and DESeq2 and also found their intersections. But finally decided to only use MAST. Also output cluster averages.
for(i in 1:3){
  assign(paste("MAST",".",i,"vall",sep=""),FindMarkers(s134,ident.1=i,test.use="MAST",group.by="seurat_clusters",min.pct=0,logfc.threshold=0))};for(i in 2:3){assign(paste("MAST",".1v",i,sep=""),FindMarkers(s134,ident.1=1,ident.2=i,test.use="MAST",group.by="seurat_clusters",min.pct=0,logfc.threshold=0))};MAST.2v3<-FindMarkers(s134,ident.1=2,ident.2=3,test.use="MAST",group.by="seurat_clusters",min.pct=0,logfc.threshold=0)
for(i in ls()[grep("MAST.*",ls())]){write.table(get(i),file=paste(i,".tsv",sep=""),quote=FALSE,row.names=TRUE,sep="\t")}
Idents(s134)<-s134@meta.data$seurat_clusters
cluster.averages<-AverageExpression(s134)
write.table(cluster.averages$SCT,file="S134.3clusteraverages.tsv",sep="\t",quote=FALSE,row.names=TRUE)




# ____ USEFUL FUNCTIONS ____
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
