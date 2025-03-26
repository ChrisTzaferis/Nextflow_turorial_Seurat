library(dplyr)
library(Seurat)
library(patchwork)
library(repr)
library(ggplot2)

# Description:
#   This function takes as an input a folder generated from 10x cellranger.
#   It performs QC filtering, Normalization, detection of HVGs, PCA,
#   clustering and UMAP visualization utilizing Seurat package workflow

# Parameters:
#   inputDir: a path pointing to a cellranger directory containing
#features.tsv (or genes.tsv), matrix.mtx, barcodes.tsv

#   outputName: the processed filed is saved in the current working
#directory as "outputName.RDS"

#   mtPattern: a prefix for mitochondrial genes use "^MT-" for
#human data or "^mt-" for mouse data

args <- commandArgs(trailingOnly = TRUE)

# Read arguments from Nextflow
input_Dir <- args[1]
output_Name <- args[2]
mt_Pattern <- args[3]

seuratAnalysis <- function(inputDir, outputName, mtPattern)
{
  print(getwd())
  #set seed
  set.seed(42)

  # Load dataset
  obj.data <- Read10X(data.dir = inputDir)

  # Initialize the Seurat object with the raw (non-normalized data).
  obj <- CreateSeuratObject(counts = obj.data, project = "seurat_proj",
                            min.cells = 3, min.features = 200)
  obj

  # The [[ operator can add columns to object metadata.
  # This is a great place to stash QC stats
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mtPattern)

  #Filtering cells
  metaData <- obj@meta.data
  head(metaData)

  median_value <- median(metaData[["nFeature_RNA"]], na.rm = TRUE)
  mad_value <- mad(metaData[["nFeature_RNA"]], na.rm = TRUE)
  # Define lower and upper bounds for nFeature_RNA
  lower_bound_nF <- median_value - 3 * mad_value
  upper_bound_nF <- median_value + 3 * mad_value

  #filtering for % of Mt reads
  median_value_mt <- median(metaData[["percent.mt"]], na.rm = TRUE)
  mad_value_mt <- mad(metaData[["percent.mt"]], na.rm = TRUE)
  # Define lower and upper bounds for percent.mt
  upper_bound_mt <- median_value_mt + 4 * mad_value_mt

  #Filter out outlier cells
  obj <- subset(obj, subset = nFeature_RNA > lower_bound_nF &
                  nFeature_RNA < upper_bound_nF &
                  percent.mt < upper_bound_mt)

  #Normalization
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000)

  #Highly variable genes
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

  #Scaling
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)

  #PCA
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))

  #selection of optimal PCs
  pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  nPCs <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
               decreasing = TRUE)[1] + 1

  #Not suggested to have less than 5PCs
  if (nPCs < 5)
    nPCs <- 5

  #Clustering
  obj <- FindNeighbors(obj, dims = 1:nPCs)
  obj <- FindClusters(obj, resolution = seq(from = 0.1, to = 2, by = 0.1))

  #UMAP visualization
  obj <- RunUMAP(obj, dims = 1:nPCs)

  #save output file
  saveRDS(obj, paste0(outputName, ".RDS"))

  #Finished!
  print("Seurat analysis completed successfully!")
}

seuratAnalysis(input_Dir, output_Name, mt_Pattern)
