# Seurat tutorial with example 10X PBMC data
# http://satijalab.org/seurat/pbmc-tutorial.html

source("https://bioconductor.org/biocLite.R")
biocLite(c("import", "tidyverse", "satijalab/seurat"))

library(tidyverse)
library(Seurat)
import::from(Matrix, colSums)
import::from(plyr, mapvalues)

download.file("https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz", "pbmc.tar.gz")

untar("pbmc.tar.gz")
pbmc.data <- Read10X("filtered_gene_bc_matrices/hg19/")

dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size

pbmc <- new("seurat", raw.data = pbmc.data)
pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, do.logNormalize = TRUE,
              total.expr = 1e4, project = "10X_PBMC")

mito.genes <- grep("^MT-", rownames(pbmc@data), value = TRUE)
percent.mito <- colSums(expm1(pbmc@data[mito.genes, ])) /
    colSums(expm1(pbmc@data))
pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
VlnPlot(pbmc, c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(pbmc, "nUMI", "percent.mito")
GenePlot(pbmc, "nUMI", "nGene")

pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = 2500)
pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.05)

# Regress, takes about 60 seconds
pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))

pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean,
                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,
                    do.contour = FALSE)
length(pbmc@var.genes)

pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = FALSE,
            pcs.print = 5, genes.print = 5)
pbmc <- ProjectPCA(pbmc, do.print = FALSE)
PrintPCA(pbmc, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
VizPCA(pbmc, 1:2)
PCAPlot(pbmc, 1, 2)

PCHeatmap(pbmc, pc.use = 1, cells.use = 100, do.balanced = TRUE)
PCHeatmap(pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE,
          label.columns = FALSE, use.full = FALSE)

# pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = FALSE)
# JackStrawPlot(pbmc, PCs = 1:12)

PCElbowPlot(pbmc)
pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = 0.6, save.SNN = TRUE)

pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(pbmc)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
print(head(cluster1.markers, 5))

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25,
                               thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_diff)

VlnPlot(pbmc, c("MS4A1","CD79A"))
VlnPlot(pbmc, c("NKG7","PF4"), use.raw = TRUE, y.log = TRUE)

FeaturePlot(pbmc, c("MS4A1", "GNLY","CD3E","CD14","FCER1A","FCGR3A", "LYZ",
                    "PPBP", "CD8A"),
            cols.use = c("grey","blue"))

pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10
DoHeatmap(pbmc, genes.use = top10$gene, order.by.ident = TRUE,
          slim.col.label = TRUE, remove.key = TRUE)

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells",
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells",
                     "Megakaryocytes")
pbmc@ident <- mapvalues(pbmc@ident,
                        from = current.cluster.ids,
                        to = new.cluster.ids)
TSNEPlot(pbmc, do.label = T, pt.size = 0.5)
