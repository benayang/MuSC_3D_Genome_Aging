library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ArchR)
library(future)

projdir = "/nas/homes/benyang/HiC/13_MultiOme/"

# Make QC plots ------------------------------------------------------
aged_rep5 = readRDS(file.path(projdir,"JAL_RNA","seurat_obj","aged_rep5.RDS"))
young_rep3 = readRDS(file.path(projdir,"JAL_RNA","seurat_obj","young_rep3.RDS"))

tms_3mo = readRDS(file.path(projdir,"TabulaMurisSenis","ADJUSTED_limb_muscle_3mo_droplet_Seurat.RDS"))
tms_24mo = readRDS(file.path(projdir,"TabulaMurisSenis","ADJUSTED_limb_muscle_24mo_droplet_Seurat.RDS"))
tms_3mo$complexity = log10(tms_3mo$nFeature_RNA)/log10(tms_3mo$nCount_RNA)
tms_24mo$complexity = log10(tms_24mo$nFeature_RNA)/log10(tms_24mo$nCount_RNA)
tms_3mo$orig.ident = "TMS_3m"
tms_24mo$orig.ident = "TMS_24m"

md_cols = c("orig.ident","nCount_RNA","nFeature_RNA","complexity")
all_qc_md = rbind(aged_rep5@meta.data[,md_cols],
                    young_rep3@meta.data[,md_cols],
                    tms_3mo@meta.data[,md_cols],
                    tms_24mo@meta.data[,md_cols])

plotlist = list()
for(i in unique(all_qc_md$orig.ident)) {
    plotlist[[i]] = ggPoint(
    x = log10(all_qc_md$nCount_RNA[all_qc_md$orig.ident==i]), 
    y = all_qc_md$nFeature_RNA[all_qc_md$orig.ident==i], 
    colorDensity = TRUE,
    xlim = c(2.5, 5),
    ylim = c(0, 7000),
    continuousSet = "sambaNight",
    xlabel = "Log10 # UMIs",
    ylabel = "# Genes"
) + geom_hline(yintercept = 500, lty = "dashed") +
    geom_hline(yintercept = 4000, lty = "dashed") +
    geom_vline(xintercept = log10(1000), lty = "dashed")

    png(file.path(projdir, "TabulaMurisSenis", paste0(i, "_qc.png")), res=300, units='in', width=5, height=5)
    print(plotlist[[i]])
    dev.off()
}

# filter objects ------------------------------------------------------
aged_rep5 = subset(aged_rep5, nCount_RNA>500 & nFeature_RNA>300 & nFeature_RNA<4000 & complexity>0.8 & percent.mt<10)
young_rep3 = subset(young_rep3, nCount_RNA>500 & nFeature_RNA>300 & nFeature_RNA<4000 & complexity>0.8 & percent.mt<10)

# preprocess datasets ------------------------------------------------------
obj_list = list(aged_rep5=aged_rep5, young_rep3=young_rep3, tms_3mo=tms_3mo, tms_24mo=tms_24mo)
#names(obj_list) = c("aged_rep1", "aged_rep2", "aged_rep5", "aged_rep6", "young_rep1", "young_rep3")

obj_list = lapply(obj_list, function(x){
  x = x %>%
    NormalizeData(normalization.method="LogNormalize", scale.factor=10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA()
})

qc_plts = lapply(1:length(obj_list), function(x) FeaturePlot(obj_list[[x]], features=c("percent.mt","nCount_RNA","nFeature_RNA","complexity"), order=T))

obj_list[c("aged_rep5","young_rep3")] = lapply(obj_list[c("aged_rep5","young_rep3")], function(x){
  x = x %>%
    NormalizeData(normalization.method="LogNormalize", scale.factor=10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(vars.to.regress = c("percent.mt","nCount_RNA")) %>%
    RunPCA()
})

obj_list[c("tms_3mo","tms_24mo")] = lapply(obj_list[c("tms_3mo","tms_24mo")], function(x){
  x = x %>%
    NormalizeData(normalization.method="LogNormalize", scale.factor=10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(vars.to.regress = "nCount_RNA") %>%
    RunPCA()
})

# Dimensionality of datasets  ---------------------------------------------

quant_pc_elbow = function(data, reduction) {
  # Determine percent of variation associated with each PC
  pct <- data[[reduction]]@stdev / sum(data[[reduction]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  co2
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  # Create a dataframe with values
  plot_df <- data.frame(pct = pct, 
                        cumu = cumu, 
                        rank = 1:length(pct))
  # Elbow plot to visualize 
  ggplot(plot_df, aes(x=cumu, y=pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    ggtitle(deparse(substitute(data))) +
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw()
}

quant_pc_elbow(obj_list[['aged_rep5']], "pca") # 21
quant_pc_elbow(obj_list[['young_rep3']], "pca") # 13
quant_pc_elbow(obj_list[['tms_3mo']], "pca") # 29
quant_pc_elbow(obj_list[['tms_24mo']], "pca") # 34

# Embedding and clustering  ---------------------------------------------

dim_list = c(21, 13, 29, 34)
obj_list = lapply(seq_along(obj_list), 
                  function(x) {
                    obj_list[[x]] = obj_list[[x]] %>%
                      RunUMAP(dims = seq_len(dim_list[x]), n.neighbors=50, min.dist=0.5, reduction.name="umap.rna", reduction.key="rnaUMAP_") %>%
                      FindNeighbors(reduction="pca", dims=seq_len(dim_list[x])) %>%
                      FindClusters(algorithm=1, resolution=0.1)
                  })

DimPlot(obj_list[[1]], label=T)
DimPlot(obj_list[[2]], label=T)
DimPlot(obj_list[[3]], label=T)
DimPlot(obj_list[[4]], label=T)

# Integrate datasets  ---------------------------------------------

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures=2000)
aging.anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
aging_RNA <- IntegrateData(anchorset = aging.anchors, normalization.method = "LogNormalize")

DefaultAssay(aging_RNA) = "integrated"

aging_RNA = aging_RNA %>%
  ScaleData() %>%
  RunPCA(npcs=50)

quant_pc_elbow(aging_RNA, "pca") # 33

aging_RNA = aging_RNA %>%
                RunUMAP(dims = 1:33, n.neighbors=30, min.dist=0.5, reduction.name="umap.rna", reduction.key="rnaUMAP_") %>%
                FindNeighbors(reduction="pca", dims=1:33)

# Clustering silhouette analysis  ---------------------------------------------

for(res in seq(0.1,1.2,0.1)) {
  aging_RNA = FindClusters(aging_RNA, graph.name="integrated_snn", algorithm=1, resolution=res)
}

library(cluster, quietly = TRUE)
silhouette_analysis = function(seurat_subset, reduction, clusters, dims) {
  tmp = list()
  dist.matrix <- dist(x = seurat_subset[[reduction]]@cell.embeddings[ ,dims])
  for(res in clusters) {
    sil <- silhouette(x = as.numeric(x = as.factor(x = unlist(seurat_subset[[paste0("integrated_snn_res.", res)]]))), 
                      dist = dist.matrix)
    tmp[[paste0("res_",res)]] = sil[ ,3]
  }
  return(tmp)  
}

aging_RNA_silhouette = silhouette_analysis(aging_RNA, "pca", seq(0.1,1.2,0.1), 1:33)

png(file.path(projdir,"JAL_RNA","QC_plots","aging_RNA_silhouette_boxplot.png"), res=300, units="in", width=4, height=4)
boxplot(aging_RNA_silhouette, xlab="Clustering Resolution", ylab="Silhouette Score")
dev.off()

saveRDS(aging_RNA, "/nas/homes/benyang/HiC/13_MultiOme/TabulaMurisSenis/aging_RNA_integrated.RDS")
saveRDS(obj_list, "/nas/homes/benyang/HiC/13_MultiOme/TabulaMurisSenis/seurat_obj_list.RDS")

# Cell cycle scoring ------------------------------------------------------

cc.genes = read.table(file.path(projdir,"TabulaMurisSenis","hbc_tinyatlas_cell_cycle_Mus_musculus.csv"), sep=",", header=T)
# data(cc.genes.updated.2019)
s.genes = cc.genes$gprofiler[cc.genes$phase=="S"]
g2m.genes = cc.genes$gprofiler[cc.genes$phase=="G2/M"]
aging_RNA = CellCycleScoring(aging_RNA, assay="RNA", s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(aging_RNA, reduction = "pca", group.by="Phase")

# Verification of clustering ------------------------------------------------------

DefaultAssay(aging_RNA) = "RNA"

plan("multicore", workers=30)

# Look at marker genes from Tabula Muris Senis dataset
lm_droplet = readRDS(file.path(projdir,"TabulaMurisSenis","ADJUSTED_limb_muscle_droplet_raw_data.RDS"))
Idents(lm_droplet) = "age"
lm_droplet_subset = subset(lm_droplet, idents=c("3m","24m"))
Idents(lm_droplet_subset) = "cell_ontology_class"
lm_droplet_markers = FindAllMarkers(lm_droplet_subset, logfc.threshold = 0.25, test.use="wilcox", min.pct=0.1)
saveRDS(lm_droplet_markers, file.path(projdir, "TabulaMurisSenis", "TMS_subset_all_markers.RDS"))

tcell.markers <- FindConservedMarkers(lm_droplet_subset, ident.1 = "T cell", grouping.var = "age", logfc.threshold=0.25, min.pct=0.1, verbose = FALSE)

Idents(aging_RNA) = "integrated_snn_res.0.1"
all_markers = FindAllMarkers(aging_RNA, logfc.threshold = 0.25, test.use="wilcox", min.pct=0.1)
saveRDS(all_markers, file.path(projdir,"TabulaMurisSenis","all_markers_integrated.RDS"))

png(file.path(projdir,"TabulaMurisSenis","integrated_clusters_UMAP.png"), res=300, units='in', width=5, height=5)
DimPlot(aging_RNA, group="integrated_snn_res.0.1", label=T) +
    labs(x="UMAP 1", y="UMAP 2") +
    coord_equal()
dev.off()

png(file.path(projdir,"TabulaMurisSenis","integrated_clusters_split_UMAP.png"), res=300, units='in', width=15, height=5)
DimPlot(aging_RNA, group="integrated_snn_res.0.1", split="orig.ident", label=T) +
    labs(x="UMAP 1", y="UMAP 2") +
    coord_equal()
dev.off()

# 0 - MuSC
# 1 - MSC
# 2 - Endothelial
# 3 - Schwann cells
# 4 - Tenocyte
# 5 - Immune (Macrophage)
# 6 - Immune (B Cells)
# 7 - Smooth Muscle
# 8 - Immune (NK Cells)
# 9 - Immune (T Cells)
# 10 - More Endothelial (probably JAL MuSC)
# 11 - More MSCs
# 12 - Skeletal muscle

Idents(aging_RNA) = "integrated_snn_res.0.1"
immune.markers <- FindConservedMarkers(aging_RNA, ident.1 = 8, ident.2 = 9, grouping.var = "orig.ident", logfc.threshold=0.25, min.pct=0.1, verbose = FALSE)
immune.markers[order(immune.markers$TMS_24m_avg_log2FC, decreasing=T), ]

c5.markers <- FindConservedMarkers(aging_RNA, ident.1 = 5, grouping.var = "orig.ident", logfc.threshold=0.25, min.pct=0.1, verbose = FALSE)
c10.markers <- FindConservedMarkers(aging_RNA, ident.1 = 10, grouping.var = "orig.ident", logfc.threshold=0.25, min.pct=0.1, verbose = FALSE)
c11.markers <- FindConservedMarkers(aging_RNA, ident.1 = 11, grouping.var = "orig.ident", logfc.threshold=0.25, min.pct=0.1, verbose = FALSE)

# Verification of clustering ------------------------------------------------------

aging_RNA$celltype = aging_RNA$integrated_snn_res.0.1
Idents(aging_RNA) = "celltype"
aging_RNA = RenameIdents(aging_RNA, '0'="MuSC")
aging_RNA = RenameIdents(aging_RNA, '1'="MSC")
aging_RNA = RenameIdents(aging_RNA, '2'="Endothelial")
aging_RNA = RenameIdents(aging_RNA, '3'="Schwann Cell")
aging_RNA = RenameIdents(aging_RNA, '4'="Tenocyte")
aging_RNA = RenameIdents(aging_RNA, '5'="Macrophage")
aging_RNA = RenameIdents(aging_RNA, '6'="B Cell")
aging_RNA = RenameIdents(aging_RNA, '7'="Smooth Muscle")
aging_RNA = RenameIdents(aging_RNA, '8'="NK Cell")
aging_RNA = RenameIdents(aging_RNA, '9'="T Cell")
aging_RNA = RenameIdents(aging_RNA, '10'="Endothelial")
aging_RNA = RenameIdents(aging_RNA, '11'="MSC")
aging_RNA = RenameIdents(aging_RNA, '12'="Skeletal Muscle")
aging_RNA$celltype = Idents(aging_RNA)
aging_RNA$celltype = factor(aging_RNA$celltype, 
                        levels=rev(c("MuSC","Skeletal Muscle","Smooth Muscle","MSC","Endothelial","Schwann Cell","Tenocyte","Macrophage","B Cell","NK Cell","T Cell")))

marker_genes = c(
    "Pax7", "Myod1", "Myf5",
    "Myoz1", "Ckm", "Myl1",
    "Acta2", "Myl9", "Myh11",
    "Pdgfra", "Smoc2", "Gsn",
    "Fabp4", "Cdh5", "Pecam1",
    "Plp1", "Ptn", "Mpz", 
    "Fmod", "Thbs4", "Tnmd",
    "Fcer1g", "Lyz1", "C1qb",
    "Cd79a", "Cd79b", "Ms4a1",
    "Cd8b1", "Nkg7", "Ccl5",
    "Il7r", "Cxcr6", "Cd3g")

png(file.path(projdir,"TabulaMurisSenis","integrated_clusters_marker_genes.png"), res=300, units="in", width=10, height=4)
DotPlot(aging_RNA, group.by="celltype", features=marker_genes) + RotatedAxis() + labs(x=NULL, y=NULL)
dev.off()

saveRDS(aging_RNA, file.path(projdir,"TabulaMurisSenis","aging_RNA_integrated_annotated.RDS"))

# Make dim plot ---------------------------------------

dimplot <- DimPlot(aging_RNA, split.by="orig.ident", group.by="celltype", ncol=2) + 
  scale_color_manual(values=DiscretePalette(n=length(unique(aging_RNA$celltype)), "glasbey")) +
  labs(x=NULL, y=NULL, title=NULL)
dimplot$data$orig.ident <- factor(dimplot$data$orig.ident, levels=c("young_rep3","aged_rep5","TMS_3m","TMS_24m"))
levels(dimplot$data$orig.ident) <- list("Young MuSC" = "young_rep3", "Aged MuSC" = "aged_rep5", "TMS (3m)" = "TMS_3m", "TMS (24m)" = "TMS_24m")

png(file.path(projdir,"TabulaMurisSenis","integrated_celltype_dimplot.png"), res=300, units="in", width=6, height=5)
print(dimplot)
dev.off()