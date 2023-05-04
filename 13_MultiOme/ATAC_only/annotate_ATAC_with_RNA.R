library(Seurat)
library(Signac)
library(patchwork)
library(ggplot2)

projdir = '/nas/homes/benyang/HiC/13_MultiOme/ATAC_only'

integrated.ATAC = readRDS(file.path(projdir, "integrated_ATAC.RDS"))
integrated.RNA = readRDS(file.path(projdir, "JAL_RNA", "aging_MuSC_integrated.RDS"))

DefaultAssay(integrated.ATAC) = "GeneActivity"
DefaultAssay(integrated.RNA) = "integrated"

integrated.RNA = FindVariableFeatures(integrated.RNA, selection.method = "vst", nfeatures = 2000)
integrated.geneactivities = GeneActivity(integrated.ATAC, assay="ATAC", process_n=5e4)
integrated.ATAC[['IntegratedGeneActivity']] = CreateAssayObject(integrated.geneactivities)
DefaultAssay(integrated.ATAC) = "IntegratedGeneActivity"
integrated.ATAC = NormalizeData(integrated.ATAC) %>% ScaleData(features=rownames(integrated.ATAC))

transfer.anchors <- FindTransferAnchors(reference = integrated.RNA, query = integrated.ATAC, features = VariableFeatures(object = integrated.RNA), reference.assay = "integrated", query.assay = "IntegratedGeneActivity", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = integrated.RNA$celltype,
    weight.reduction = integrated.ATAC[["integrated_lsi"]], dims = 2:40)

integrated.ATAC <- AddMetaData(integrated.ATAC, metadata = celltype.predictions)

DimPlot(integrated.ATAC, group.by="predicted.id")

plt_list = lapply(paste0("prediction.score.",c("MuSC","NMJ","Immune","Smooth.Muscle","MSC","Tenocyte","Endothelial","Myonuclei")), function(x) FeaturePlot(integrated.ATAC, features=x))
wrap_plots(plt_list)
ggsave(file.path(projdir,"Figures","ATAC_transfer_RNA_predicted_scores.png"))