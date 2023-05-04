library(TabulaMurisSenisData)
library(Seurat)
library(dplyr)
library(tidyr)
library(rhdf5)

projdir = "/nas/homes/benyang/HiC/13_MultiOme/TabulaMurisSenis/"

# aws s3api get-object --bucket czb-tabula-muris-senis --key Data-objects/tabula-muris-senis-droplet-processed-official-annotations-Limb_Muscle.h5ad tabula-muris-senis-droplet-processed-official-annotations-Limb_Muscle.h5ad
# {
#     "AcceptRanges": "bytes",
#     "LastModified": "Thu, 05 Dec 2019 18:48:03 GMT",
#     "ContentLength": 826884650,
#     "ETag": "\"89d0eee0ccba97cebc7fb7d02ac8e4f3-99\"",
#     "ContentType": "binary/octet-stream",
#     "Metadata": {}
# }

## USE SCANPY IN PYTHON TO FIX H5AD FOR IMPORT INTO SEURAT
# import scanpy
# data = scanpy.read_h5ad("/nas/homes/benyang/HiC/13_MultiOme/TabulaMurisSenis/tabula-muris-senis-droplet-processed-official-annotations-Limb_Muscle.h5ad")
# data.write("/nas/homes/benyang/HiC/13_MultiOme/TabulaMurisSenis/ADJUSTED_tabula-muris-senis-droplet-processed-official-annotations-Limb_Muscle.h5ad")
# data.obs.to_csv("/nas/homes/benyang/HiC/13_MultiOme/TabulaMurisSenis/ADJUSTED_tabula-muris-senis-droplet-processed-official-annotations-Limb_Muscle.csv")


## CONVERT H5AD ANNDATA TO SEURAT 
SeuratDisk::Convert(file.path(projdir,"ADJUSTED_tabula-muris-senis-droplet-processed-official-annotations-Limb_Muscle.h5ad"), dest = "h5seurat", overwrite = TRUE)
lm_droplet = SeuratDisk::LoadH5Seurat(file.path(projdir,"ADJUSTED_tabula-muris-senis-droplet-processed-official-annotations-Limb_Muscle.h5seurat"), misc=F, meta.data=F)

lm_droplet_obs = read.table(file.path(projdir, "ADJUSTED_tabula-muris-senis-droplet-processed-official-annotations-Limb_Muscle.csv"), sep=",", header=T)
rownames(lm_droplet_obs) = lm_droplet_obs$index

identical(Cells(lm_droplet), lm_droplet_obs$index)

saveRDS(lm_droplet, file.path(projdir,"ADJUSTED_limb_muscle_droplet_raw_data.RDS"))

## Subset Seurat object
Idents(lm_droplet) = "age"
lm_droplet_3m = subset(lm_droplet, idents="3m")
lm_droplet_24m = subset(lm_droplet, idents="24m")

saveRDS(lm_droplet_3m, file.path(projdir,"ADJUSTED_limb_muscle_3mo_droplet_Seurat.RDS"))
saveRDS(lm_droplet_24m, file.path(projdir,"ADJUSTED_limb_muscle_24mo_droplet_Seurat.RDS"))

# TabulaMurisSenisData::listTabulaMurisSenisTissues("Droplet")
# lm_droplet = TabulaMurisSenisDroplet(
#     tissues = "Limb_Muscle",
#     processedCounts = F,
#     reducedDims = F,
#     infoOnly = F
# )[[1]]

# saveRDS(lm_droplet, file.path(projdir,"limb_muscle_droplet_raw_data.RDS"))

# lm_coldata = colData(lm_droplet)
# cellnames_3mo = rownames(lm_coldata)[lm_coldata$age == "3m"]
# cellnames_24mo = rownames(lm_coldata)[lm_coldata$age == "24m"]

# lm_droplet_3m = lm_droplet[,cellnames_3mo]
# lm_droplet_3m_counts = as.matrix(counts(lm_droplet_3m))
# lm_droplet_24m = lm_droplet[,cellnames_24mo]
# lm_droplet_24m_counts = as.matrix(counts(lm_droplet_24m))

# lm_droplet_3m_md = as.data.frame(lm_coldata)[cellnames_3mo, !(colnames(lm_coldata) %in% c("n_genes","n_counts","leiden","louvain"))]
# lm_droplet_24m_md = as.data.frame(lm_coldata)[cellnames_24mo, !(colnames(lm_coldata) %in% c("n_genes","n_counts","leiden","louvain"))]

# lm_seurat_3m = CreateSeuratObject(counts = lm_droplet_3m_counts, names.delim = "_", meta.data = lm_droplet_3m_md, project="TMS_3mo")
# lm_seurat_3m$age = as.character(lm_seurat_3m$age)
# lm_seurat_3m$orig.ident = "TMS_3mo"

# lm_seurat_24m = CreateSeuratObject(counts = lm_droplet_24m_counts, names.delim = "_", meta.data = lm_droplet_24m_md, project="TMS_24mo")
# lm_seurat_24m$age = as.character(lm_seurat_24m$age)
# lm_seurat_24m$orig.ident = "TMS_24mo"

# saveRDS(lm_seurat_3m, file.path(projdir,"limb_muscle_3mo_droplet_Seurat.RDS"))
# saveRDS(lm_seurat_24m, file.path(projdir,"limb_muscle_24mo_droplet_Seurat.RDS"))