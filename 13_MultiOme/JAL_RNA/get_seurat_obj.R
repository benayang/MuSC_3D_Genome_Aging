library(Seurat)

projdir = '/nas/homes/jlarouche/Aging/SingleCell/Seurat/10XDatasets'

# Aged
MuSC_d0_Aged_1_data <- Read10X(file.path(projdir,"MuSC_d0_Aged_1/mm10")) #v2
MuSC_d0_Aged_2_data <- Read10X(file.path(projdir,"MuSC_d0_Aged_2/mm10")) #v2
MuSC_d0_Aged_5_data <- Read10X(file.path(projdir,"MuSC_d0_Aged_5/mm10")) #v3
MuSC_d0_Aged_6_data <- Read10X(file.path(projdir,"Sample_Geriatric_MuSC_1")) #v3
# Young
MuSC_d0_Young_1_data <- Read10X(file.path(projdir,"MuSC_d0_Young_1/mm10")) #v2
MuSC_d0_Young_2_data <- Read10X(file.path(projdir,"MuSC_d0_Young_2/mm10")) #v2
MuSC_d0_Young_3_data <- Read10X(file.path(projdir,"MuSC_d0_Young_3/mm10")) #v3

data_list = list(
    aged_rep1 = MuSC_d0_Aged_1_data, 
    aged_rep2 = MuSC_d0_Aged_2_data, 
    aged_rep5 = MuSC_d0_Aged_5_data, 
    aged_rep6 = MuSC_d0_Aged_6_data, 
    young_rep1 = MuSC_d0_Young_1_data, 
    young_rep2 = MuSC_d0_Young_2_data, 
    young_rep3 = MuSC_d0_Young_3_data
    )

chemistry_list = list(
    aged_rep1 = "v2", 
    aged_rep2 = "v2", 
    aged_rep5 = "v3", 
    aged_rep6 = "v3", 
    young_rep1 = "v2", 
    young_rep2 = "v2", 
    young_rep3 = "v3"
    )

for (n in names(data_list)) {
    print(n)
    obj = CreateSeuratObject(data_list[[n]], project=n)
    obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")
    obj$chemistry = chemistry_list[[n]]
    obj$complexity = log10(obj$nFeature_RNA)/log10(obj$nCount_RNA)
    saveRDS(obj, file.path('/nas/homes/benyang/HiC/13_MultiOme/JAL_RNA/seurat_obj',paste0(n,".RDS")))
}