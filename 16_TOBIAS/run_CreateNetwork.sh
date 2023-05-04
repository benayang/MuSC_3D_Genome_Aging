find /nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/ -name *_bound.bed > "/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/TFBS_bound_file_list.txt"

cp $(cat merged_peakset_annotated_expressed_TF/TFBS_bound_file_list.txt | grep "aged") merged_peakset_annotated_expressed_TF_network/bound_files/aged
cp $(cat merged_peakset_annotated_expressed_TF/TFBS_bound_file_list.txt | grep "young") merged_peakset_annotated_expressed_TF_network/bound_files/young

TOBIAS CreateNetwork --TFBS \
/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF_network/bound_files/young/* \
--origin JASPAR_to_Ensembl.txt --outdir merged_peakset_annotated_expressed_TF_network/young

TOBIAS CreateNetwork --TFBS \
/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF_network/bound_files/aged/* \
--origin JASPAR_to_Ensembl.txt --outdir merged_peakset_annotated_expressed_TF_network/aged

# TOBIAS CreateNetwork --TFBS \
# "/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/CTCF_MA0139.1/beds/CTCF_MA0139.1_young_merged_ATAC_bound.bed" \
# "/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/CTCF_MA1929.1/beds/CTCF_MA1929.1_young_merged_ATAC_bound.bed" \
# "/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/CTCF_MA1930.1/beds/CTCF_MA1930.1_young_merged_ATAC_bound.bed" \
# --origin JASPAR_to_Ensembl.txt --outdir merged_peakset_annotated_expressed_TF_network/CTCF/young

# TOBIAS CreateNetwork --TFBS \
# "/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/CTCF_MA0139.1/beds/CTCF_MA0139.1_aged_merged_ATAC_bound.bed" \
# "/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/CTCF_MA1929.1/beds/CTCF_MA1929.1_aged_merged_ATAC_bound.bed" \
# "/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/CTCF_MA1930.1/beds/CTCF_MA1930.1_aged_merged_ATAC_bound.bed" \
# --origin JASPAR_to_Ensembl.txt --outdir merged_peakset_annotated_expressed_TF_network/CTCF/aged

# TOBIAS CreateNetwork --TFBS \
# "/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/MYOD1_MA0499.2/beds/MYOD1_MA0499.2_aged_merged_ATAC_bound.bed" \
# --origin JASPAR_to_Ensembl.txt --outdir merged_peakset_annotated_expressed_TF_network/Myod1/aged

# TOBIAS CreateNetwork --TFBS \
# "/nas/homes/benyang/HiC/16_TOBIAS/merged_peakset_annotated_expressed_TF/MYOD1_MA0499.2/beds/MYOD1_MA0499.2_young_merged_ATAC_bound.bed" \
# --origin JASPAR_to_Ensembl.txt --outdir merged_peakset_annotated_expressed_TF_network/Myod1/young
