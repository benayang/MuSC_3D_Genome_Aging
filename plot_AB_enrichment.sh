prefix=/nas/homes/benyang/HiC

#fanc compartments $prefix/02_HIC/aged.merged/aged.merged.hic@100kb@KR \
#$prefix/04_FANC/without_KR_normalization/aged.merged_100kb_noKR.ab \
fanc compartments $prefix/08_HiCExplorer/aged.merged_100kb_KR.cool \
$prefix/04_FANC/without_KR_normalization/aged.merged_100kb_noKR_noNAN.ev \
-g "/nas/homes/benyang/JC_H3K27me3/encode_genome_data/mm10_no_alt_analysis_set_ENCODE.fasta.gz" \
-m $prefix/04_FANC/without_KR_normalization/AB_enrichment_profiles/cool_aged_merged_100kb_KR_ab_profile.mat \
--compartment-strength $prefix/04_FANC/without_KR_normalization/AB_enrichment_profiles/cool_aged_merged_100kb_KR_ab_strength.mat \
-e cool_aged_merged_100kb_KR_ab_profile.png \
-x chrM chrY \
-s 0 \
-f \
--enrichment-min -1 \
--enrichment-max 1 \
-p 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100

#fanc compartments $prefix/02_HIC/young.merged/young.merged.hic@100kb@KR \
#$prefix/04_FANC/without_KR_normalization/young.merged_100kb_noKR.ab \
fanc compartments $prefix/08_HiCExplorer/young.merged_100kb_KR.cool \
$prefix/04_FANC/without_KR_normalization/young.merged_100kb_noKR_noNAN.ev \
-g "/nas/homes/benyang/JC_H3K27me3/encode_genome_data/mm10_no_alt_analysis_set_ENCODE.fasta.gz" \
-m $prefix/04_FANC/without_KR_normalization/AB_enrichment_profiles/cool_young_merged_100kb_KR_ab_profile.mat \
--compartment-strength $prefix/04_FANC/without_KR_normalization/AB_enrichment_profiles/cool_young_merged_100kb_KR_ab_strength.mat \
-e cool_young_merged_100kb_KR_ab_profile.png \
-s 0 \
-f \
-x chrM chrY \
--enrichment-min -1 \
--enrichment-max 1 \
-p 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100