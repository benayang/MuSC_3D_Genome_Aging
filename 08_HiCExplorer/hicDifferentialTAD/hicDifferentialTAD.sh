prefix='/nas/homes/benyang/HiC/08_HiCExplorer'

hicDifferentialTAD --targetMatrix $prefix/hicDifferentialTAD/aged.merged_40000_minNorm_KR.cool \
--controlMatrix $prefix/hicDifferentialTAD/young.merged_40000_minNorm_KR.cool \
--tadDomains $prefix/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
--outFileNamePrefix $prefix/hicDifferentialTAD/aged_allMode_differential_TAD \
--pValue 0.05 \
--mode all \
--modeReject all \
--threads 30

hicDifferentialTAD --targetMatrix $prefix/hicDifferentialTAD/young.merged_40000_minNorm_KR.cool \
--controlMatrix $prefix/hicDifferentialTAD/aged.merged_40000_minNorm_KR.cool \
--tadDomains $prefix/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
--outFileNamePrefix $prefix/hicDifferentialTAD/young_allMode_differential_TAD \
--pValue 0.05 \
--mode all \
--modeReject all \
--threads 30