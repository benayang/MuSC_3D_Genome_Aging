# hicDir="/nas/homes/benyang/HiC/02_HIC"

# hicConvertFormat -m "${hicDir}/young.merged/young.merged.hic" \
# -o /nas/homes/benyang/HiC/08_HiCExplorer/normalized_matrices/young.merged.cool \
# --inputFormat hic \
# --outputFormat cool \
# --resolutions 40000

# hicConvertFormat -m "${hicDir}/aged.merged/aged.merged.40kb.hic" \
# -o /nas/homes/benyang/HiC/08_HiCExplorer/normalized_matrices/aged.merged.cool \
# --inputFormat hic \
# --outputFormat cool \
# --resolutions 40000

hicNormalize -m young.merged_40000.cool aged.merged_40000.cool \
--normalize smallest \
--outFileName young.merged_40000.normalized.cool aged.merged_40000.normalized.cool