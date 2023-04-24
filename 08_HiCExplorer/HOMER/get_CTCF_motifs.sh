mainDir=/nas/homes/benyang/HiC/08_HiCExplorer

findMotifsGenome.pl \
$mainDir/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
mm10 \
$mainDir/homer/young.merged \
-bg $mainDir/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
-size 200 \
-p 50

findMotifsGenome.pl \
$mainDir/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
mm10 \
$mainDir/homer/aged.merged \
-bg $mainDir/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
-size 200 \
-p 50 > $mainDir/homer/aged.merged_40kb_0.01FDR.bed


