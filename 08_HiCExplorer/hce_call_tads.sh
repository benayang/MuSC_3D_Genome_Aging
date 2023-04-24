mainDir="/nas/homes/benyang/HiC/08_HiCExplorer"
sample="aged.merged"

hicFindTADs --matrix "$mainDir/${sample}_100kb_KR.cool" \
--outPrefix "$mainDir/${sample}/100kb/${sample}_100kb_min300kb_max3mb_step100kb_thresh0.01_delta_0.01_fdr" \
--minDepth 300000 \
--maxDepth 3000000 \
--step 100000 \
--thresholdComparisons 0.01 \
--delta 0.01 \
--correctForMultipleTesting fdr \
--chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX \
-p 45

for t in 0.1 0.05 0.005 0.001
do
hicFindTADs --matrix "$mainDir/${sample}_100kb_KR.cool" \
--outPrefix "$mainDir/${sample}/100kb/${sample}_100kb_min300kb_max3mb_step100kb_thresh${t}_delta_0.01_fdr" \
--TAD_sep_score_prefix "$mainDir/${sample}/100kb/${sample}_100kb_min300kb_max3mb_step100kb_thresh0.01_delta_0.01_fdr" \
--thresholdComparisons $t \
--delta 0.01 \
--correctForMultipleTesting fdr \
--chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX \
-p 45
done
