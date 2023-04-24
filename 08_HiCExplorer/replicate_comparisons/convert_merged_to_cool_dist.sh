hicDir="/home/benjy/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C/02_HIC"

# get raw data from hic file
#hicConvertFormat -m "${hicDir}/young.merged.hic" -o young.merged_250kb.cool --inputFormat hic --outputFormat cool --resolutions 250000
#hicConvertFormat -m "${hicDir}/aged.merged_paired_alignment_dedup.fixed.hic" -o aged.merged_250kb.cool --inputFormat hic --outputFormat cool --resolutions 250000

# apply KR correction
#hicConvertFormat -m young.merged_250kb_250000.cool --inputFormat cool --outputFormat cool -o young.merged_250kb_KR.cool --correction_name KR
#hicConvertFormat -m aged.merged_250kb_250000.cool --inputFormat cool --outputFormat cool -o aged.merged_250kb_KR.cool --correction_name KR

hicPlotDistVsCounts --matrices 4140-KS-2_250kb_250000_KR.cool \
4140-KS-3_250kb_250000_KR.cool \
4140-KS-1_250kb_250000_KR.cool \
4504-KS-1_250kb_250000_KR.cool \
young.merged_250kb_KR.cool \
aged.merged_250kb_KR.cool \
--labels 'young.r1' 'young.r2' 'aged.r1' 'aged.r2' 'young.merged' 'aged.merged' \
--maxdepth 30000000 \
--plotFile ./young_aged_reps_250kb_KR_dist.png \
--outFileData ./young_aged_reps_250kb_KR_dist.txt \
--plotsize 5 4.2