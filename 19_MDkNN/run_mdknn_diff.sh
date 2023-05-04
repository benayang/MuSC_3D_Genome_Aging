# conda deactivate
# conda activate py36

cool_dir="/nas/homes/benyang/HiC/08_HiCExplorer"
outdir="/nas/homes/benyang/HiC/19_MDkNN"

for age in young aged
do

    mkdir $outdir/${age}_TAD_diff

    python ~/MDkNN/mdknn.py \
    -p $cool_dir/aged.merged_10kb_KR.cool \
    -p2 $cool_dir/young.merged_10kb_KR.cool \
    -t "/nas/homes/benyang/HiC/08_HiCExplorer/$age.merged/40kb/$age.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed" \
    -O $outdir/${age}_TAD_diff/MDkNN_${age}_TAD_diff_10kb.txt \
    --out-long-range $outdir/${age}_TAD_diff/MDkNN_${age}_TAD_diff_10kb_longrange.txt \
    -k 3 --ww 5 --pw 2

done