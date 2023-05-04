# conda deactivate
# conda activate py36

cool_dir="/nas/homes/benyang/HiC/08_HiCExplorer"
outdir="/nas/homes/benyang/HiC/19_MDkNN"

for age in young 
do

    mkdir $outdir/${age}_individual

    python ~/MDkNN/mdknn.py \
    -p $cool_dir/$age.merged_10kb_KR.cool \
    -t "/nas/homes/benyang/HiC/08_HiCExplorer/$age.merged/40kb/$age.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed" \
    -O $outdir/${age}_individual/MDkNN_${age}_10kb_individual.txt \
    --out-long-range $outdir/${age}_individual/MDkNN_${age}_10kb_individual_longrange.txt \
    -k 3 --ww 5 --pw 2

done