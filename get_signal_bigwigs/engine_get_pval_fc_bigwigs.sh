#conda activate encode-chip-seq-pipeline-python2
# bedClip --> http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

sample_tagAlign=$1
ctrl_tagAlign=$2
fraglen=$3
outdir=$4

prefix=/nas/homes/benyang/HiC/HistoneTracks
chrom_sizes=/nas/homes/benyang/JC_H3K27me3/encode_genome_data/mm10_no_alt.chrom.sizes.tsv
bedClip=/nas/homes/benyang/HiC/bedClip

mkdir $outdir

# get fraglen from call-fraglen_mean
# Tracks from pooled Aged H3K27me3
~/chip-seq-pipeline2/src/encode_task_macs2_signal_track_chip.py $sample_tagAlign $ctrl_tagAlign \
--gensz mm --chrsz $chrom_sizes \
--fraglen $fraglen --pval-thresh 0.01 --mem-gb 30 \
--out-dir $outdir \
--log-level INFO 2> $outdir/stderr.txt 1> $outdir/stdout.txt

bedtools slop -i $outdir/rep.pooled_x_ctl.pooled_FE.bdg -g $chrom_sizes -b 0 | awk '{if ($3 != -1) print $0}' | $bedClip stdin $chrom_sizes "$outdir/rep.pooled_x_ctl.pooled.fc.signal.bedgraph"

LC_COLLATE=C sort -k1,1 -k2,2n -S 10000M $outdir/rep.pooled_x_ctl.pooled.fc.signal.bedgraph | awk 'BEGIN{OFS="\t"}{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > $outdir/rep.pooled_x_ctl.pooled.fc.signal.srt.bedgraph

rm -f $outdir/rep.pooled_x_ctl.pooled.fc.signal.bedgraph

bedGraphToBigWig $outdir/rep.pooled_x_ctl.pooled.fc.signal.srt.bedgraph $chrom_sizes $outdir/rep.pooled_x_ctl.pooled.fc.signal.bigwig

rm -f $outdir/rep.pooled_x_ctl.pooled.fc.signal.srt.bedgraph

scale_factor=`zcat $sample_tagAlign -f  | wc -l`

macs2 bdgcmp -t "$outdir/rep.pooled_x_ctl.pooled_treat_pileup.bdg" -c "$outdir/rep.pooled_x_ctl.pooled_control_lambda.bdg" --o-prefix $outdir/rep.pooled_x_ctl.pooled -m ppois -S $(awk -v s=$scale_factor 'BEGIN{print s/1e6}') 2>> $outdir/stderr.txt 1>> $outdir/stdout.txt

bedtools slop -i "$outdir/rep.pooled_x_ctl.pooled_ppois.bdg" -g $chrom_sizes -b 0 | awk '{if ($3 != -1) print $0}' | $bedClip stdin $chrom_sizes "$outdir/rep.pooled_x_ctl.pooled.pval.signal.bedgraph"

LC_COLLATE=C sort -k1,1 -k2,2n -S 10000M $outdir/rep.pooled_x_ctl.pooled.pval.signal.bedgraph | awk 'BEGIN{OFS="\t"}{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > $outdir/rep.pooled_x_ctl.pooled.pval.signal.srt.bedgraph

rm -f $outdir/rep.pooled_x_ctl.pooled.pval.signal.bedgraph

bedGraphToBigWig $outdir/rep.pooled_x_ctl.pooled.pval.signal.srt.bedgraph $chrom_sizes $outdir/rep.pooled_x_ctl.pooled.pval.signal.bigwig

rm -f $outdir/rep.pooled_x_ctl.pooled.pval.signal.srt.bedgraph

#rm -f rep.pooled_x__1_1_NoAbYoung__2_1_NoAbYoung__3_1_NoAbYoung__R1.fastq.srt.nodup_*