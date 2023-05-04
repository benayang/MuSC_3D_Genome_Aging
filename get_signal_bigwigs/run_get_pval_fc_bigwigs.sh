chip_prefix=/nas/homes/benyang/HiC/HistoneTracks/Tom_Rando/chip
croo_prefix=/nas/homes/benyang/HiC/HistoneTracks/Tom_Rando/croo
prefix=/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/macs2_signal

# Young H3K27me3
# ./engine_get_pval_fc_bigwigs.sh \
# $croo_prefix/croo_Young_H3K27me3/align/pooled-rep/rep.pooled.tagAlign.gz \
# $croo_prefix/croo_Young_H3K27me3/align/pooled-ctl/ctl.pooled.tagAlign.gz \
# `cat $chip_prefix/5dcd131e-e292-41ad-a145-3d7453ce50fb/call-fraglen_mean/execution/tmp.txt` \
# $prefix/Young_H3K27me3

# Young H3K4me3
# ./engine_get_pval_fc_bigwigs.sh \
# $croo_prefix/croo_Young_H3K4me3/align/pooled-rep/rep.pooled.tagAlign.gz \
# $croo_prefix/croo_Young_H3K4me3/align/pooled-ctl/ctl.pooled.tagAlign.gz \
# `cat $chip_prefix/151f676b-0fcb-4f13-929f-0893582c197d/call-fraglen_mean/execution/tmp.txt` \
# $prefix/Young_H3K4me3

# Aged H3K4me3
# ./engine_get_pval_fc_bigwigs.sh \
# $croo_prefix/croo_Aged_H3K4me3/align/pooled-rep/rep.pooled.tagAlign.gz \
# $croo_prefix/croo_Aged_H3K4me3/align/pooled-ctl/ctl.pooled.tagAlign.gz \
# `cat $chip_prefix/cc20afd8-1343-41eb-99d1-c392accf6c94/call-fraglen_mean/execution/tmp.txt` \
# $prefix/Aged_H3K4me3

# H3K27ac
./engine_get_pval_fc_bigwigs.sh \
"/nas/homes/benyang/HiC/HistoneTracks/H3K27ac/chip/eafea56b-4971-4718-aaa9-b37b24febaba/call-pool_ta/execution/rep.pooled.tagAlign.gz" \
"/nas/homes/benyang/HiC/HistoneTracks/H3K27ac/chip/eafea56b-4971-4718-aaa9-b37b24febaba/call-pool_ta_ctl/execution/ctl.pooled.tagAlign.gz" \
`cat "/nas/homes/benyang/HiC/HistoneTracks/H3K27ac/chip/eafea56b-4971-4718-aaa9-b37b24febaba/call-fraglen_mean/execution/tmp.txt"` \
$prefix/H3K27ac