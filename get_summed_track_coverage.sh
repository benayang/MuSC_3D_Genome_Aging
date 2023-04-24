prefix=/nas/homes/benyang/HiC
trackPrefix=$prefix/HistoneTracks/get_signal_bigwigs
outdir=$prefix/04_FANC/track_coverage/sum_coverage

for age in young aged
do
    for f in 'A_to_B' 'B_to_A' 'StaticA' 'StaticB'
    do
        python3 get_summed_track_coverage.py --BW ${age}_bigwig_list.txt \
        --BED $prefix/04_FANC/compartmentExpression/compartmentBed/$f.bed \
        --outf $outdir/${age}.$f.sumCoverage.txt
    done
done

for age in young aged
do
    for f in A B
    do
        python3 get_summed_track_coverage.py --BW ${age}_bigwig_list.txt \
        --BED $prefix/04_FANC/compartmentExpression/compartmentBed/$age.$f.bed \
        --outf $outdir/${age}.$f.sumCoverage.txt
    done
done
