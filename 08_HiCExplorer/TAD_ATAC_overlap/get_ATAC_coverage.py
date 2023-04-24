import pyBigWig

mainDir="/nas/homes/benyang/HiC/"

young_ATAC_bw=pyBigWig.open('/nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-macs2_signal_track_pooled/execution/MuSC_d0_Y_10K_R1_r1.merged.nodup.tn5.pooled.fc.signal.bigwig','r')
aged_ATAC_bw=pyBigWig.open('/nas/homes/annashch/Age_ATAC/caper_out/atac/b4245fba-248f-4201-bd66-907bb8e97a85/call-macs2_signal_track_pooled/execution/MuSC_d0_Old_10K_R2_r1.merged.nodup.tn5.pooled.fc.signal.bigwig','r')

young_promoter_boundary=open(mainDir+'08_HiCExplorer/TAD expression/young.merged.TAD.promoter.boundaries.bed','r').read().strip().split('\n')
aged_promoter_boundary=open(mainDir+'08_HiCExplorer/TAD expression/aged.merged.TAD.promoter.boundaries.bed','r').read().strip().split('\n')
young_non_promoter_boundary=open(mainDir+'08_HiCExplorer/TAD expression/young.merged.TAD.nonpromoter.boundaries.bed','r').read().strip().split('\n')
aged_non_promoter_boundary=open(mainDir+'08_HiCExplorer/TAD expression/aged.merged.TAD.nonpromoter.boundaries.bed','r').read().strip().split('\n')

young_promoter_outf=open('young.boundary.promoter.ATAC.coveraged.txt','w')
aged_promoter_outf=open('aged.boundary.promoter.ATAC.coveraged.txt','w')
young_non_promoter_outf=open('young.boundary.nonpromoter.ATAC.coveraged.txt','w')
aged_non_promoter_outf=open('aged.boundary.nonpromoter.ATAC.coveraged.txt','w')

young_promoter_outf.write('chr\tstart\tend\tATAC_Coverage\tGene\n')
for line in young_promoter_boundary:
    tokens=line.split('\t')
    chrom=tokens[0]
    start=tokens[1]
    end=tokens[2]
    gene=tokens[3]
    ATAC_sum=sum(young_ATAC_bw.values(chrom,int(start),int(end)))
    young_promoter_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(ATAC_sum)+'\t'+gene+'\n')
young_promoter_outf.close()

aged_promoter_outf.write('chr\tstart\tend\tATAC_Coverage\tGene\n')
for line in aged_promoter_boundary:
    tokens=line.split('\t')
    chrom=tokens[0]
    start=tokens[1]
    end=tokens[2]
    gene=tokens[3]
    ATAC_sum=sum(aged_ATAC_bw.values(chrom,int(start),int(end)))
    aged_promoter_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(ATAC_sum)+'\t'+gene+'\n')
aged_promoter_outf.close()

young_non_promoter_outf.write('chr\tstart\tend\tATAC_Coverage\tGene\n')
for line in young_non_promoter_boundary:
    tokens=line.split('\t')
    chrom=tokens[0]
    start=tokens[1]
    end=tokens[2]
    gene=tokens[3]
    ATAC_sum=sum(young_ATAC_bw.values(chrom,int(start),int(end)))
    young_non_promoter_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(ATAC_sum)+'\t'+gene+'\n')
young_non_promoter_outf.close()

aged_non_promoter_outf.write('chr\tstart\tend\tATAC_Coverage\tGene\n')
for line in aged_non_promoter_boundary:
    tokens=line.split('\t')
    chrom=tokens[0]
    start=tokens[1]
    end=tokens[2]
    gene=tokens[3]
    ATAC_sum=sum(aged_ATAC_bw.values(chrom,int(start),int(end)))
    aged_non_promoter_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(ATAC_sum)+'\t'+gene+'\n')
aged_non_promoter_outf.close()