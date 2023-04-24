import pyBigWig

mainDir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C"
h3k27me3Dir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/JC_H3K27me3"

# get histone marks
young_h3k27me3_fc_bw=pyBigWig.open(h3k27me3Dir+'/rep.pooled_x__1_1_NoAbYoung__2_1_NoAbYoung__3_1_NoAbYoung__R1.fastq.srt.nodup.fc.signal.bigwig','r')
aged_h3k27me3_fc_bw=pyBigWig.open(h3k27me3Dir+'/rep.pooled_x__1_5_NoAbAged__2_5_NoAbAged__3_5_NoAbAged__R1.fastq.srt.nodup.fc.signal.bigwig','r')
young_h3k27me3_rpkm_bw=pyBigWig.open(h3k27me3Dir+'/H3K27me3Young.counts.rpkm.track.bw','r')
aged_h3k27me3_rpkm_bw=pyBigWig.open(h3k27me3Dir+'/H3K27me3Aged.counts.rpkm.track.bw','r')
young_GSE_h3k27me3_bw=pyBigWig.open(mainDir+'/11_HistoneMods/GSE129749_yqsc_h3k27me3_fold-change_signal.bw','r')
young_5p_h4k20me1=pyBigWig.open(mainDir+'/11_HistoneMods/5p.counts.Young.bw','r')
aged_5p_h4k20me1=pyBigWig.open(mainDir+'/11_HistoneMods/5p.counts.Old.bw','r')
young_h4k20me1=pyBigWig.open(mainDir+'/11_HistoneMods/h4k20me1.counts.rpkm.track.Young.bw','r')
aged_h4k20me1=pyBigWig.open(mainDir+'/11_HistoneMods/h4k20me1.counts.rpkm.track.Old.bw','r')

# get gene bodies
a_to_b=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/A_to_B.genebodies.uniq.bed','r').read().strip().split('\n')
b_to_a=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/B_to_A.genebodies.uniq.bed','r').read().strip().split('\n')
static=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/static.genebodies.uniq.bed','r').read().strip().split('\n')
young_a=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/young.A.genebodies.uniq.bed','r').read().strip().split('\n')
young_b=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/young.B.genebodies.uniq.bed','r').read().strip().split('\n')
aged_a=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/aged.A.genebodies.uniq.bed','r').read().strip().split('\n')
aged_b=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/aged.B.genebodies.uniq.bed','r').read().strip().split('\n')

# create output files
a_to_b_outf=open(mainDir+'/08_HiCExplorer/TAD histone overlap/coverage/'+'H4K20me1.GSE129749.H3K27me3.coveraged.A_to_B.tsv','w')
b_to_a_outf=open(mainDir+'/08_HiCExplorer/TAD histone overlap/coverage/'+'H4K20me1.GSE129749.H3K27me3.coveraged.B_to_A.tsv','w')
static_outf=open(mainDir+'/08_HiCExplorer/TAD histone overlap/coverage/'+'H4K20me1.GSE129749.H3K27me3.coveraged.static.tsv','w')
young_a_outf=open(mainDir+'/08_HiCExplorer/TAD histone overlap/coverage/'+'H4K20me1.GSE129749.H3K27me3.coveraged.young.A.tsv','w')
young_b_outf=open(mainDir+'/08_HiCExplorer/TAD histone overlap/coverage/'+'H4K20me1.GSE129749.H3K27me3.coveraged.young.B.tsv','w')
aged_a_outf=open(mainDir+'/08_HiCExplorer/TAD histone overlap/coverage/'+'H4K20me1.GSE129749.H3K27me3.coveraged.aged.A.tsv','w')
aged_b_outf=open(mainDir+'/08_HiCExplorer/TAD histone overlap/coverage/'+'H4K20me1.GSE129749.H3K27me3.coveraged.aged.B.tsv','w')

for s in ['A_to_B','B_to_A','static','young.A','young.B','aged.A','aged.B']:
    outf=open(mainDir+'/08_HiCExplorer/TAD histone overlap/coverage/H4K20me1.GSE129749.H3K27me3.coveraged.{}.tsv'.format(s),'w')
    genebodies=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/{}.genebodies.uniq.bed'.format(s),'r').read().strip().split('\n')
    
    print(s)
    outf.write('chr\tstart\tend\tyoungGSEH3K27me3\tyoungFC_H3K27me3\tagedFC_H3K27me3\tyoungRPKM_H3K27me3\tagedRPKM_H3K27me3\tyoung5pH4K20me1\taged5pH4K20me1\tyoungH4K20me1\tagedH4K20me1\tGene\n')

    for line in genebodies:
        tokens=line.split('\t')
        chrom=tokens[0]
        start=tokens[1]
        end=tokens[2]
        gene=tokens[3]

        y_GSE_h3k27me3_sum=sum(young_GSE_h3k27me3_bw.values(chrom,int(start),int(end)))
        y_h3k27me3_fc_sum=sum(young_h3k27me3_fc_bw.values(chrom,int(start),int(end)))
        a_h3k27me3_fc_sum=sum(aged_h3k27me3_fc_bw.values(chrom,int(start),int(end)))
        y_h3k27me3_rpkm_sum=sum(young_h3k27me3_rpkm_bw.values(chrom,int(start),int(end)))
        a_h3k27me3_rpkm_sum=sum(aged_h3k27me3_rpkm_bw.values(chrom,int(start),int(end)))
        y_h4k20me1_5p_sum=sum(young_5p_h4k20me1.values(chrom,int(start),int(end)))
        a_h4k20me1_5p_sum=sum(aged_5p_h4k20me1.values(chrom,int(start),int(end)))
        y_h4k20me1_sum=sum(young_h4k20me1.values(chrom,int(start),int(end)))
        a_h4k20me1_sum=sum(aged_h4k20me1.values(chrom,int(start),int(end)))

        outf_line = '\t'.join([chrom, start, end,
        str(y_GSE_h3k27me3_sum),
        str(y_h3k27me3_fc_sum), str(a_h3k27me3_fc_sum), 
        str(y_h3k27me3_rpkm_sum), str(a_h3k27me3_rpkm_sum), 
        str(y_h4k20me1_5p_sum), str(a_h4k20me1_5p_sum), 
        str(y_h4k20me1_sum), str(a_h4k20me1_sum), gene]) + '\n'
        outf.write(outf_line)
    outf.close()


# print('a to b')
# a_to_b_outf.write('chr\tstart\tend\tH3K27me3\tyoungH4K20me1\tagedH4K20me1\tGene\n')
# for line in a_to_b:
#     tokens=line.split('\t')
#     chrom=tokens[0]
#     start=tokens[1]
#     end=tokens[2]
#     gene=tokens[3]
#     k27me3_sum=sum(young_h3k27me3_bw.values(chrom,int(start),int(end)))
#     y_k20me1_sum=sum(young_5p_h4k20me1.values(chrom,int(start),int(end)))
#     a_k20me1_sum=sum(aged_5p_h4k20me1.values(chrom,int(start),int(end)))
#     a_to_b_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(k27me3_sum)+'\t'+str(y_k20me1_sum)+'\t'+str(a_k20me1_sum)+'\t'+gene+'\n')
# a_to_b_outf.close()

# print('b to a')
# b_to_a_outf.write('chr\tstart\tend\tH3K27me3\tyoungH4K20me1\tagedH4K20me1\tGene\n')
# for line in b_to_a:
#     tokens=line.split('\t')
#     chrom=tokens[0]
#     start=tokens[1]
#     end=tokens[2]
#     gene=tokens[3]
#     k27me3_sum=sum(young_h3k27me3_bw.values(chrom,int(start),int(end)))
#     y_k20me1_sum=sum(young_5p_h4k20me1.values(chrom,int(start),int(end)))
#     a_k20me1_sum=sum(aged_5p_h4k20me1.values(chrom,int(start),int(end)))
#     b_to_a_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(k27me3_sum)+'\t'+str(y_k20me1_sum)+'\t'+str(a_k20me1_sum)+'\t'+gene+'\n')
# b_to_a_outf.close()

# print('static')
# static_outf.write('chr\tstart\tend\tH3K27me3\tyoungH4K20me1\tagedH4K20me1\tGene\n')
# for line in static:
#     tokens=line.split('\t')
#     chrom=tokens[0]
#     start=tokens[1]
#     end=tokens[2]
#     gene=tokens[3]
#     k27me3_sum=sum(young_h3k27me3_bw.values(chrom,int(start),int(end)))
#     y_k20me1_sum=sum(young_5p_h4k20me1.values(chrom,int(start),int(end)))
#     a_k20me1_sum=sum(aged_5p_h4k20me1.values(chrom,int(start),int(end)))
#     static_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(k27me3_sum)+'\t'+str(y_k20me1_sum)+'\t'+str(a_k20me1_sum)+'\t'+gene+'\n')
# static_outf.close()

# print('young A')
# young_a_outf.write('chr\tstart\tend\tH3K27me3\tyoungH4K20me1\tagedH4K20me1\tGene\n')
# for line in young_a:
#     tokens=line.split('\t')
#     chrom=tokens[0]
#     start=tokens[1]
#     end=tokens[2]
#     gene=tokens[3]
#     k27me3_sum=sum(young_h3k27me3_bw.values(chrom,int(start),int(end)))
#     y_k20me1_sum=sum(young_5p_h4k20me1.values(chrom,int(start),int(end)))
#     a_k20me1_sum=sum(aged_5p_h4k20me1.values(chrom,int(start),int(end)))
#     young_a_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(k27me3_sum)+'\t'+str(y_k20me1_sum)+'\t'+str(a_k20me1_sum)+'\t'+gene+'\n')
# young_a_outf.close()

# print('young B')
# young_b_outf.write('chr\tstart\tend\tH3K27me3\tyoungH4K20me1\tagedH4K20me1\tGene\n')
# for line in young_b:
#     tokens=line.split('\t')
#     chrom=tokens[0]
#     start=tokens[1]
#     end=tokens[2]
#     gene=tokens[3]
#     k27me3_sum=sum(young_h3k27me3_bw.values(chrom,int(start),int(end)))
#     y_k20me1_sum=sum(young_5p_h4k20me1.values(chrom,int(start),int(end)))
#     a_k20me1_sum=sum(aged_5p_h4k20me1.values(chrom,int(start),int(end)))
#     young_b_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(k27me3_sum)+'\t'+str(y_k20me1_sum)+'\t'+str(a_k20me1_sum)+'\t'+gene+'\n')
# young_b_outf.close()

# print('aged A')
# aged_a_outf.write('chr\tstart\tend\tH3K27me3\tyoungH4K20me1\tagedH4K20me1\tGene\n')
# for line in aged_a:
#     tokens=line.split('\t')
#     chrom=tokens[0]
#     start=tokens[1]
#     end=tokens[2]
#     gene=tokens[3]
#     k27me3_sum=sum(young_h3k27me3_bw.values(chrom,int(start),int(end)))
#     y_k20me1_sum=sum(young_5p_h4k20me1.values(chrom,int(start),int(end)))
#     a_k20me1_sum=sum(aged_5p_h4k20me1.values(chrom,int(start),int(end)))
#     aged_a_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(k27me3_sum)+'\t'+str(y_k20me1_sum)+'\t'+str(a_k20me1_sum)+'\t'+gene+'\n')
# aged_a_outf.close()

# print('aged B')
# aged_b_outf.write('chr\tstart\tend\tH3K27me3\tyoungH4K20me1\tagedH4K20me1\tGene\n')
# for line in aged_b:
#     tokens=line.split('\t')
#     chrom=tokens[0]
#     start=tokens[1]
#     end=tokens[2]
#     gene=tokens[3]
#     k27me3_sum=sum(young_h3k27me3_bw.values(chrom,int(start),int(end)))
#     y_k20me1_sum=sum(young_5p_h4k20me1.values(chrom,int(start),int(end)))
#     a_k20me1_sum=sum(aged_5p_h4k20me1.values(chrom,int(start),int(end)))
#     aged_b_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(k27me3_sum)+'\t'+str(y_k20me1_sum)+'\t'+str(a_k20me1_sum)+'\t'+gene+'\n')
# aged_b_outf.close()