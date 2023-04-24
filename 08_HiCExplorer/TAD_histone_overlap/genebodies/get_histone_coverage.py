import pyBigWig

mainDir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C"

young_h3k27me3_bw=pyBigWig.open(mainDir+'/11_HistoneMods/GSE129749_yqsc_h3k27me3_fold-change_signal.bw','r')
a_to_b=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/A_to_B.genebodies.bed','r').read().strip().split('\n')
b_to_a=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/B_to_A.genebodies.bed','r').read().strip().split('\n')
static=open(mainDir+'/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies/static.genebodies.bed','r').read().strip().split('\n')

a_to_b_outf=open('young.GSE129749.H3K27me3.coveraged.a_to_b.txt','w')
b_to_a_outf=open('young.GSE129749.H3K27me3.coveraged.b_to_a.txt','w')
static_outf=open('young.GSE129749.H3K27me3.coveraged.static.txt','w')

a_to_b_outf.write('chr\tstart\tend\tH3K27me3Coverage\tGene\n')
for line in a_to_b:
    tokens=line.split('\t')
    chrom=tokens[0]
    start=tokens[1]
    end=tokens[2]
    gene=tokens[3]
    k27me3_sum=sum(young_h3k27me3_bw.values(chrom,int(start),int(end)))
    a_to_b_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(k27me3_sum)+'\t'+gene+'\n')
a_to_b_outf.close()

b_to_a_outf.write('chr\tstart\tend\tH3K27me3Coverage\tGene\n')
for line in b_to_a:
    tokens=line.split('\t')
    chrom=tokens[0]
    start=tokens[1]
    end=tokens[2]
    young_sum=sum(young_bw.values(chrom,int(start),int(end)))
    aged_sum=sum(aged_bw.values(chrom,int(start),int(end)))
    b_to_a_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(aged_sum)+'\t'+str(young_sum)+'\n')
b_to_a_outf.close()

static_outf.write('chr\tstart\tend\tH3K27me3Coverage\tGene\n')
for line in static:
    tokens=line.split('\t')
    chrom=tokens[0]
    start=tokens[1]
    end=tokens[2]
    young_sum=sum(young_bw.values(chrom,int(start),int(end)))
    aged_sum=sum(aged_bw.values(chrom,int(start),int(end)))
    static_outf.write(chrom+'\t'+start+'\t'+end+'\t'+str(aged_sum)+'\t'+str(young_sum)+'\n')
static_outf.close()