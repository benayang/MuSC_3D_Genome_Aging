import argparse
import pyBigWig

def parse_args(): 
    parser=argparse.ArgumentParser(description="get sum of BigWig values over set of BED regions")
    parser.add_argument("--BW")
    parser.add_argument("--BW_labels")
    parser.add_argument("--BED")
    parser.add_argument("--outf")
 
    return parser.parse_args()

def main(): 
    args=parse_args()
    
    regions=open(args.BED,'r').read().strip().split('\n')

    bw_list=open(args.BW,'r').read().strip().split('\n')
    bw_labels=[r.split('\t')[0] for r in bw_list]
    bw_files=[r.split('\t')[1] for r in bw_list]
    #bw_labels=[l.strip().split(' ') for l in args.BW_labels]
    assert len(bw_files) == len(bw_labels)
    num_bw=len(bw_files)

    bw_dict={}
    for i in range(num_bw):
        bw_dict[bw_labels[i]] = pyBigWig.open(bw_files[i],'r')

    outf=open(args.outf,'w')
    outf.write('chr\tstart\tend\t' + '\t'.join(bw_labels) + '\n')

    for line in regions:
        tokens=line.split('\t')
        chrom=tokens[0]
        start=tokens[1]
        end=tokens[2]

        bw_sums=[]
        for i in range(num_bw):
            bw_sums.append(sum(bw_dict[bw_labels[i]].values(chrom,int(start),int(end))))

        outf_line = '\t'.join([chrom, start, end]) + '\t' + '\t'.join([str(s) for s in bw_sums]) + '\n'
        outf.write(outf_line)
    outf.close()

if __name__=="__main__": 
    main() 