import argparse
import math 
def parse_args():
    parser=argparse.ArgumentParser(description="aggregate significant p-values")
    parser.add_argument("--inputs",nargs="+")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    data_dict=dict()
    thresh=0.05
    for inputf in args.inputs:
        print(str(inputf))
        data=open(inputf,'r').read().strip().split('\n')
        task=inputf.split('/')[1]
        data_dict[task]=dict() 
        for line in data[1::]:
            tokens=line.split('\t')
            motif=tokens[0].split('(')[0]
            pval=float(tokens[2])
            if pval>thresh:
                continue 
            pval=-10*math.log10(pval)
            data_dict[task][motif]=str(pval)
    outf=open(args.outf,'w')
    for task in data_dict:
        for motif in data_dict[task]:
            outf.write(task+'\t'+motif+'\t'+data_dict[task][motif]+'\n')
if __name__== "__main__":
    main()