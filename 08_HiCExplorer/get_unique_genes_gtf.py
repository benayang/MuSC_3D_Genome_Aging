import pandas as pd
import csv

data=pd.read_csv("/nas/homes/benyang/Genome_References/gencode.vM25.basic.annotation.gtf",header=None,sep='\t',skiprows=5)
data=data[data.iloc[:,2]=="gene"]
data.to_csv("/nas/homes/benyang/HiC/08_HiCExplorer/gencode.vM25.basic.annotation.genes.bed", sep='\t', header=None, index=False)

# cat ./gtf_header.txt \
# ./gencode.vM25.basic.annotation.genes.bed > \
# ./gencode.vM25.basic.annotation.genes.gtf
