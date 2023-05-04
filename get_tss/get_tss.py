
import pandas as pd
#gtf file source: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz
data=pd.read_csv("/media/labuser/BenPassport/GenomeRefs/gencode.vM25.basic.annotation.gtf",header=None,sep='\t',skiprows=7)
gene_names=[] 
for index,row in data.iterrows():
    cur_gene_name=None
    gene_info=str(row[8]).split(';')
    for entry in gene_info:
        entry=entry.strip()
        if entry.startswith('gene_name'):
            cur_gene_name=entry.split(' ')[1].strip('"')
            break
    gene_names.append(cur_gene_name)
data['gene']=gene_names 
data=data[data.iloc[:,2]=="gene"]
data.to_csv("data.csv", sep="\t", header=False, index=False)
plus_strand=data.loc[data.iloc[:,6]=="+"]
minus_strand=data.loc[data.iloc[:,6]=='-']
plus_strand['chrom']=plus_strand[0]
plus_strand['start']=plus_strand[3].astype(int)
# plus_strand[3].apply(float.is_integer).all() # test if floats are integer values
plus_strand['end']=(plus_strand[3]+1).astype(int)
minus_strand['chrom']=minus_strand[0]
minus_strand['end']=minus_strand[4].astype(int)
minus_strand['start']=(minus_strand[4]-1).astype(int)
#data=data[data[2]=="gene"]
#plus_strand=data[data[6]=="+"]
#minus_strand=data[data[6]=='-']
#plus_strand['chrom']=plus_strand[0]
#plus_strand['start']=plus_strand[3]
#plus_strand['end']=plus_strand[3]+1
#minus_strand['chrom']=minus_strand[0]
#minus_strand['end']=minus_strand[4]
#minus_strand['start']=minus_strand[4]-1
plus_subset=plus_strand[['chrom','start','end','gene',5,6]]
minus_subset=minus_strand[['chrom','start','end','gene',5,6]]
merged=pd.concat((plus_subset,minus_subset),axis=0)
merged.to_csv("tss.gencode.vM25.basic.annotation.bed",sep='\t',header=False,index=False)

# gene_names = [g.strip('"') for g in gene_names]
# data_subset = data.iloc[:,[0,3,4,9,5,6,7]].drop_duplicates()
# data_subset.loc[~data_subset.duplicated(subset=['gene'])].to_csv("genebodies.genecode.vM25.basic.annotation.bed", sep='\t', header=False, index=False)


