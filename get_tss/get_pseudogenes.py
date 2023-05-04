import pandas as pd
#gtf file source: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz
data=pd.read_csv("/nas/homes/benyang/H4K20me1/get_tss/gencode.vM25.basic.annotation.gtf",header=None,sep='\t',skiprows=7)
gene_names=[] 
for index,row in data.iterrows():
    cur_gene_name=None
    gene_info=str(row[8]).split(';')
    for entry in gene_info:
        entry=entry.strip()
        if not entry.contains('pseudogene'):
            break
        if entry.startswith('gene_name'):
            cur_gene_name=entry.split('=')[1]
            break
    gene_names.append(cur_gene_name)
data.to_csv("pseudogene_data.csv", sep="\t", header=False, index=False)
# data['gene']=gene_names 
# data=data[data.iloc[:,2]=="gene"]
# plus_strand=data[data.iloc[:,6]=="+"]
# minus_strand=data[data.iloc[:,6]=='-']
# plus_strand.loc[:,'chrom']=plus_strand[0]
# plus_strand.loc[:,'start']=plus_strand[3].astype(int)
# # plus_strand[3].apply(float.is_integer).all() # test if floats are integer values
# plus_strand.loc[:,'end']=(plus_strand[3]+1).astype(int)
# minus_strand.loc[:,'chrom']=minus_strand[0]
# minus_strand.loc[:,'end']=minus_strand[4].astype(int)
# minus_strand.loc[:,'start']=(minus_strand[4]-1).astype(int)
#data=data[data[2]=="gene"]
#plus_strand=data[data[6]=="+"]
#minus_strand=data[data[6]=='-']
#plus_strand['chrom']=plus_strand[0]
#plus_strand['start']=plus_strand[3]
#plus_strand['end']=plus_strand[3]+1
#minus_strand['chrom']=minus_strand[0]
#minus_strand['end']=minus_strand[4]
#minus_strand['start']=minus_strand[4]-1
# plus_subset=plus_strand[['chrom','start','end','gene']]
# minus_subset=minus_strand[['chrom','start','end','gene']]
# merged=pd.concat((plus_subset,minus_subset),axis=0)
# merged.to_csv("tss.gencode.vM25.basic.annotation.bed",sep='\t',header=False,index=False)




