data=pd.read_csv("/nas/home/benyang/Genome_References/gencode.vM25.basic.annotation.gtf",header=None,sep='\t',skiprows=5)
gene_names=[]
for index,row in data.iterrows():
    cur_gene_name=None
    gene_info=str(row[8]).split(';')
    for entry in gene_info:
        entry=entry.strip()
        if entry.startswith('gene_name'):
            cur_gene_name=entry.split('=')[1]
            break
    gene_names.append(cur_gene_name)
data=data[data.iloc[:,2]=="gene"]
data['gene']=gene_names
data.to_csv("/nas/home/benyang/HiC/08_HiCExplorer/gencode.vM25.basic.annotation.genes.gtf", sep="\t", header=False, index=False)
