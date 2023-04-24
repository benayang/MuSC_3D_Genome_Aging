import pandas as pd

# young_A=pd.read_csv("/nas/homes/benyang/HiC/08_HiCExplorer/TAD expression/promoter_distance/young.promoter.A.boundaries.knownGenes.distance.bed",header=None,sep='\t')
# aged_A=pd.read_csv("/nas/homes/benyang/HiC/08_HiCExplorer/TAD expression/promoter_distance/aged.promoter.A.boundaries.knownGenes.distance.bed",header=None,sep='\t')
# young_B=pd.read_csv("/nas/homes/benyang/HiC/08_HiCExplorer/TAD expression/promoter_distance/young.promoter.B.boundaries.knownGenes.distance.bed",header=None,sep='\t')
# aged_B=pd.read_csv("/nas/homes/benyang/HiC/08_HiCExplorer/TAD expression/promoter_distance/aged.promoter.B.boundaries.knownGenes.distance.bed",header=None,sep='\t')

for f in ["young.promoter.A","aged.promoter.A","young.promoter.B","aged.promoter.B"]:
    data=pd.read_csv("/nas/homes/benyang/HiC/08_HiCExplorer/TAD expression/promoter_distance/"+f+".boundaries.knownGenes.distance.bed",header=None,sep='\t')
    boundary=data[data[10]<=5000]
    nonBoundary=data[data[10]>5000]
    boundary.iloc[:,0:4].drop_duplicates().to_csv(f+".boundary.knownGenes.bed",sep='\t',header=False,index=False)
    nonBoundary.iloc[:,0:4].drop_duplicates().to_csv(f+".nonBoundary.knownGenes.bed",sep='\t',header=False,index=False)