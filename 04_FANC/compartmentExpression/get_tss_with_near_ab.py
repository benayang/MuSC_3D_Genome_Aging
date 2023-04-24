import pandas as pd

aged=pd.read_csv("./compartmentBed/100kb/tss/aged.ab.gene.bfilt.distance.bed",header=None,sep='\t')
aged=aged[aged[10]<500]
aged_coord=aged[[4,5,6,7,3]]
aged_coord=aged_coord.drop_duplicates() 
aged_coord.to_csv("./compartmentBed/100kb/tss/aged.tss.near.ab.bed",header=False,index=False,sep='\t')

young=pd.read_csv("./compartmentBed/100kb/tss/young.ab.gene.bfilt.distance.bed",header=None,sep='\t')
young=young[young[10]<500]
young_coord=young[[4,5,6,7,3]]
young_coord=young_coord.drop_duplicates() 
young_coord.to_csv("./compartmentBed/100kb/tss/young.tss.near.ab.bed",header=False,index=False,sep='\t')

a_to_b=pd.read_csv("./compartmentBed/100kb/tss/A_to_B.gene.bfilt.distance.bed",header=None,sep='\t')
a_to_b=a_to_b[a_to_b[9]<500]
a_to_b_coord=a_to_b[[3,4,5,6]]
a_to_b_coord=a_to_b_coord.drop_duplicates() 
a_to_b_coord.to_csv("./compartmentBed/100kb/tss/A_to_B.tss.near.ab.bed",header=False,index=False,sep='\t')

b_to_a=pd.read_csv("./compartmentBed/100kb/tss/B_to_A.gene.bfilt.distance.bed",header=None,sep='\t')
b_to_a=b_to_a[b_to_a[9]<500]
b_to_a_coord=b_to_a[[3,4,5,6]]
b_to_a_coord=b_to_a_coord.drop_duplicates() 
b_to_a_coord.to_csv("./compartmentBed/100kb/tss/B_to_A.tss.near.ab.bed",header=False,index=False,sep='\t')

static=pd.read_csv("./compartmentBed/100kb/tss/static.gene.bfilt.distance.bed",header=None,sep='\t')
static=static[static[9]<500]
static_coord=static[[3,4,5,6]]
static_coord=static_coord.drop_duplicates() 
static_coord.to_csv("./compartmentBed/100kb/tss/static.tss.near.ab.bed",header=False,index=False,sep='\t')