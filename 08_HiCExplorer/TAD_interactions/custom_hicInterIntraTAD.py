import hicstraw
import pandas as pd
import numpy as np
from tqdm import tqdm

def main():
    prefix='/nas/homes/benyang/HiC'

    aged_hic = hicstraw.HiCFile(prefix+'/02_HIC/aged.merged/aged.merged.40kb.hic')
    young_hic = hicstraw.HiCFile(prefix+'/02_HIC/young.merged/young.merged.hic')

    young_TAD = pd.read_csv(prefix+'/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed', sep='\t', header=None)
    aged_TAD = pd.read_csv(prefix+'/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed', sep='\t', header=None)

    chrom_name=[]
    chrom_length=[]
    for chrom in aged_hic.getChromosomes():
        chrom_name.append(chrom.name)
        chrom_length.append(chrom.length)
    chrom_sizes = pd.DataFrame({'chr':chrom_name, 'length':chrom_length})

    get_TAD_data_by_row(young_hic, young_TAD, chrom_sizes, 40000, prefix, 'young')
    get_TAD_data_by_row(aged_hic, aged_TAD, chrom_sizes, 40000, prefix, 'aged')

def get_TAD_data_by_row(hic, TAD, chrom_sizes, res, prefix, age):
    TAD_data_list = []
    for i in tqdm(range(TAD.shape[0])):
        coords = TAD.iloc[i,0:3]
        mzd = hic.getMatrixZoomData(coords[0],coords[0],'observed','KR','BP',res)
        TAD_data_list.append(get_TAD_data(mzd, coords, chrom_sizes))
    TAD_data = pd.DataFrame(TAD_data_list)
    TAD_data = pd.concat([TAD.iloc[range(TAD.shape[0]), 0:5], TAD_data], axis=1)
    TAD_data.columns = ['chr','start','end','ID','Score','intraTAD_sum','intraTAD_count_nnz','interTAD_sum','interTAD_count_nnz']
    TAD_data.to_csv(prefix+'/08_HiCExplorer/TAD_interactions/custom_hicInterIntraTAD/'+age+'_hicInterIntraTAD.csv',sep=",",index=False)

def get_TAD_data(mzd, coords, chrom_sizes):
    intraTAD_data = get_intraTAD_records(mzd, coords, chrom_sizes)
    interTAD_data = get_interTAD_records(mzd, coords, chrom_sizes)
    return [np.nansum(intraTAD_data.counts), np.count_nonzero(intraTAD_data.counts), np.nansum(interTAD_data.counts), np.count_nonzero(interTAD_data.counts)]

def get_intraTAD_records(mzd,coords,chrom_sizes):
    binX=[]
    binY=[]
    counts=[]
    intra_TAD_records=mzd.getRecords(coords[1],coords[2],coords[1],coords[2])
    for r in intra_TAD_records:
        binX.append(int(r.binX))
        binY.append(int(r.binY))
        counts.append(r.counts)
    return pd.DataFrame({'X':binX, 'Y':binY, 'counts':counts})

def get_interTAD_records(mzd,coords,chrom_sizes):
    binX=[]
    binY=[]
    counts=[]
    
    left_TAD_records=mzd.getRecords(coords[1],coords[2],0,coords[1])
    for r in left_TAD_records:
        binX.append(int(r.binX))
        binY.append(int(r.binY))
        counts.append(r.counts)
    
    right_TAD_records=mzd.getRecords(coords[2],int(chrom_sizes.loc[chrom_sizes['chr']==coords[0],'length']),coords[1],coords[2])
    for r in right_TAD_records:
        binX.append(int(r.binX))
        binY.append(int(r.binY))
        counts.append(r.counts)
    return pd.DataFrame({'X':binX, 'Y':binY, 'counts':counts})

if __name__ == "__main__":
    main()