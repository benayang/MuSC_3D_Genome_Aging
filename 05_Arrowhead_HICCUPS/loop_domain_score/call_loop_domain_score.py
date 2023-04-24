import hicstraw
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description = "Call loop domain score.")
    parser.add_argument('--input_loops')
    parser.add_argument('--input_hic')
    parser.add_argument('--counts')
    parser.add_argument('--outfile')

    return parser.parse_args()

def main():
    args = parse_args()

    hic = hicstraw.HiCFile(args.input_hic)

    loopdomains = pd.read_csv(args.input_loops, sep='\t', header=None)

    chrom_name=[]
    chrom_length=[]
    for chrom in hic.getChromosomes():
        chrom_name.append(chrom.name)
        chrom_length.append(chrom.length)
    chrom_sizes = pd.DataFrame({'chr':chrom_name, 'length':chrom_length})

    get_data_by_row(hic, loopdomains, chrom_sizes, 5000, args.counts, args.outfile)

def get_data_by_row(hic, regions, chrom_sizes, res, counts, outfile):
    data_list = []
    for i in tqdm(range(regions.shape[0])):
        coords = regions.iloc[i,0:3]
        mzd = hic.getMatrixZoomData(coords[0],coords[0],counts,'KR','BP',res)
        data_list.append(get_region_data(mzd, coords, chrom_sizes))
    data = pd.DataFrame(data_list)
    data = pd.concat([regions.iloc[range(regions.shape[0])], data], axis=1)
    data.columns = ['chr','start','end','intraRegion_sum','intraRegion_count_nnz','interRegion_sum','interRegion_count_nnz']
    data.to_csv(outfile,sep=",",index=False)

def get_region_data(mzd, coords, chrom_sizes):
    intraRegion_data = get_intraRegion_records(mzd, coords, chrom_sizes)
    interRegion_data = get_interRegion_records(mzd, coords, chrom_sizes)
    return [np.nansum(intraRegion_data.counts), np.count_nonzero(intraRegion_data.counts), np.nansum(interRegion_data.counts), np.count_nonzero(interRegion_data.counts)]

def get_intraRegion_records(mzd,coords,chrom_sizes):
    binX=[]
    binY=[]
    counts=[]
    intraRegion_records=mzd.getRecords(coords[1],coords[2],coords[1],coords[2])
    for r in intraRegion_records:
        binX.append(int(r.binX))
        binY.append(int(r.binY))
        counts.append(r.counts)
    return pd.DataFrame({'X':binX, 'Y':binY, 'counts':counts})

def get_interRegion_records(mzd,coords,chrom_sizes):
    binX=[]
    binY=[]
    counts=[]
    
    left_region_records=mzd.getRecords(coords[1],coords[2],0,coords[1])
    for r in left_region_records:
        binX.append(int(r.binX))
        binY.append(int(r.binY))
        counts.append(r.counts)
    
    right_region_records=mzd.getRecords(coords[2],int(chrom_sizes.loc[chrom_sizes['chr']==coords[0],'length']),coords[1],coords[2])
    for r in right_region_records:
        binX.append(int(r.binX))
        binY.append(int(r.binY))
        counts.append(r.counts)
    return pd.DataFrame({'X':binX, 'Y':binY, 'counts':counts})

if __name__ == "__main__":
    main()