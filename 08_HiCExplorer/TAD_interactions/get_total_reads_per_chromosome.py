import hicstraw
import pandas as pd
import numpy as np

aged_list = []
young_list = []

for chr in ['chr'+str(i+1) for i in range(19)] + ['chrM','chrX','chrY']:
    print(chr)

    aged_result = hicstraw.straw("observed","KR","/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.40kb.hic",chr,chr,"BP",40000)
    young_result = hicstraw.straw("observed","KR","/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic",chr,chr,"BP",40000)

    aged_count_sum = np.nansum([r.counts for r in aged_result])
    young_count_sum = np.nansum([r.counts for r in young_result])

    aged_list.append([chr, aged_count_sum])
    young_list.append([chr, young_count_sum])

aged_df = pd.DataFrame(aged_list)
aged_df.columns = ['chr','count_sum']
young_df = pd.DataFrame(young_list)
young_df.columns = ['chr','count_sum']

aged_df.to_csv('aged_total_counts_per_chromosome_40kb.txt',sep="\t",index=False)
young_df.to_csv('young_total_counts_per_chromosome_40kb.txt',sep="\t",index=False)