import pandas as pd
import numpy as np
import pybedtools
from tqdm import tqdm

def get_binned_gene(gene,res,tss_slop_df):
  if str(gene) in tss_slop_df['name'].to_list():
    row = np.array(tss_slop_df.loc[tss_slop_df['name']==str(gene)].iloc[0])
    start = int(row[1]/res)*res
    end = int(row[2]/res)*res
    return [row[0],start,end]
  else:
    return "empty"

def get_binned_region(row,res):
  start = int(row[1]/res)*res
  end = int(row[2]/res)*res
  return [row[0],start,end]

def get_interaction_records(mzd,coords,chrom_sizes):
  binX=[]
  binY=[]
  counts=[]
  records=mzd.getRecords(coords[1],coords[2],0,int(chrom_sizes.loc[chrom_sizes['chr']==coords[0],'length']))
  for r in records:
    if (coords[1] <= r.binX <= coords[2]) or ((coords[1] <= r.binY <= coords[2])):
    #if r.binX==coords[1] or r.binX==coords[2] or r.binY==coords[1] or r.binY==coords[2]:
      binX.append(int(r.binX))
      binY.append(int(r.binY))
      counts.append(r.counts)
  return pd.DataFrame({'X':binX, 'Y':binY, 'counts':counts})

def bin_class(start,end,x):
  'Promoter' if (x>=start and x<=end) else 'Non-promoter'

def classify_interactions(interaction_df,tss_slop,young_ATAC,young_H3K4me3,aged_ATAC,aged_H3K4me3,young_TAD,aged_TAD):
  # need to define tss_slop as a string to file for pybedtools annotate
  #tss_slop = '/content/drive/MyDrive/HiC/TAD_interactions/tss_slop_df.bed'

  # convert upstream and downstream bins to BED files
  # young_interaction_binX_bed = pybedtools.BedTool.from_dataframe(interaction_df[['chr','Young_binX','Young_binX_end','ID']]).saveas()
  # young_interaction_binY_bed = pybedtools.BedTool.from_dataframe(interaction_df[['chr','Young_binY','Young_binY_end','ID']]).saveas()
  # aged_interaction_binX_bed = pybedtools.BedTool.from_dataframe(interaction_df[['chr','Aged_binX','Aged_binX_end','ID']]).saveas()
  # aged_interaction_binY_bed = pybedtools.BedTool.from_dataframe(interaction_df[['chr','Aged_binY','Aged_binY_end','ID']]).saveas()

  binX_cols = interaction_df.filter(regex='binX').columns.to_list()
  binY_cols = interaction_df.filter(regex='binY').columns.to_list()
  binX_bed = pybedtools.BedTool.from_dataframe(interaction_df[['chr']+binX_cols+['ID']]).saveas()
  binY_bed = pybedtools.BedTool.from_dataframe(interaction_df[['chr']+binY_cols+['ID']]).saveas()

  # Annotate the BED files with ATAC and H3K4me3 peaksets then convert to Pandas dataframes
  # young_interaction_binX_bed = binX_bed.annotate(files=[young_H3K4me3,young_ATAC,tss_slop],counts=True).to_dataframe(names=['chr','Young_binX','Young_binX_end','ID','Young_binX_H3K4me3','Young_binX_ATAC','Young_binX_tss_1kb'])
  # young_interaction_binY_bed = binY_bed.annotate(files=[young_H3K4me3,young_ATAC,tss_slop],counts=True).to_dataframe(names=['chr','Young_binY','Young_binY_end','ID','Young_binY_H3K4me3','Young_binY_ATAC','Young_binY_tss_1kb'])
  # aged_interaction_binX_bed = binX_bed.annotate(files=[aged_H3K4me3,aged_ATAC,tss_slop],counts=True).to_dataframe(names=['chr','Aged_binX','Aged_binX_end','ID','Aged_binX_H3K4me3','Aged_binX_ATAC','Aged_binX_tss_1kb'])
  # aged_interaction_binY_bed = binY_bed.annotate(files=[aged_H3K4me3,aged_ATAC,tss_slop],counts=True).to_dataframe(names=['chr','Aged_binY','Aged_binY_end','ID','Aged_binY_H3K4me3','Aged_binY_ATAC','Aged_binY_tss_1kb'])

  young_interaction_binX_bed = binX_bed.annotate(files=[young_H3K4me3,young_ATAC,tss_slop],counts=True)
  young_interaction_binY_bed = binY_bed.annotate(files=[young_H3K4me3,young_ATAC,tss_slop],counts=True)
  aged_interaction_binX_bed = binX_bed.annotate(files=[aged_H3K4me3,aged_ATAC,tss_slop],counts=True)
  aged_interaction_binY_bed = binY_bed.annotate(files=[aged_H3K4me3,aged_ATAC,tss_slop],counts=True)

  filter_cols = [0,1,2,3,4,5,6,8,9,10,11]
  young_interaction_binX_bed = young_interaction_binX_bed.intersect(young_TAD, wao=True).to_dataframe(disable_auto_names=True, header=None, low_memory=False).iloc[:,filter_cols]
  young_interaction_binX_bed.columns = ['chr','Young_binX','Young_binX_end','ID','Young_binX_H3K4me3','Young_binX_ATAC','Young_binX_tss_1kb','Young_binX_TAD_start','Young_binX_TAD_end','Young_binX_TAD_name','Young_binX_TAD_score']
  young_interaction_binY_bed = young_interaction_binY_bed.intersect(young_TAD, wao=True).to_dataframe(disable_auto_names=True, header=None, low_memory=False).iloc[:,filter_cols]
  young_interaction_binY_bed.columns = ['chr','Young_binY','Young_binY_end','ID','Young_binY_H3K4me3','Young_binY_ATAC','Young_binY_tss_1kb','Young_binY_TAD_start','Young_binY_TAD_end','Young_binY_TAD_name','Young_binY_TAD_score']
  aged_interaction_binX_bed = aged_interaction_binX_bed.intersect(aged_TAD, wao=True).to_dataframe(disable_auto_names=True, header=None, low_memory=False).iloc[:,filter_cols]
  aged_interaction_binX_bed.columns = ['chr','Aged_binX','Aged_binX_end','ID','Aged_binX_H3K4me3','Aged_binX_ATAC','Aged_binX_tss_1kb','Aged_binX_TAD_start','Aged_binX_TAD_end','Aged_binX_TAD_name','Aged_binX_TAD_score']
  aged_interaction_binY_bed = aged_interaction_binY_bed.intersect(aged_TAD, wao=True).to_dataframe(disable_auto_names=True, header=None, low_memory=False).iloc[:,filter_cols]
  aged_interaction_binY_bed.columns = ['chr','Aged_binY','Aged_binY_end','ID','Aged_binY_H3K4me3','Aged_binY_ATAC','Aged_binY_tss_1kb','Aged_binY_TAD_start','Aged_binY_TAD_end','Aged_binY_TAD_name','Aged_binY_TAD_score']
  # bind columns between young and aged datasets 
  bin_interaction = young_interaction_binX_bed.merge(young_interaction_binY_bed,on=['chr','ID']).merge(aged_interaction_binX_bed,on=['chr','ID']).merge(aged_interaction_binY_bed,on=['chr','ID'])

  # since pybedtools rearranges row order when creating from DataFrame, need to merge peakset data with original interaction DataFrame 
  output = interaction_df.merge(bin_interaction,on=['ID'],suffixes=['','_y'],how='left')

  # verify merge is correct
  # assert output['Young_binX'].equals(output['Young_binX_y'])
  # assert output['Young_binY'].equals(output['Young_binY_y'])
  # assert output['Young_binX_end'].equals(output['Young_binX_end_y'])
  # assert output['Young_binY_end'].equals(output['Young_binY_end_y'])
  # assert output['Aged_binX'].equals(output['Aged_binX_y'])
  # assert output['Aged_binY'].equals(output['Aged_binY_y'])
  # assert output['Aged_binX_end'].equals(output['Aged_binX_end_y'])
  # assert output['Aged_binY_end'].equals(output['Aged_binY_end_y'])

  # Drop the duplicate chr columns and rename as chr
  output.drop(output.filter(regex='_y$').columns.to_list(),axis=1, inplace=True)
  output.drop(output.filter(regex='binX_end$').columns.to_list(),axis=1, inplace=True)
  output.drop(output.filter(regex='binY_end$').columns.to_list(),axis=1, inplace=True)

  return output

def make_bedpe(data,chrom,res):
  data_start=data['X']
  data_end=data['Y']
  data_scores=data['counts']
  bedpe=pd.concat([pd.Series([chrom]*len(data)),data_start,data_start+res,
             pd.Series([chrom]*len(data)),data_end,data_end+res,
             data_scores],axis=1)
  bedpe.columns = ['chr1','start1','end1','chr2','start2','end2','counts']
  return bedpe

###### Find all interactions involving bins covering gene promoters
def get_summary_stats(young_hic,aged_hic,tss_slop_df,chrom_sizes,genes,res,datatype):
  coords_list=[]
  count_num=[]
  count_sum=[]
  count_mean=[]
  genes_list=[]
  interactions_list=[]

  for g in tqdm(genes):
    coords = get_binned_gene(g,res,tss_slop_df)
    if coords == "empty":
      continue
    young_mzd = young_hic.getMatrixZoomData(coords[0],coords[0],datatype,'KR','BP',res)
    aged_mzd = aged_hic.getMatrixZoomData(coords[0],coords[0],datatype,'KR','BP',res)
    young_data = get_interaction_records(young_mzd,coords,chrom_sizes)
    aged_data = get_interaction_records(aged_mzd,coords,chrom_sizes)

    coords_list.append(coords)
    genes_list.append(g)
    count_num.append([young_data.shape[0], aged_data.shape[0]]) # get non-zero counts
    count_sum.append([np.nansum(young_data['counts']), np.nansum(aged_data['counts'])]) # skip NaN
    count_mean.append([np.nanmean(young_data['counts']), np.nanmean(aged_data['counts'])])

    max_nrow = max(young_data.shape[0],aged_data.shape[0])
    interactions_list.append(pd.concat([pd.DataFrame([coords]*max_nrow),pd.Series([g]*max_nrow),young_data,aged_data],axis=1,ignore_index=True))
  
  coords_df = pd.DataFrame(coords_list)
  count_num_df = pd.DataFrame(count_num)
  count_sum_df = pd.DataFrame(count_sum)
  count_mean_df = pd.DataFrame(count_mean)
  output_stat = pd.concat([coords_df,pd.Series(genes_list),count_num_df,count_sum_df,count_mean_df], axis=1)
  output_stat.columns = ["chr","start","end","gene","Young_num","Aged_num","Young_sum","Aged_sum","Young_mean","Aged_mean"]
  
  interaction_df = pd.concat(interactions_list)
  interaction_df.columns = ['chr','start','end','gene','Young_binX','Young_binY','Young_counts','Aged_binX','Aged_binY','Aged_counts']
  interaction_df = interaction_df.assign(Young_binX_end=lambda x: x.Young_binX + res,
                                         Young_binY_end=lambda x: x.Young_binY + res,
                                         Aged_binX_end=lambda x: x.Aged_binX + res,
                                         Aged_binY_end=lambda x: x.Aged_binY + res,
                                         ID=[i for i in range(interaction_df.shape[0])]) 

  return output_stat, interaction_df

###### Add gene expression data to interaction data
def merge_genes_with_interactions(stat_sum_df, gene_list):
  stat_sum_df.index=stat_sum_df.genes.values

  stat_merged=stat_sum_df.join(gene_list,how='left',lsuffix='_left', rsuffix='_right')
  stat_merged['color'] = stat_merged['logFC'].apply(lambda x: 'blue' if x>0 else 'red')
  stat_merged['direction'] = stat_merged['logFC'].apply(lambda x: 'Young' if x>0 else 'Aged')
  
  return stat_merged