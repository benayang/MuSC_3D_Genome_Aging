import pandas as pd
import numpy as np
import hicstraw
import pybedtools
import matplotlib.pyplot as plt
from tqdm import tqdm
from TAD_interactions_functions import *

def main():
    prefix='/nas/homes/benyang/HiC'

    aged_hic = hicstraw.HiCFile(prefix+'/02_HIC/aged.merged/aged.merged.hic')
    young_hic = hicstraw.HiCFile(prefix+'/02_HIC/young.merged/young.merged.hic')

    tss = pybedtools.BedTool(prefix+'/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.bed')
    tss_slop = tss.slop(b=1000,genome="mm10")
    tss_slop.saveas('tss_slop_df.bed')
    tss_slop_df = tss_slop.to_dataframe()

    #top_2000_promoters = pd.read_csv('top2000_expressed_genes_per_group.txt',sep="\t",header=0)
    top_2000_promoters = pd.read_csv('DEtop2000_expressed_genes_per_group.txt',sep="\t",header=0)

    young_H3K4me3 = prefix+'/HistoneTracks/Tom_Rando/mm10LiftOver_GSM1148110_yQ_K4m3_peaks_w_input_control_peaks.bed'
    aged_H3K4me3 = prefix+'/HistoneTracks/Tom_Rando/mm10LiftOver_GSM1148118_oQSC_K4m3_MACSwInputControl.bed_peaks.bed'
    young_ATAC = prefix+'/HistoneTracks/atac.d0.Young.optimal.narrowPeak.gz'
    aged_ATAC = prefix+'/HistoneTracks/atac.d0.Aged.optimal.narrowPeak.gz'

    young_TAD = pd.read_csv(prefix+'/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed', sep='\t', header=None)
    young_TAD = pybedtools.BedTool.from_dataframe(young_TAD).saveas()
    aged_TAD = pd.read_csv(prefix+'/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed', sep='\t', header=None)
    aged_TAD = pybedtools.BedTool.from_dataframe(aged_TAD).saveas()

    chrom_name=[]
    chrom_length=[]
    for chrom in aged_hic.getChromosomes():
        chrom_name.append(chrom.name)
        chrom_length.append(chrom.length)
    chrom_sizes = pd.DataFrame({'chr':chrom_name, 'length':chrom_length})

    age = 'young'
    stat_df, interaction_df = get_summary_stats(young_hic,aged_hic,tss_slop_df,chrom_sizes,
    top_2000_promoters.loc[(top_2000_promoters['Age']==age) & (top_2000_promoters['Group']=="H3K4me3+/ATAC+")]['gene'].to_list(),5000,'observed')
    
    stat_df.to_csv(age+'DETop2000_posATAC_posH3K4me3_1kbPromoter_interaction_obs_counts.csv',sep=",",index=False)
    interaction_df.to_csv(age+'DETop2000_bothHiC_posATACposH3K4me3_1kb_obs_expanded_counts.csv',sep=",",index=False)

    # Need to separate Aged and Young interactions with gene promoter b/c may have different numbers of interactions
    young_cols = ['chr','start','end','gene'] + interaction_df.filter(regex='Young').columns.to_list() + ['ID']
    aged_cols = ['chr','start','end','gene'] + interaction_df.filter(regex='Aged').columns.to_list() + ['ID']
    # Can now drop the NaN values generated when pasting young and aged interactions together
    young_interaction_df = interaction_df[young_cols].dropna()
    aged_interaction_df = interaction_df[aged_cols].dropna()
    # Convert all bin indices to integers
    cols=young_interaction_df.filter(regex='_bin').columns.to_list() + ['start','end']
    young_interaction_df[cols] = young_interaction_df[cols].astype('int')
    cols=aged_interaction_df.filter(regex='_bin').columns.to_list() + ['start','end']
    aged_interaction_df[cols] = aged_interaction_df[cols].astype('int')
    # Save intermediate data tables
    young_interaction_df.to_csv(age+'DETop2000_youngHiC_posATACposH3K4me3_1kb_obs_expanded_counts.csv',sep=",",index=False)
    aged_interaction_df.to_csv(age+'DETop2000_agedHiC_posATACposH3K4me3_1kb_obs_expanded_counts.csv',sep=",",index=False)
    # young_interaction_df = pd.read_csv(age+'DETop2000_youngHiC_posATACposH3K4me3_1kb_obs_expanded_counts.csv',sep=",")
    # aged_interaction_df = pd.read_csv(age+'DETop2000_agedHiC_posATACposH3K4me3_1kb_obs_expanded_counts.csv',sep=",")
    # Separately classify Young and Aged interactions
    print('Classify young interactions.')
    young_mark_interaction_df = classify_interactions(young_interaction_df,'tss_slop_df.bed',young_ATAC,young_H3K4me3,aged_ATAC,aged_H3K4me3,young_TAD,aged_TAD)
    young_mark_interaction_df.to_csv(age+'DETop2000_youngHiC_posATACposH3K4me3_TAD_1kb_obs_counts_with_marks.csv',sep=",",index=False)
    print('Classify aged interactions.')
    aged_mark_interaction_df = classify_interactions(aged_interaction_df,'tss_slop_df.bed',young_ATAC,young_H3K4me3,aged_ATAC,aged_H3K4me3,young_TAD,aged_TAD)
    aged_mark_interaction_df.to_csv(age+'DETop2000_agedHiC_posATACposH3K4me3_TAD_1kb_obs_counts_with_marks.csv',sep=",",index=False)

    # avg_tpm = pd.concat([tpm.iloc[:,0:3].mean(axis=1), tpm.iloc[:,3:9].mean(axis=1)],axis=1)
    # avg_tpm.columns = ['Aged_TPM','Young_TPM']
    # avg_tpm = avg_tpm.applymap(lambda x: x if x>0 else 0.001)
    # avg_tpm = avg_tpm.assign(FC = lambda x: x['Young_TPM']/x['Aged_TPM'], logFC = lambda x: np.log2(x['FC']))

    # DE_genes = pd.read_csv(prefix+'/04_FANC/compartmentExpression/DE_0.05padj_log2(1.5)LFC.csv',sep="\t")
    # DE_genes.index = [s.split("_")[1] for s in DE_genes.index]
    # DE_genes = DE_genes.loc[np.logical_and((DE_genes.index.str.contains("Rik")==False),(DE_genes.index.str.contains("Gm")==False))]

    # DE_stat_df = get_summary_stats(young_hic,aged_hic,tss_slop_df,DE_genes.index.to_list(),5000,'obs')
    # DE_stat_df.to_csv(prefix+'/08_HiCExplorer/TAD_interactions/DE_genes_3kbPromoter_interaction_obs_counts.csv',sep=",",index=False)

if __name__ == "__main__":
    main()