import pandas as pd
import numpy as np
import hicstraw
import pyBigWig
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt 

def main():
    prefix='/nas/homes/benyang/HiC'

    aged_hic = hicstraw.HiCFile(prefix+'/02_HIC/aged.merged/aged.merged.hic')
    young_hic = hicstraw.HiCFile(prefix+'/02_HIC/young.merged/young.merged.hic')

    chrom_name=[]
    chrom_length=[]
    for chrom in aged_hic.getChromosomes():
        chrom_name.append(chrom.name)
        chrom_length.append(chrom.length)
    chrom_sizes = pd.DataFrame({'chr':chrom_name, 'length':chrom_length})

    ### aged bigwigs
    bw1 = pyBigWig.open(prefix+"/HistoneTracks/get_signal_bigwigs/macs2_signal/Aged_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig")
    bw2 = pyBigWig.open("/nas/homes/annashch/Age_ATAC/caper_out/atac/b4245fba-248f-4201-bd66-907bb8e97a85/call-macs2_signal_track_pooled/execution/MuSC_d0_Old_10K_R2_r1.merged.nodup.tn5.pooled.fc.signal.bigwig")
    bw3 = pyBigWig.open(prefix+"/HistoneTracks/get_signal_bigwigs/macs2_signal/Aged_H3K27me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig")
    bw4 = pyBigWig.open(prefix+"/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bw")
    bw5 = pyBigWig.open(prefix+"/04_FANC/without_KR_normalization/aged.merged_100kb_noKR.bw")
    ### young bigwigs
    # bw1 = pyBigWig.open(prefix+"/HistoneTracks/get_signal_bigwigs/macs2_signal/Young_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig")
    # bw2 = pyBigWig.open("/nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-macs2_signal_track_pooled/execution/MuSC_d0_Y_10K_R1_r1.merged.nodup.tn5.pooled.fc.signal.bigwig")
    # bw3 = pyBigWig.open(prefix+"/HistoneTracks/get_signal_bigwigs/macs2_signal/Young_H3K27me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig")
    # bw4 = pyBigWig.open(prefix+"/04_FANC/without_KR_normalization/young.merged_100kb_noKR.bw")
    # bw5 = pyBigWig.open(prefix+"/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bw")

    for chrom in (['chr'+str(i+1) for i in range(19)] + ['chrX']):
        chr = str(chrom)
        matrix_oe = aged_hic.getMatrixZoomData(chr, chr, "oe", "KR", "BP", 250000)
        start = 3000000
        end = int(chrom_sizes.loc[chrom_sizes['chr']==chr,'length'])
        numpy_matrix = matrix_oe.getRecordsAsMatrix(start, end, start, end)

        signal1 = np.array(bw1.stats(chr, start, end, type="mean", nBins=numpy_matrix.shape[0])).astype(float)
        signal2 = np.array(bw2.stats(chr, start, end, type="mean", nBins=numpy_matrix.shape[0])).astype(float)
        signal3 = np.array(bw3.stats(chr, start, end, type="mean", nBins=numpy_matrix.shape[0])).astype(float)
        signal4 = np.array(bw4.stats(chr, start, end, type="mean", nBins=numpy_matrix.shape[0])).astype(float)
        signal5 = np.array(bw5.stats(chr, start, end, type="mean", nBins=numpy_matrix.shape[0])).astype(float)

        signal_list = [signal1, signal2, signal3, signal4, signal5]
        signal_name_list = ['H3K4me3','ATAC','H3K27me3','TAD_Score','PC1']

        plot_2dhic_1dtrack_map(np.corrcoef(numpy_matrix), signal_list, signal_name_list, numpy_matrix.shape[0], 
        'hic_1d_2d_plots/aged/aged_'+chr+'_2d_1d_mat.png')

        for n, s in enumerate(signal_list):
            indices, numpy_matrix2 = get_sorted_indices(numpy_matrix, s)
            signal_sorted_list = [signal1[indices], signal2[indices], signal3[indices], signal4[indices], signal5[indices]]
            plot_2dhic_1dtrack_map(np.corrcoef(numpy_matrix2), signal_sorted_list, signal_name_list, numpy_matrix2.shape[0], 
            'hic_1d_2d_plots/aged/aged_'+chr+'_2d_1d_'+signal_name_list[n]+'_sorted_mat.png')

def get_sorted_indices(matrix, signal):
    indices = np.argsort(signal)
    numpy_matrix2 = matrix.copy()
    numpy_matrix2 = numpy_matrix2[:,indices]
    numpy_matrix2 = numpy_matrix2[indices,:]
    return indices,numpy_matrix2

def plot_2dhic_1dtrack_map(dense_matrix, signal, signal_names, nBins, fname):
    #d2 = np.log(dense_matrix)
    d2 = dense_matrix
    d2[np.isnan(d2)] = 0
    d2[np.isinf(d2)] = 0
    fig = plt.figure()
    nrow = 1+len(signal)
    #fig, axs = plt.subplots(nrow,1)

    fig.set_figheight(12+1.25*len(signal))
    fig.set_figwidth(8)
    spec = gridspec.GridSpec(ncols=1, nrows=nrow+1,
    width_ratios=[1], wspace=1,
    hspace=.15, height_ratios=[6]*len(signal) + [30,1])
    
    cols = ['red','blue','green','orange','brown','pink','yellow','black','grey']
    
    for i in range(len(signal)-1):
        ax = fig.add_subplot(spec[i])
        ax.plot(np.arange(len(signal[i])), signal[i], color=cols[i])
        plt.xlim(0,nBins)
        plt.ylabel(signal_names[i])
    # need separate code to plot PC1 data so we can add a dashed line in background
    ax = fig.add_subplot(spec[i+1])
    plt.axhline(y=0, color='black', linestyle='--')
    ax.plot(np.arange(len(signal[i+1])), signal[i+1], color=cols[i+1])
    plt.xlim(0,nBins)
    plt.ylabel(signal_names[i+1])

    ax = fig.add_subplot(spec[nrow-1])
    mat = ax.matshow(d2, cmap='bwr', vmin=-1, vmax=1, aspect='equal')
    ax.xaxis.set_ticks_position('bottom')
    plt.xlim(0,nBins)
    plt.ylim(nBins,0)
    
    colorAx = fig.add_subplot(spec[nrow])
    cb = fig.colorbar(mat, cax = colorAx, orientation="horizontal")
    cb.set_label('Pearson Correlation')
    #fig.colorbar(mat,orientation="horizontal")
    plt.savefig(fname,dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()