prefix='/nas/homes/benyang/HiC/13_MultiOme/RawData/'

Rscript --vanilla plot_barcode_rank.R \
--indir $prefix/Aged/10x_analysis_5422-JL/Sample_5422-JL-1 \
--ATAC "/nas/homes/benyang/HiC/13_MultiOme/RawData/aged_raw_seurat.RDS" \
--outdir $prefix/barcode_rank_plots \
--maxlines 5e7 \
--prefix aged \
--threads 45

Rscript --vanilla plot_barcode_rank.R \
--indir $prefix/Young/10x_analysis_5644-JL/Sample_5644-JL-1 \
--ATAC "/nas/homes/benyang/HiC/13_MultiOme/RawData/young_raw_seurat.RDS" \
--outdir $prefix/barcode_rank_plots \
--maxlines NULL \
--prefix young \
--threads 45

Rscript --vanilla plot_barcode_rank.R \
--indir /nas/datasets/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1 \
--ATAC "/nas/homes/benyang/HiC/13_MultiOme/RawData/young_v2_raw_seurat.RDS" \
--outdir $prefix/barcode_rank_plots \
--maxlines 5e7 \
--prefix youngv2 \
--threads 45