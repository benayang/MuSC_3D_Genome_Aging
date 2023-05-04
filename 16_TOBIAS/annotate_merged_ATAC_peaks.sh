uropa --bed young_aged_merged_ATAC_peaks.bed --gtf "/nas/homes/benyang/Genome_References/gencode.vM25.basic.annotation.gtf" --show_attributes gene_id gene_name --filter_attribute gene_type --attribute_values protein_coding --feature_anchor start --distance 10000 1000 --feature gene --threads 50

cut -f 1-6,16-17 young_aged_merged_ATAC_peaks_finalhits.txt | head -n 1 > young_aged_merged_ATAC_peaks_annotated_header.txt
cut -f 1-6,16-17 young_aged_merged_ATAC_peaks_finalhits.txt | tail -n +2 > young_aged_merged_ATAC_peaks_annotated.bed

Rscript --vanilla remove_version_number_gene_id.R young_aged_merged_ATAC_peaks_annotated.bed young_aged_merged_ATAC_peaks_annotated_parsed.bed

uropa --bed young_aged_merged_ATAC_peaks_merged_loop_anchors.bed --gtf "/nas/homes/benyang/Genome_References/gencode.vM25.basic.annotation.gtf" --show_attributes gene_id gene_name --filter_attribute gene_type --attribute_values protein_coding --feature_anchor start --distance 10000 1000 --feature gene --threads 50

cut -f 1-6,16-17 young_aged_merged_ATAC_peaks_merged_loop_anchors_finalhits.txt | head -n 1 > young_aged_merged_ATAC_peaks_merged_loop_anchors_annotated_header.txt
cut -f 1-6,16-17 young_aged_merged_ATAC_peaks_merged_loop_anchors_finalhits.txt | tail -n +2 > young_aged_merged_ATAC_peaks_merged_loop_anchors_annotated.bed