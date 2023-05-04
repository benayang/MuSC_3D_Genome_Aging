require(MASS)
require(mclust)
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(GenomicRanges)
require(GenomicInteractions)
require(doParallel)
require(parallel)
require(maxmatching)

#-----------------------------------------------------------------------------------------------------------------------------------
# ATAC QC plot functions
#-----------------------------------------------------------------------------------------------------------------------------------

######################################
# get density coloring for QC metrics
######################################
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#################################################
# Plot log10(nFrags) against TSSEnrichment for QC
#################################################
plot_qc <- function(df, nfrag_thresh=0, tss_thresh=0) {
    p <- ggPoint(
            x = log10(df[['nFrags']] + 1), 
            y = log10(df[['TSSEnrichment']] + 1), 
            colorDensity = TRUE,
            continuousSet = "sambaNight",
            xlabel = "Log10 Unique Fragments",
            ylabel = "Log10 TSS Enrichment",
            rastr = TRUE,
            dpi = 300) + 
        geom_hline(yintercept = log10(tss_thresh + 1), lty = "dashed") + 
        geom_vline(xintercept = log10(nfrag_thresh + 1), lty = "dashed") +
        guides(color = guide_colorbar(title="Density")) +
        theme(legend.position = "top", 
            legend.text = element_text(size=10), 
            legend.title = element_text(size=12),
            #legend.key.width = unit(1, "cm"),
            axis.text = element_text(size=10), 
            axis.title = element_text(size=12))
    return(p)
}

plot_qc_facet <- function(df, facet_var, log10nfrag_thresh, tss_thresh) {
    plt_list <- lapply(unique(df[[facet_var]]), function(x) plot_qc(df[df[[facet_var]] == x, ], log10nfrag_thresh, tss_thresh) + ggtitle(x))
    return(plt_list)
}

#######################
# define GMM clustering
#######################
gmm_clustering <- function(archr_obj, minFrags=0, outdir) {
    df <- getCellColData(archr_obj, select = c("Age", "Sample", "nFrags", "TSSEnrichment", "BlacklistRatio", "NucleosomeRatio"))
    for (i in unique(projAging1$Sample)){
        plt_df <- df[df$Sample==i,]
        proj.i <- df[df$Sample==i & df$nFrags>minFrags,]

        # GMM for fragments per cell
        depth.clust <- Mclust(log10(proj.i$nFrags + 1), G = 2)
        proj.i$depth.cluster <- depth.clust$classification
        proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty

        # GMM for TSS per cell
        TSS.clust <- Mclust(log10(proj.i$TSSEnrichment + 1), G = 2)
        proj.i$TSS.cluster <- TSS.clust$classification
        proj.i$TSS.cluster.uncertainty <- TSS.clust$uncertainty

        # Define data frames for thresholds and clusters for depth and TSS enrichment
        df.TSS <- data.frame(cellNames = rownames(proj.i),
                             TSS.cluster = proj.i$TSS.cluster,
                             TSS.cluster.uncertainty = proj.i$TSS.cluster.uncertainty,
                             TSSEnrichment = proj.i$TSSEnrichment,
                             Sample = i)
        saveRDS(df.TSS, file.path(outdir, paste0("df_TSS_",i,".rds")))

        df.depth <- data.frame(cellNames = rownames(proj.i),
                               depth.cluster = proj.i$depth.cluster,
                               depth.cluster.uncertainty = proj.i$depth.cluster.uncertainty,
                               nFrags = proj.i$nFrags,
                               Sample = i)
        saveRDS(df.depth, file.path(outdir, paste0("df_depth_",i,".rds")))

        # Filter cells by GMM thresholds for plotting
        df.TSS <- dplyr::filter(df.TSS, TSS.cluster == 2)
        df.TSS <- dplyr::filter(df.TSS, TSS.cluster.uncertainty <= 0.05)
        df.depth <- dplyr::filter(df.depth, depth.cluster == 2)
        df.depth <- dplyr::filter(df.depth, depth.cluster.uncertainty <= 0.05)

        ggPoint(
            x = log10(plt_df$nFrags + 1),
            y = log10(plt_df$TSSEnrichment + 1),
            size = 0.1,
            colorDensity = TRUE,
            continuousSet = "sambaNight",
            xlabel = "log10(unique fragments + 1)",
            ylabel = "log10(TSS Enrichment + 1)",
            textFamily = "Arial",
            rastr = FALSE
        ) + 
        geom_hline(yintercept = log10(min(df.TSS$TSSEnrichment + 1)), linetype = "dashed") +
        geom_vline(xintercept = min(log10(df.depth$nFrags + 1)), linetype = "dashed") +
        ggtitle(paste0("QC thresholds:\n",i)) 
        ggsave(file.path(outdir, paste0(i,"_QC.png")), dpi=300, width = 5,height = 5)
    }
}

##########################################################################
# Seurat-inspired dot plot for gene activity and integrated RNA expression
##########################################################################                       
custom_DotPlot <- function(gene_list, genedata, celldata, genescoredata, clusters, cluster_levels, legend_label) {
    cluster_size <- as.data.frame(table(celldata[[clusters]]))
    colnames(cluster_size) <- c("Cluster","Cluster_size")
    
    # Get row indices of gene score activity matrix from gene list
    idx <- as.numeric(rownames(genedata[genedata$name %in% gene_list, ]))

    genescoredataframe <- list()
    for(i in unique(celldata[[clusters]])) {
        tmp <- data.frame(avg_exp = Matrix::rowMeans(genescoredata[idx, rownames(celldata)[celldata[[clusters]]==i]]),
                          pct_exp = Matrix::rowSums(genescoredata[idx, rownames(celldata)[celldata[[clusters]]==i]]>0) / cluster_size$Cluster_size[cluster_size$Cluster == i] * 100,
                          gene = genedata[genedata$name %in% gene_list, "name"],
                          Cluster = i)
        genescoredataframe[[i]] <- tmp
    }
    genescoredataframe <- bind_rows(genescoredataframe, .id="Cluster")
    genescoredataframe <- genescoredataframe %>% mutate(gene = factor(gene, levels=gene_list), Cluster = factor(Cluster, levels=cluster_levels), log1p_exp = log(avg_exp + 1))

    plt <- ggplot(data = genescoredataframe, mapping = aes(x = gene, y = Cluster)) +
            geom_point(mapping = aes(size = pct_exp, color = log1p_exp)) +
            scale_colour_gradient(limits = c(0, max(genescoredataframe$log1p_exp)), low="lightgrey", high="blue") +
            scale_size(range=c(0,6), limits = c(0, max(genescoredataframe$pct_exp))) +
            cowplot::theme_cowplot() +
            labs(
                x = NULL,
                y = NULL
            ) +
            guides(color = guide_colorbar(title=legend_label), size = guide_legend(title="Percent\nExpressed")) +
            theme(axis.text.x = element_text(size=14, angle=35, hjust=1),
                    axis.text.y = element_text(size=14))

    return(plt)
}

                       
#-----------------------------------------------------------------------------------------------------------------------------------
# Peak to Gene related functions
#-----------------------------------------------------------------------------------------------------------------------------------

##################################################################################
# Get overlapping genes in Reactome term gene set for peak-to-gene linkage heatmap
##################################################################################               
get_reactome_overlap <- function(reactome_list, name, term) {
    tmp <- sapply(reactome_list, function(df) {
        df %>%
            filter(get(name) %in% term) %>% # use 'get' to use variable in filter
            pull(userId) %>% 
            sapply(function(x) unlist(strsplit(x,split=';',fixed=TRUE))) %>% 
            unlist() %>% 
            unname() %>% 
            unique()})
    return(unique(unlist(tmp)))    
}

###############################################
# Prepare ATAC and RNA matrices for P2G heatmap
###############################################               

apply_limit = function(mat, limits = c(-2,2)) {
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)

    mat <- (mat - min(limits)) / (max(limits) - min(limits))

    return(mat)  
}

#########################################################
# Convert P2G linkages to bedpe for plotgardener plotting
#########################################################  
p2g_to_bedpe = function(links) {
    links_bedpe <- links %>% 
                    as.data.frame() %>% 
                    dplyr::select(seqnames, start, end, value) %>%
                    mutate(start1 = start, end1 = start + 1,
                            start2 = end, end2 = end + 1)
    links_bedpe <- cbind(links_bedpe[,c("seqnames","start1","end1")], links_bedpe[,c("seqnames","start2","end2","value")])
    colnames(links_bedpe)[1] <- "chrom1"
    colnames(links_bedpe)[4] <- "chrom2"
    links_bedpe_GI <- GenomicInteractions(GRanges(setNames(links_bedpe[,1:3], c("chrom","start","end"))),
                                            GRanges(setNames(links_bedpe[,4:6], c("chrom","start","end"))),
                                            counts = links_bedpe[,7],
                                            mode = "strict")
    return(list(bedpe = links_bedpe, bedpe_GI = links_bedpe_GI))
}
                  
#-----------------------------------------------------------------------------------------------------------------------------------
# Cicero related functions
#-----------------------------------------------------------------------------------------------------------------------------------

####################################################
# Make sparse peak matrix for Cicero (make_atac_cds)
####################################################               
make_sparse_peakmatrix = function(mat) {
    peakmatrix_sparse = as.data.frame(summary(assay(mat)))
    peakmatrix_sparse = cbind(peakmatrix_sparse, 
                                    as.data.frame(rowRanges(mat))[peakmatrix_sparse$i, ],
                                    cell = rownames(colData(mat))[peakmatrix_sparse$j])
    peakmatrix_sparse = peakmatrix_sparse %>% mutate(range_string = paste(seqnames, start, end, sep="_"))
    return(peakmatrix_sparse)
}
                  
##########################################
# Prepare cell data set objects for Cicero
##########################################               
prepare_for_cicero = function(cds, proj, embedding) {
    cds <- cds[Matrix::rowSums(exprs(cds)) != 0,] 
    cds <- detect_genes(cds)
    cds <- estimate_size_factors(cds)
    reducedDims(cds)[[embedding]] <- getEmbedding(proj, embedding)
    return(cds)
}

##################################################
# Examine effect of varying k parameter for Cicero 
##################################################
k_param_vis <- function(cds, coords, k_range=seq(10,80,5)) {
    cds.cicero.klist <- lapply(k_range, function(x) make_cicero_cds(cds, reduced_coordinates = coords, k = x, return_agg_info = T))
    names(cds.cicero.klist) <- k_range

    cds.cicero.klist.df <- lapply(cds.cicero.klist, function(x) as.data.frame(x[2]))
    cds.cicero.klist.df <- lapply(cds.cicero.klist.df, function(x) table(x$cell, x$agg_cell))
    cds.cicero.klist.df <- lapply(cds.cicero.klist.df, function(x) rowMeans(x, na.rm=T))
    names(cds.cicero.klist.df) <- k_range
                                  
    plt <- boxplot(cds.cicero.klist.df, ylab="Fraction of Total Aggregates Per Cell", xlab="# Cells Per Aggregate", main="Aggregate Cell Diversity")
    return(list(klist = cds.cicero.klist, boxplt = plt))
}
                                  
############
# Run Cicero 
############
run_cicero = function(cicero_obj, genome.df, age, outdir) {
    conns <- cicero::run_cicero(cicero_obj, genomic_coords = genome.df, window = 5e5, sample_num = 100)
    ccans <- cicero::generate_ccans(conns)

    save(conns, ccans, file=file.path(outdir, paste0(age,"_MuSC_cicero_conns_ccans.RData")))
}

########################################
# Annotate Cicero connections with genes 
########################################
add_genes_to_conns <- function(cds, conns) {
    conns$Peak1_gene <- fData(cds)$gene[match(conns$Peak1, fData(cds)$site_name)]
    conns$Peak2_gene <- fData(cds)$gene[match(conns$Peak2, fData(cds)$site_name)]
    conns$CCAN <- young_ccans$CCAN[match(conns$Peak1, young_ccans$Peak)]
    return(conns)
}
                                  

############################################
# Rank genes by number of Cicero connections 
############################################
get_DORCs <- function(conns, coaccess_thresh) {
    # only using Peak1_genes ensures no mirrored duplicates from conn file
    DORC <- na.omit(conns[which(conns$coaccess>coaccess_thresh),'Peak1_gene'])
    DORC_multiple <- grep(",", DORC, fixed=T, value=T)
    print(sprintf('%d out of %d promoters marked by peaks are marked by the same peak', 
                length(table(DORC_multiple)), 
                length(table(DORC)))) # 1,303 out of 11,895 promoters marked by peaks are marked by the same peak
    DORC <- DORC[!(DORC %in% DORC_multiple)]
    DORC_multiple <- unname(unlist(sapply(DORC_multiple, function(x) unlist(strsplit(x, split=",", fixed=T)))))
    DORC <- c(DORC, DORC_multiple)
    DORC_df <- as.data.frame(sort(table(DORC)))
    return(DORC_df)
}
                                          
###################################
# Make Cicero connections GI object 
###################################                                 
make_conns_gi <- function(conns) {
    conns_peak1 <- Signac::StringToGRanges(conns$Peak1, sep = c("_", "_"))
    conns_peak2 <- Signac::StringToGRanges(conns$Peak2, sep = c("_", "_"))
    conns_gi <- GInteractions(conns_peak1, conns_peak2, mode="strict")
    conns_gi$coaccess <- conns$coaccess
    conns_gi <- unique(conns_gi)
    conns_gi$distance <- pairdist(conns_gi, type="mid")
    conns_gi$idx <- seq_along(conns_gi)

    conns_gi$Peak1_idx <- anchorOne(conns_gi) %>% as.data.frame() %>% 
        mutate(peak_idx = as.numeric(factor(paste(seqnames,start,end,sep="_")))) %>% 
        pull(peak_idx)
    conns_gi$Peak2_idx <- anchorTwo(conns_gi) %>% as.data.frame() %>% 
        mutate(peak_idx = as.numeric(factor(paste(seqnames,start,end,sep="_")))) %>% 
        pull(peak_idx)

    return(conns_gi)
}
 
#########################################
# Get density for Cicero connections plot
#########################################
getDensity <- function(x = NULL, y = NULL, n = 100, sample = NULL, densityMax = 0.95){
  #modified from http://slowkow.com/notes/ggplot2-color-by-density/
  df <- data.frame(x=x,y=y)
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  df$Density <- dens$z[ii]
  df$Density[df$Density > quantile(unique(df$Density),densityMax)] <- quantile(unique(df$Density),densityMax) #make sure the higher end doesnt bias colors
  if(!is.null(sample)){
    df <- df[sample(nrow(df), min(sample,nrow(df))),]
  }
  return(df)
}

#############################
# Annotate Cicero connections 
#############################                              
add_link_annotation <- function(conns_gi, thresh) {
    coaccess_idx <- which(conns_gi$coaccess>thresh)

    peak1_anno <- annotatePeak(anchorOne(conns_gi)[coaccess_idx,], 
                            tssRegion = c(-1000,1000), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")
    peak2_anno <- annotatePeak(anchorTwo(conns_gi)[coaccess_idx,], 
                            tssRegion = c(-1000,1000), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

    conns_gi$Peak1_anno <- "None"
    conns_gi$Peak2_anno <- "None"

    conns_gi$Peak1_anno[coaccess_idx] <- sapply(peak1_anno@anno$annotation, function(x) ifelse(grepl("Intron",x,fixed=T),"Intron",
                                                                                         ifelse(grepl("Exon",x,fixed=T),"Exon",x)))
    conns_gi$Peak2_anno[coaccess_idx] <- sapply(peak2_anno@anno$annotation, function(x) ifelse(grepl("Intron",x,fixed=T),"Intron",
                                                                                         ifelse(grepl("Exon",x,fixed=T),"Exon",x)))

    conns_gi$link_anno <- "Other"
    distal_anno_names <- c("Distal Intergenic", "3' UTR", "5' UTR", "Downstream (<=300bp)")

    conns_gi$link_anno[conns_gi$Peak1_anno=="Promoter" & conns_gi$Peak2_anno=="Promoter"] <- "Promoter-Promoter"
    conns_gi$link_anno[(conns_gi$Peak1_anno %in% distal_anno_names) & (conns_gi$Peak2_anno %in% distal_anno_names)] <- "Distal-Distal"
    conns_gi$link_anno[((conns_gi$Peak1_anno %in% distal_anno_names) & (conns_gi$Peak2_anno=="Promoter")) | 
                    ((conns_gi$Peak2_anno %in% distal_anno_names) & (conns_gi$Peak1_anno=="Promoter"))] <- "Promoter-Distal"

    return(conns_gi)
}

###############################################################################
# Annotate Cicero connection anchors by whether they sit completely within TADs 
###############################################################################
add_conns_TAD <- function(conns_gi, anchor, TAD_domains) {
    if(anchor=="first") {
        gr_a <- first(conns_gi)
    } else if(anchor=="second") {
        gr_a <- second(conns_gi)
    }
    gr_b <- TAD_domains
    
    # generate paired object for all connections
    hits <- findOverlaps(gr_a, gr_b, type='any')
    pairs <- Pairs(gr_a, gr_b, hits = hits)
    olap <- pintersect(pairs)
    # measure overlap between anchor and TAD domain
    mcols(pairs)$overlap_width <- width(olap)
    mcols(pairs)$overlap_fraction <- width(olap)/width(first(pairs))
    mcols(pairs)$idx <- mcols(conns_gi)$idx[from(hits)]
    mcols(pairs)$TAD <- gr_b$ID[to(hits)]
    # add "None" as TAD to anchors not completely sitting in a TAD domain    
    non_olap_idx <- which(mcols(pairs)$overlap_fraction < 1)
    mcols(pairs)$TAD[non_olap_idx] <- "None"
    
    return(pairs)
}
                                                
# Now count the number of TADs per anchor
tally_TADs_per_anchor <- function(conns_gi, conns_gi_a1, conns_gi_a2) {
    conns_gi_a1_agg <- mcols(conns_gi_a1) %>% as.data.frame() %>% dplyr::filter(TAD!="None") %>% group_by(idx) %>% summarise(TAD_agg = TAD)
    conns_gi_a2_agg <- mcols(conns_gi_a2) %>% as.data.frame() %>% dplyr::filter(TAD!="None") %>% group_by(idx) %>% summarise(TAD_agg = TAD)

    mcols(conns_gi)$a1_TAD <- "None"
    mcols(conns_gi)$a2_TAD <- "None"

    mcols(conns_gi)$a1_TAD[match(conns_gi_a1_agg$idx, mcols(conns_gi)$idx)] <- conns_gi_a1_agg$TAD_agg
    mcols(conns_gi)$a2_TAD[match(conns_gi_a2_agg$idx, mcols(conns_gi)$idx)] <- conns_gi_a2_agg$TAD_agg

    mcols(conns_gi)$same_TAD <- mcols(conns_gi) %>% as.data.frame() %>% mutate(same_TAD = ifelse((a1_TAD==a2_TAD) & (a1_TAD!="None") & (a2_TAD!="None"), TRUE, FALSE)) %>% pull(same_TAD)
  
    return(conns_gi)
}

####################################
# Calculate intra TAD ratio of conns 
####################################                                            
calculate_intra_tad_ratio <- function(thresh_range, bins, conns_gi) {
    # create matrix for Intra-TAD enrichment per coaccess threshold
    intra_tad_enrichment <- matrix(nrow = length(thresh_range), ncol = length(bins))
    # 
    for(thresh in seq_along(thresh_range)) {
        for(bin in seq_along(bins)) {
            intra_tad_enrichment[thresh, bin] <- 
                mcols(conns_gi) %>%
                as.data.frame() %>%
                dplyr::filter(coaccess>thresh_range[thresh] & distance_bin==bins[bin]) %>%
                summarise(intra_tad_ratio = sum(same_TAD)/sum(!same_TAD)) %>%
                pull(intra_tad_ratio)
        }  
    }
    rownames(intra_tad_enrichment) <- paste('thresh', thresh_range, sep='_')
    colnames(intra_tad_enrichment) <- bins 

    intra_tad_enrichment <- intra_tad_enrichment %>%
        as.data.frame() %>%
        tibble::rownames_to_column('Coaccess') %>% 
        pivot_longer(cols = !Coaccess, names_to="distance_bins", values_to="ratio") %>%
        mutate(distance_bins = factor(distance_bins, levels=bins))

    return(intra_tad_enrichment)
}

# Set up parallel processing
setup_parallel_processing <- function(n_cores) {
    my.cluster <- makeCluster(n_cores, type = "FORK")
    registerDoParallel(cl = my.cluster)
    return(my.cluster)
}    
                                                
close_parallel_processing <- function(cluster) {
    stopCluster(cluster)
}
                                                
# Shuffle Cicero connections to create a background for analysis
permute_conns <- function(conns_gi, TAD_domains, thresh_range, n_perm, n_cores) {
    my.cluster <- setup_parallel_processing(n_cores)
    stopifnot(getDoParRegistered())
    
    shuffled_ratio_mat_replicated <- foreach(
      i = 1:n_perm, 
      .combine = 'rbind'
    ) %dopar% {
        
        shuffled_ratio_mat <- matrix(nrow = length(levels(conns_gi$distance_bin)), ncol = length(thresh_range))
        
        for(thresh in seq_along(thresh_range)) {
            print(sprintf('Permutation: %d; Thresh: %f', i, thresh_range[thresh]))
            
            subset_conns <- conns_gi[which(conns_gi$coaccess>thresh_range[thresh]), ]
            subset_conns$idx <- 1:length(subset_conns)
            shuffled_subset <- mcols(subset_conns) %>%
                as.data.frame() %>%
                mutate(group_left = paste(seqnames(first(subset_conns)), distance_bin, sep="_"),
                       group_right = paste(seqnames(second(subset_conns)), distance_bin, sep="_"))
            shuffled_subset_idx_left <- shuffled_subset %>%
                group_by(group_left) %>%
                mutate(shuffle_idx = sample(idx, size = n(), replace = FALSE)) %>%
                pull(shuffle_idx)
            shuffled_subset_idx_right <- shuffled_subset %>%
                group_by(group_right) %>%
                mutate(shuffle_idx = sample(idx, size = n(), replace = FALSE)) %>%
                pull(shuffle_idx)
            shuffled_gi <- GInteractions(anchor1=anchorOne(subset_conns)[shuffled_subset_idx_left, ],
                                         anchor2=anchorTwo(subset_conns)[shuffled_subset_idx_right, ])
            shuffled_gi <- GInteractions(anchor1=anchorOne(subset_conns)[shuffled_subset_idx_left, ],
                                         anchor2=anchorTwo(subset_conns)[shuffled_subset_idx_right, ])
            mcols(shuffled_gi)$idx <- 1:length(shuffled_gi)
            mcols(shuffled_gi)$coaccess <- 1 # set fake coaccess for calculate_intra_tad_ratio function
            mcols(shuffled_gi)$distance_bin_left <- as.character(shuffled_subset$distance_bin)[shuffled_subset_idx_left]
            mcols(shuffled_gi)$distance_bin_right <- as.character(shuffled_subset$distance_bin)[shuffled_subset_idx_right]
            mcols(shuffled_gi)$distance_bin_equal <- mcols(shuffled_gi)$distance_bin_left == mcols(shuffled_gi)$distance_bin_right
            shuffled_gi <- shuffled_gi[mcols(shuffled_gi)$distance_bin_equal, ]

            mcols(shuffled_gi)$distance_bin <- mcols(shuffled_gi)$distance_bin_left
            mcols(shuffled_gi)$actual_distance <- pairdist(shuffled_gi, type="mid")

            shuffled_gi_a1 <- add_conns_TAD(shuffled_gi, "first", TAD_domains)
            shuffled_gi_a2 <- add_conns_TAD(shuffled_gi, "second", TAD_domains)
            shuffled_gi <- tally_TADs_per_anchor(shuffled_gi, shuffled_gi_a1, shuffled_gi_a2)

            shuffled_ratio <- calculate_intra_tad_ratio(0, levels(conns_gi$distance_bin), shuffled_gi)
            shuffled_ratio_mat[,thresh] <- shuffled_ratio$ratio
        }

        rownames(shuffled_ratio_mat) <- levels(conns_gi$distance_bin)
        colnames(shuffled_ratio_mat) <- thresh_range

        shuffled_ratio_mat %>%
            t() %>%
            as.data.frame() %>%
            tibble::rownames_to_column('thresh') %>%
            pivot_longer(cols = !thresh, names_to="distance_bin", values_to="ratio") %>%
            mutate(replicate = i)
    }
    return(shuffled_ratio_mat_replicated)
    
    close_parallel_processing(my.cluster)
}

#######################
# Max-matching of CCANs 
#######################  
          
# https://github.com/cole-trapnell-lab/cicero/issues/15

find_best_ccan_match <- function(ccans1, ccans2) {
    names(ccans1) <- c("Peak", "CCAN1")
    names(ccans2) <- c("Peak", "CCAN2")
    comps_mem <- merge(ccans1, ccans2, all=TRUE)
    comps_mem$value <- 1

    comps_mem$CCAN1 <- paste0("C1_", comps_mem$CCAN1)
    comps_mem$CCAN2 <- paste0("C2_", comps_mem$CCAN2)

    comps_mem_temp <- subset(comps_mem, CCAN1 != "C1_NA" | CCAN2 != "C2_NA")
    comps_mem_dc <- reshape2::dcast(CCAN1 ~ CCAN2, value.var = "value", data = comps_mem_temp)

    comps_mem_dc <- subset(comps_mem_dc, CCAN1 != "C1_NA")
    comps_mem_dc <- comps_mem_dc[,1:(ncol(comps_mem_dc) - 1)]
    row.names(comps_mem_dc) <- comps_mem_dc$CCAN1
    comps_mem_dc$CCAN1 <- NULL

    comps_mem_dc_melt <- comps_mem_dc
    comps_mem_dc_melt$CCAN1 <- row.names(comps_mem_dc_melt)
    comps_mem_dc_melt <- reshape2::melt(comps_mem_dc_melt, id.vars = "CCAN1")
    names(comps_mem_dc_melt) <- c("from", "to", "weight")
    comps_mem_dc_melt <- subset(comps_mem_dc_melt, from != to)
    grph <- igraph::graph_from_data_frame(subset(comps_mem_dc_melt, weight > 0), directed = FALSE)
    igraph::V(grph)$type <- grepl("orig", igraph::V(grph))
    mbm <- maxmatching::maxmatching(grph)

    matching <- data.frame(CCAN1 = mbm$matching, CCAN2 = names(mbm$matching))
    matching <- subset(matching, grepl("C1", CCAN1) & grepl("C2", CCAN2))
    matching$CCAN1 <- gsub("C1_", "", matching$CCAN1)
    matching$CCAN2 <- gsub("C2_", "", matching$CCAN2)
    return(matching)
}

#################################
# Get genes that lie within CCANs 
#################################  
get_CCAN_genes <- function(gene_anno_ranges, CCAN_list, ccans) {
  output <- list()
  for(i in CCAN_list) {
    tmp <- ccans[ccans$CCAN==i,] %>% 
      tidyr::separate(Peak, into=c("seqnames","start","end"), sep="_", remove=F) %>%
      dplyr::select(seqnames,start,end) %>%
      GRanges() %>%
      sort()
    output[[i]] <- subsetByOverlaps(gene_anno_ranges, 
                                   ranges = GRanges(data.frame(seqnames=as.character(unique(seqnames(tmp))),
                                                                                 start=min(start(tmp)),
                                                                                 end=max(end(tmp)))))
  }
  return(output)
}

#################################
# CCANs within provided features
#################################  
ccans_in_feature <- function(ccans_ranges, domain_ranges) {
    gr_a <- ccans_ranges
    gr_b <- domain_ranges
    hits <- findOverlaps(gr_a, gr_b, type="within", ignore.strand = TRUE)
    gr_a$ID <- NA
    gr_a$ID[from(hits)] <- gr_b$ID[to(hits)]
    gr_a$CCAN <- ccans_ranges$CCAN
    return(gr_a)
}
                                                
#########################################
# Get peaks within target CCAN for plotting
#########################################  
get_ccan_subset <- function(conns, ccan, target_ccan) {
    ccan_subset <- conns[conns$Peak1 %in% ccan$Peak[ccan$CCAN==target_ccan], ] %>%
        separate(Peak1, into=c("seqnames1","start1","end1"), sep="_", remove=F) %>%
        separate(Peak2, into=c("seqnames2","start2","end2"), sep="_", remove=F)
    ccan_subset_GI = GenomicInteractions(anchor1 = GRanges(setNames(ccan_subset[,c("seqnames1","start1","end1")], c("seqnames","start","end"))),
                                         anchor2 = GRanges(setNames(ccan_subset[,c("seqnames2","start2","end2")], c("seqnames","start","end"))),
                                         counts = ccan_subset$coaccess,
                                         mode = "strict")
    ccan_subset_GI$distance <- pairdist(ccan_subset_GI)
    ccan_subset_GI$h <- ccan_subset_GI$distance / max(ccan_subset_GI$distance, na.rm=T)
    return(ccan_subset_GI)
}

###############################
# Plot peaks within target CCAN
###############################  
make_ccans_plt <- function(ccan_subset, coaccess_ylim, p2g_palette, params, y, flip) {
    tmp <- as.data.frame(ccan_subset) %>% 
        dplyr::select(seqnames1, start1, end1, seqnames2, start2, end2, counts, distance, h) %>%
        mutate(seqnames1 = as.character(seqnames1), seqnames2 = as.character(seqnames2))
    
    ccans_plt <- plotPairsArches(tmp, height=1, y=y, flip = flip,
                                 fill=colorby("counts", palette=colorRampPalette(p2g_palette),
                                              range = c(0.1, coaccess_ylim)),
                                   archHeight = 'h', linecolor = "fill", alpha=0.6, params=params)

    annoHeatmapLegend(
      plot = ccans_plt, fontcolor = "black",
      ticks = FALSE, digits=2,
      params=params,
      x=as.numeric(gsub("inches","",ccans_plt$width,fixed=TRUE)) + 0.5,
      y=as.numeric(gsub("inches","",ccans_plt$y,fixed=TRUE)) + 0.125,
      width = 0.10, height = 0.75, fontsize = 10)
    
    return(ccans_plt)
}
