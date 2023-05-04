suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(GenomicInteractions))
suppressPackageStartupMessages(library(HelloRanges))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(svMisc))

suppressPackageStartupMessages(library(progressr))
handlers(global = TRUE)
handlers("progress", "beepr")

registerDoParallel(50)

add_conns_TAD = function(conns_gi, a, TAD_domains) {
  conns_TAD_pairs = eval({
    genome <- Seqinfo(genome = "mm10")
    if(a=="first") {
      gr_a <- anchorOne(conns_gi)
    } else if(a=="second") {
      gr_a <- anchorTwo(conns_gi)
    }
    gr_a$idx = mcols(conns_gi)$idx
    gr_b <- TAD_domains
    hits <- findOverlaps(gr_a, gr_b, type = "within", ignore.strand = TRUE)
    pairs <- pair(gr_a, gr_b, hits, all.x = TRUE)
    olap <- pintersect(pairs, ignore.strand = TRUE)
    mcols(pairs)$overlap_width <- width(olap)
    mcols(pairs)$overlap_fraction <- width(olap)/width(first(pairs))
    mcols(pairs)$idx = first(pairs)$idx
    mcols(pairs)$TAD = second(pairs)$ID
    pairs
  })
  non_olap_idx = which(mcols(pairs)$overlap_fraction < 1)
  mcols(pairs)$TAD[non_olap_idx] = "None"
  return(pairs)
}

calculate_intra_tad_ratio = function(thresh_range, bins, conns_gi) {
  # create matrix for Intra-TAD enrichment per coaccess threshold
  intra_tad_enrichment = matrix(nrow = length(thresh_range), ncol = length(bins))
  # 
  for(thresh in seq_along(thresh_range)) {
    for(bin in seq_along(bins)) {
      intra_tad_enrichment[thresh, bin] = 
        mcols(conns_gi) %>%
        as.data.frame() %>%
        dplyr::filter(coaccess>thresh_range[thresh] & distance_bin==bins[bin]) %>%
        summarise(intra_tad_ratio = sum(same_TAD)/sum(!same_TAD)) %>%
        pull(intra_tad_ratio)
    }  
  }
  rownames(intra_tad_enrichment) = paste('thresh', thresh_range, sep='_')
  colnames(intra_tad_enrichment) = bins 
  
  intra_tad_enrichment = intra_tad_enrichment %>%
    as.data.frame() %>%
    tibble::rownames_to_column('Coaccess') %>% 
    pivot_longer(cols = !Coaccess, names_to="distance_bins", values_to="ratio") %>%
    mutate(distance_bins = factor(distance_bins, levels=bins))
  
  return(intra_tad_enrichment)
}

permute_conns = function(conns_gi, TADs, thresh_range, iters) {
  #p = progressor(along = iters)
  shuffled_ratio_mat_replicated = 
    foreach(i = iters, .combine='rbind') %dopar% {
      #Sys.sleep(6.0-i)
      #p(sprintf("i=%g", i))
      #if(i==101) cat("Done!\n")
      #shuffled_ratio_mat = matrix(nrow = length(levels(conns_gi$distance_bin)), ncol = length(thresh_range))
      print(i)
      foreach(thresh = seq_along(thresh_range), .combine='cbind') %do% {
        #print(paste(thresh_range[thresh], i, sep="_"))

        shuffled_subset = as.data.frame(mcols(conns_gi))
        shuffled_subset = shuffled_subset[which(shuffled_subset$coaccess>thresh_range[thresh]), ]
        
        shuffled_subset_idx_left = shuffled_subset %>%
          group_by(chr, distance_bin) %>%
          mutate(shuffle_idx = sample(idx, size = n(), replace = F)) %>%
          pull(shuffle_idx)
        shuffled_subset_idx_right = shuffled_subset %>%
          group_by(chr, distance_bin) %>%
          mutate(shuffle_idx = sample(idx, size = n(), replace = F)) %>%
          pull(shuffle_idx)    
        shuffled_gi = GenomicInteractions(anchor1=anchorOne(conns_gi)[shuffled_subset_idx_left, ],
                                          anchor2=anchorTwo(conns_gi)[shuffled_subset_idx_right, ])
        mcols(shuffled_gi)$idx = seq_len(length(shuffled_gi))
        mcols(shuffled_gi)$coaccess = 1 # set fake coaccess for calculate_intra_tad_ratio function
        mcols(shuffled_gi)$distance_bin_left = mcols(conns_gi)$distance_bin[shuffled_subset_idx_left]
        mcols(shuffled_gi)$distance_bin_right = mcols(conns_gi)$distance_bin[shuffled_subset_idx_right]
        mcols(shuffled_gi)$distance_bin_equal = identical(mcols(shuffled_gi)$distance_bin_left, mcols(shuffled_gi)$distance_bin_right)
        
        tryCatch(all(mcols(shuffled_gi)$distance_bin_equal), finally=sprintf("%s\n%s",mcols(shuffled_gi)$distance_bin_left,mcols(shuffled_gi)$distance_bin_right))
        #shuffled_gi = shuffled_gi[mcols(shuffled_gi)$distance_bin_equal, ]
        
        mcols(shuffled_gi)$distance_bin = mcols(shuffled_gi)$distance_bin_left[mcols(shuffled_gi)$distance_bin_equal]
        mcols(shuffled_gi)$actual_distance = pairdist(shuffled_gi, type="mid")
        
        shuffled_gi_a1 = add_conns_TAD(shuffled_gi, "first", TADs)
        shuffled_gi_a2 = add_conns_TAD(shuffled_gi, "second", TADs)
        shuffled_gi = anno_conns(shuffled_gi, shuffled_gi_a1, shuffled_gi_a2)

        print(shuffled_gi)

        #shuffled_ratio = calculate_intra_tad_ratio(0, levels(conns_gi$distance_bin), shuffled_gi)
        calculate_intra_tad_ratio(0, levels(conns_gi$distance_bin), shuffled_gi)$ratio
        #shuffled_ratio_mat[,thresh] = shuffled_ratio$ratio
      }
      rownames(shuffled_ratio_mat) = levels(conns_gi$distance_bin)
      colnames(shuffled_ratio_mat) = thresh_range

      shuffled_ratio_mat %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column('thresh') %>%
        pivot_longer(cols = !thresh, names_to="distance_bin", values_to="ratio") %>%
        mutate(replicate = i)
  }
  stopImplicitCluster()
  return(shuffled_ratio_mat_replicated)
}

aged.domain = read.table("/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed")
aged.domain.ranges = GRanges(aged.domain[,1:5] %>% setNames(c("chr","start","end","ID","score")))

young.domain = read.table("/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed")
young.domain.ranges = GRanges(young.domain[,1:5] %>% setNames(c("chr","start","end","ID","score")))

aged_conns_gi = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR/aged_MuSC/aged_conns_gi_peakmatrix.RDS")
young_conns_gi = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR/young_MuSC/young_conns_gi_peakmatrix.RDS")

thresh_range = seq(0.05,0.25,by=0.05)

shuffled_ratio_mat_permute_aged = permute_conns(aged_conns_gi, aged.domain.ranges, thresh_range, 1:100)
shuffled_ratio_mat_permute_young = permute_conns(young_conns_gi, young.domain.ranges, thresh_range, 1:1000)

saveRDS(shuffled_ratio_mat_permute_aged, "/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR/shuffled_ratio_mat_permute_aged.RDS")
saveRDS(shuffled_ratio_mat_permute_young, "/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR/shuffled_ratio_mat_permute_young.RDS")