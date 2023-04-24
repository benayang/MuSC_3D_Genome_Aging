library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(rstatix)
library(ggpubr)
library(ggsci)
library(pheatmap)
library(RColorBrewer)

projdir="C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer\\TAD_interactions\\promoter_H3K4me3_ATAC"

# Get expression data -----------------------------------------------------

tpm.data = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC", 
                                "compartmentExpression", "RNA_transformed_tpm_minus_sva_contribs.txt"))
gene.sym = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
gene.sym.tbl = data.frame(table(gene.sym)) # find duplicate gene symbols
duplicate.genes = gene.sym.tbl[gene.sym.tbl$Freq>1,"gene.sym"]

tpm.data = tpm.data[-which(gene.sym %in% duplicate.genes), 1:8] # isolate day0 data for young and aged and remove gene duplicates
rownames(tpm.data) = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
tpm.data = tpm.data[complete.cases(tpm.data), ] # remove NA rows
tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) all(x>=0))), ] # keep only genes with TPM>=0 in either young or aged

DE.genes = unlist(read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC", 
                                       "compartmentExpression", "DE_0.05padj_log2(1.5)LFC_genenames.csv"), sep="\t", as.is=T))
DE.genes.name = unlist(lapply(DE.genes, function(x) unlist(strsplit(x,"_",fixed=T))[2]), use.names=F)


# Get gene lists ----------------------------------------------------------

# EXCLUDING 2 YOUNG REPLICATES
filter.tpm.data = tpm.data[,c(1:3,5:7)]
avg.tpm.data = data.frame(A=rowMeans(filter.tpm.data[,c(1:3)]), Y=rowMeans(filter.tpm.data[,c(4:6)])) %>% tibble::rownames_to_column("gene")

tasks = c("tss.1kb.negH3K4me3.bed", "tss.1kb.posH3K4me3.posATAC.bed", "tss.1kb.posH3K4me3.negATAC.bed")

gene_list = list()
for (age in c("young","aged")) {
  for (taskname in tasks) {
    tmp = read.table(file.path(projdir,paste(age,taskname,sep=".")))
    colnames(tmp) = c("chr","start","end","gene")
    tmp$Age = age
    tmp$Group = taskname
    gene_list[[taskname]] = rbind(gene_list[[taskname]], tmp)
  }
}
gene_list_df = do.call(rbind,gene_list)
gene_list_df_expressed = gene_list_df %>% filter(gene %in% rownames(tpm.data))

gene_list_df_expressed = cbind(gene_list_df_expressed, avg.tpm.data[gene_list_df_expressed$gene,])
gene_list_df_expressed = gene_list_df_expressed %>% mutate(log2A = log2(A+1), 
                                                           log2Y = log2(Y+1),
                                                           Group = sapply(Group,function(x) ifelse(x=="tss.1kb.negH3K4me3.bed","H3K4me3-",
                                                                                                   ifelse(x=="tss.1kb.posH3K4me3.posATAC.bed","H3K4me3+/ATAC+",
                                                                                                          ifelse(x=="tss.1kb.posH3K4me3.negATAC.bed","H3K4me3+/ATAC-",NA)))))

# Filter data frame so that ages of H3K4me3 and ATAC peak set match RNA data 
gene_list_df_expressed_filtered = gene_list_df_expressed %>% 
  pivot_longer(cols=c("A","Y"), names_to="TPM.age", values_to="TPM") %>%
  filter((TPM.age=="A" & Age=="aged") | (TPM.age=="Y" & Age=="young")) %>%
  mutate(Age = factor(Age, levels=c("young","aged"))) 
gene_test = gene_list_df_expressed_filtered %>% group_by(Group) %>% wilcox_test(TPM~Age)
gene_test = gene_test %>% add_xy_position(x="Group",fun="mean_se")

ggplot(gene_list_df_expressed_filtered, aes(x=Group, y=TPM)) +
  stat_summary(aes(group=Age), fun.data="mean_se", geom="errorbar", width=0.5, position=position_dodge(0.9)) +
  stat_summary(aes(fill=Age), fun.data="mean_se", geom="bar", color="black", position=position_dodge(0.9)) +
  stat_pvalue_manual(gene_test, tip.length=0.0001,label.size=5) +
  scale_x_discrete(labels=c("H3K4me3-","H3K4me3+\nATAC-","H3K4me3+\nATAC+")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme_bw() +
  labs(x=NULL,y="TPM") +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,color="black",face="bold"),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.position = "top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
ggsave(file.path(projdir,"group_TPM.png"),dpi=300,width=4.5,height=5)

gene_list_df_expressed_filtered %>%
  group_by(Age,Group) %>%
  summarise(count=n())

# Get top 2000 expressed genes per age per promoter group
top2000 = gene_list_df_expressed_filtered %>%
  group_by(Age,Group) %>%
  arrange(desc(TPM), .by_group=T) %>%
  dplyr::slice(1:2000) %>%
  distinct() %>%
  drop_na()
top2000 %>% group_by(Age,Group) %>% summarise(count=n())
write.table(top2000,file.path(projdir,"top2000_expressed_genes_per_group.txt"),sep="\t",quote=F,row.names=F)

# look only at DE genes
DEtop2000 = gene_list_df_expressed_filtered %>%
  filter(gene %in% DE.genes.name) %>%
  filter(!grepl("Gm",gene,fixed=T) & !grepl("Rik",gene,fixed=T)) %>%
  group_by(Age,Group) %>%
  arrange(desc(TPM), .by_group=T) %>%
  dplyr::slice(1:2000) %>%
  distinct() %>%
  drop_na()
DEtop2000 %>% group_by(Age,Group) %>% summarise(count=n())
write.table(DEtop2000,file.path(projdir,"DEtop2000_expressed_genes_per_group.txt"),sep="\t",quote=F,row.names=F)


# Get interaction data ----------------------------------------------------

# interaction_list = list()
# for (age in c("young","aged")) {
#   for (taskname in interaction.tasks) {
#     tmp = read.table(file.path(projdir,"promoter_interactions",paste(age,taskname,sep="_")), header=T, sep=",")
#     tmp$Age = age
#     tmp$Group = ifelse(grepl("negH3K4me3",taskname,fixed=T),"H3K4me3-",
#                        ifelse(grepl("posATAC",taskname,fixed=T),"H3K4me3+/ATAC+",
#                               ifelse(grepl("negATAC",taskname,fixed=T),"H3K4me3+/ATAC-",NA)))
#     interaction_list[[taskname]] = rbind(interaction_list[[taskname]], tmp)
#   }
# }
# interaction_list_df = do.call(rbind,interaction_list) %>% distinct() %>% drop_na()
# head(interaction_list_df)

tasks = c("agedDETop2000_posATAC_posH3K4me3_1kbPromoter_interaction_obs_counts.csv",
          "youngDETop2000_posATAC_posH3K4me3_1kbPromoter_interaction_obs_counts.csv")
interaction_stat_list = list()
for(taskname in tasks) {
  fname = gsub("_posATAC_posH3K4me3_1kbPromoter_interaction_obs_counts.csv","",taskname,fixed=T,)
  tmp = read.table(file.path(projdir,"promoter_interactions","obs_counts",taskname), sep=",", header=T)
  tmp$Group = fname
  interaction_stat_list[[fname]] = tmp
}
interaction_stat_df = bind_rows(interaction_stat_list)
interaction_stat_df = left_join(interaction_stat_df, avg.tpm.data, by="gene")
#interaction_stat_df = cbind(interaction_stat_df, avg.tpm.data[interaction_stat_df$gene,])
#head(interaction_stat_df)

# see if there's a relationship between interaction number/strength and gene expression
ggplot(interaction_stat_df, aes(x=log10(Young_mean+1),y=log2(Y+1))) +
  geom_hex(bins=50)
ggplot(interaction_stat_df, aes(x=Young_num,y=log2(Y+1))) +
  geom_hex(bins=50)
ggplot(interaction_stat_df, aes(x=log10(Aged_mean+1),y=log2(A+1))) +
  geom_hex(bins=50)
ggplot(interaction_stat_df, aes(x=Aged_num,y=log2(A+1))) +
  geom_hex(bins=50)

ggplot(interaction_stat_df, aes(x=Young_num,y=Aged_num)) +
  facet_grid(cols=vars(Group), labeller=labeller(Group=c("agedDETop2000"="Aged Genes","youngDETop2000"="Young Genes"))) +
  geom_hex(bins=50) +
  geom_abline(slope=1,intercept=0,lty=2,lwd=0.5) +
  scale_fill_distiller(palette="RdBu") +
  scale_x_continuous(limits=c(0,1600)) +
  scale_y_continuous(limits=c(0,1600)) +
  coord_equal() +
  theme_bw() + 
  labs(x="# Young Interactions", y="# Aged Interactions") +
  guides(fill=guide_colorbar(title="Count")) +
  theme(legend.position = "right",
        #legend.text = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text=element_text(size=12,face="bold",color="black"),
        axis.title=element_text(size=12,face="bold",color="black"),
        strip.text=element_text(size=12,face="bold"))
ggsave(file.path(projdir,"promoter_interactions","num_interactions_top2000_DEgenes.png"),dpi=300,width=5,height=4)

ggplot(interaction_stat_df, aes(x=log10(Young_mean),y=log10(Aged_mean))) +
  facet_grid(cols=vars(Group), labeller=labeller(Group=c("agedDETop2000"="Aged Genes","youngDETop2000"="Young Genes"))) +
  geom_hex(bins=50) +
  geom_abline(slope=1,intercept=0,lty=2,lwd=0.5) +
  scale_fill_distiller(palette="RdBu") +
  scale_x_continuous(limits=c(-0.5,1.25)) +
  scale_y_continuous(limits=c(-0.5,1.25)) +
  coord_equal() +
  theme_bw() +
  labs(x=expression(bold(paste("lo","g"["10"],"(Mean Young Counts)"))), y=expression(bold(paste("lo","g"["10"],"(Mean Aged Counts)")))) +
  guides(fill=guide_colorbar(title="Count")) +
  theme(legend.position = "right",
        #legend.text = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text=element_text(size=12,face="bold",color="black"),
        axis.title=element_text(size=12,face="bold",color="black"),
        strip.text=element_text(size=12,face="bold"))
ggsave(file.path(projdir,"promoter_interactions","mean_interactions_top2000_DEgenes.png"),dpi=300,width=5,height=4)


# Get interaction data with H3K4me3/ATAC marks ----------------------------


tasks = c("agedTop2000_agedHiC_posATACposH3K4me3_TAD_1kb_obs_counts_with_marks.csv",
          "agedTop2000_youngHiC_posATACposH3K4me3_TAD_1kb_obs_counts_with_marks.csv",
          "youngTop2000_agedHiC_posATACposH3K4me3_TAD_1kb_obs_counts_with_marks.csv",
          "youngTop2000_youngHiC_posATACposH3K4me3_TAD_1kb_obs_counts_with_marks.csv")
interaction_list = list()
for(taskname in tasks) {
  fname = gsub("_posATACposH3K4me3_TAD_1kb_obs_counts_with_marks.csv","",taskname,fixed=T,)
  tmp = read.table(file.path(projdir,"promoter_interactions","obs_counts",taskname), sep=",", header=T)
  tmp$Group = fname
  tmp = tmp %>% separate(Group, into=c("promoterAge","HiCAge"), sep="_", remove=F)
  interaction_list[[fname]] = tmp
}
young_hic_interactions = bind_rows(interaction_list[c("agedTop2000_youngHiC","youngTop2000_youngHiC")])
aged_hic_interactions = bind_rows(interaction_list[c("agedTop2000_agedHiC","youngTop2000_agedHiC")])

# determine if binX or binY are in the given promoter region
bin_class = function(start,end,x) ifelse(x>=start & x<=end, 'Promoter', 'Non-promoter')
young_hic_interactions$Young_binX_class = bin_class(young_hic_interactions['start'],young_hic_interactions['end'],young_hic_interactions['Young_binX'])
young_hic_interactions$Young_binY_class = bin_class(young_hic_interactions['start'],young_hic_interactions['end'],young_hic_interactions['Young_binY'])
young_hic_interactions$Aged_binX_class = bin_class(young_hic_interactions['start'],young_hic_interactions['end'],young_hic_interactions['Aged_binX'])
young_hic_interactions$Aged_binY_class = bin_class(young_hic_interactions['start'],young_hic_interactions['end'],young_hic_interactions['Aged_binY'])
aged_hic_interactions$Young_binX_class = bin_class(aged_hic_interactions['start'],aged_hic_interactions['end'],aged_hic_interactions['Young_binX'])
aged_hic_interactions$Young_binY_class = bin_class(aged_hic_interactions['start'],aged_hic_interactions['end'],aged_hic_interactions['Young_binY'])
aged_hic_interactions$Aged_binX_class = bin_class(aged_hic_interactions['start'],aged_hic_interactions['end'],aged_hic_interactions['Aged_binX'])
aged_hic_interactions$Aged_binY_class = bin_class(aged_hic_interactions['start'],aged_hic_interactions['end'],aged_hic_interactions['Aged_binY'])

# saveRDS(young_hic_interactions, file.path(projdir,"promoter_interactions","young_hic_DEgenes_obs_interactions.RDS"))
# saveRDS(aged_hic_interactions, file.path(projdir,"promoter_interactions","aged_hic_DEgenes_obs_interactions.RDS"))

# young_hic_interactions = readRDS(file.path(projdir,"promoter_interactions","young_hic_DEgenes_interactions.RDS"))
# aged_hic_interactions = readRDS(file.path(projdir,"promoter_interactions","aged_hic_DEgenes_interactions.RDS"))

# now exclude bins that are in the promoter region. We're only interested in the regions interacting with promoters, not the promoters themselves.
aged_hic_interactions = aged_hic_interactions %>% filter((Aged_binX_class != Aged_binY_class) & (Young_binX_class != Young_binY_class))
young_hic_interactions = young_hic_interactions %>% filter((Aged_binX_class != Aged_binY_class) & (Young_binX_class != Young_binY_class))

# Get all data for non-promoter bins
cols = c(colnames(aged_hic_interactions)[1:4], "ID", "Group", "promoterAge", "HiCAge")
aged_hic_df = bind_rows(aged_hic_interactions %>% filter(Aged_binX_class=="Non-promoter") %>%
                          dplyr::select(c(all_of(cols),contains("Aged_binX"),contains("counts")) & !contains("class")) %>% 
                          mutate(markAge="Aged") %>% 
                          dplyr::rename(H3K4me3=contains("H3K4me3"), ATAC=contains("ATAC"), tss_1kb=contains("tss_1kb"), non_promoter_TAD=contains("TAD_name"), 
                                        TAD_start=contains("TAD_start"), TAD_end=contains("TAD_end"), TAD_score=contains("TAD_score"), Counts = contains("counts"),
                                        non_promoter_hic_bin=Aged_binX),
                        aged_hic_interactions %>% filter(Aged_binY_class=="Non-promoter") %>%
                          dplyr::select(c(all_of(cols),contains("Aged_binY"),contains("counts")) & !contains("class")) %>% 
                          mutate(markAge="Aged") %>% 
                          dplyr::rename(H3K4me3=contains("H3K4me3"), ATAC=contains("ATAC"), tss_1kb=contains("tss_1kb"), non_promoter_TAD=contains("TAD_name"), 
                                        TAD_start=contains("TAD_start"), TAD_end=contains("TAD_end"), TAD_score=contains("TAD_score"), Counts = contains("counts"), 
                                        non_promoter_hic_bin=Aged_binY),
                        aged_hic_interactions %>% filter(Young_binX_class=="Non-promoter") %>%
                          dplyr::select(c(all_of(cols),contains("Young_binX"),contains("counts")) & !contains("class")) %>% 
                          mutate(markAge="Young") %>% 
                          dplyr::rename(H3K4me3=contains("H3K4me3"), ATAC=contains("ATAC"), tss_1kb=contains("tss_1kb"), non_promoter_TAD=contains("TAD_name"), 
                                        TAD_start=contains("TAD_start"), TAD_end=contains("TAD_end"), TAD_score=contains("TAD_score"), Counts = contains("counts"),
                                        non_promoter_hic_bin=Young_binX),
                        aged_hic_interactions %>% filter(Young_binY_class=="Non-promoter") %>%
                          dplyr::select(c(all_of(cols),contains("Young_binY"),contains("counts")) & !contains("class")) %>% 
                          mutate(markAge="Young") %>% 
                          dplyr::rename(H3K4me3=contains("H3K4me3"), ATAC=contains("ATAC"), tss_1kb=contains("tss_1kb"), non_promoter_TAD=contains("TAD_name"), 
                                        TAD_start=contains("TAD_start"), TAD_end=contains("TAD_end"), TAD_score=contains("TAD_score"), Counts = contains("counts"),
                                        non_promoter_hic_bin=Young_binY)) %>%
  arrange(ID)

cols = c(colnames(young_hic_interactions)[1:4], "ID", "Group", "promoterAge", "HiCAge")
young_hic_df = bind_rows(young_hic_interactions %>% filter(Aged_binX_class=="Non-promoter") %>%
                           dplyr::select(c(all_of(cols),contains("Aged_binX"),contains("counts")) & !contains("class")) %>% 
                           mutate(markAge="Aged") %>% 
                           dplyr::rename(H3K4me3=contains("H3K4me3"), ATAC=contains("ATAC"), tss_1kb=contains("tss_1kb"), non_promoter_TAD=contains("TAD_name"), 
                                  TAD_start=contains("TAD_start"), TAD_end=contains("TAD_end"), TAD_score=contains("TAD_score"), Counts = contains("counts"),
                                  non_promoter_hic_bin=Aged_binX),
                         young_hic_interactions %>% filter(Aged_binY_class=="Non-promoter") %>%
                           dplyr::select(c(all_of(cols),contains("Aged_binY"),contains("counts")) & !contains("class")) %>% 
                           mutate(markAge="Aged") %>% 
                           dplyr::rename(H3K4me3=contains("H3K4me3"), ATAC=contains("ATAC"), tss_1kb=contains("tss_1kb"), non_promoter_TAD=contains("TAD_name"), 
                                  TAD_start=contains("TAD_start"), TAD_end=contains("TAD_end"), TAD_score=contains("TAD_score"), Counts = contains("counts"),
                                  non_promoter_hic_bin=Aged_binY),
                         young_hic_interactions %>% filter(Young_binX_class=="Non-promoter") %>%
                           dplyr::select(c(all_of(cols),contains("Young_binX"),contains("counts")) & !contains("class")) %>% 
                           mutate(markAge="Young") %>% 
                           dplyr::rename(H3K4me3=contains("H3K4me3"), ATAC=contains("ATAC"), tss_1kb=contains("tss_1kb"), non_promoter_TAD=contains("TAD_name"), 
                                  TAD_start=contains("TAD_start"), TAD_end=contains("TAD_end"), TAD_score=contains("TAD_score"), Counts = contains("counts"),
                                  non_promoter_hic_bin=Young_binX),
                         young_hic_interactions %>% filter(Young_binY_class=="Non-promoter") %>%
                           dplyr::select(c(all_of(cols),contains("Young_binY"),contains("counts")) & !contains("class")) %>% 
                           mutate(markAge="Young") %>% 
                           dplyr::rename(H3K4me3=contains("H3K4me3"), ATAC=contains("ATAC"), tss_1kb=contains("tss_1kb"), non_promoter_TAD=contains("TAD_name"), 
                                  TAD_start=contains("TAD_start"), TAD_end=contains("TAD_end"), TAD_score=contains("TAD_score"), Counts = contains("counts"),
                                  non_promoter_hic_bin=Young_binY)) %>%
  arrange(ID)

# Now do the same for all promoter bins, focusing on TADs
cols = c(colnames(aged_hic_interactions)[1:4], "ID", "Group", "promoterAge", "HiCAge")
aged_hic_promoter_df = bind_rows(aged_hic_interactions %>% filter(Aged_binX_class=="Promoter") %>% 
                                   dplyr::select(c(all_of(cols), "Aged_binX", "Aged_binX_TAD_name")) %>% 
                                   mutate(markAge="Aged") %>% dplyr::rename(promoter_hic_bin=Aged_binX, promoter_TAD=Aged_binX_TAD_name),
                                 aged_hic_interactions %>% filter(Aged_binY_class=="Promoter") %>% 
                                   dplyr::select(c(all_of(cols), "Aged_binY", "Aged_binY_TAD_name")) %>% 
                                   mutate(markAge="Aged") %>% dplyr::rename(promoter_hic_bin=Aged_binY, promoter_TAD=Aged_binY_TAD_name),
                                 aged_hic_interactions %>% filter(Young_binX_class=="Promoter") %>% 
                                   dplyr::select(c(all_of(cols), "Young_binX", "Young_binX_TAD_name")) %>% 
                                   mutate(markAge="Young") %>% dplyr::rename(promoter_hic_bin=Young_binX, promoter_TAD=Young_binX_TAD_name),
                                 aged_hic_interactions %>% filter(Young_binY_class=="Promoter") %>% 
                                   dplyr::select(c(all_of(cols), "Young_binY", "Young_binY_TAD_name")) %>% 
                                   mutate(markAge="Young") %>% dplyr::rename(promoter_hic_bin=Young_binY, promoter_TAD=Young_binY_TAD_name)) %>%
  arrange(ID)

cols = c(colnames(young_hic_interactions)[1:4], "ID", "Group", "promoterAge", "HiCAge")
young_hic_promoter_df = bind_rows(young_hic_interactions %>% filter(Aged_binX_class=="Promoter") %>% 
                                    dplyr::select(c(all_of(cols), "Aged_binX", "Aged_binX_TAD_name")) %>% 
                                    mutate(markAge="Aged") %>% dplyr::rename(promoter_hic_bin=Aged_binX, promoter_TAD=Aged_binX_TAD_name),
                                  young_hic_interactions %>% filter(Aged_binY_class=="Promoter") %>% 
                                    dplyr::select(c(all_of(cols), "Aged_binY", "Aged_binY_TAD_name")) %>% 
                                    mutate(markAge="Aged") %>% dplyr::rename(promoter_hic_bin=Aged_binY, promoter_TAD=Aged_binY_TAD_name),
                                  young_hic_interactions %>% filter(Young_binX_class=="Promoter") %>% 
                                    dplyr::select(c(all_of(cols), "Young_binX", "Young_binX_TAD_name")) %>% 
                                    mutate(markAge="Young") %>% dplyr::rename(promoter_hic_bin=Young_binX, promoter_TAD=Young_binX_TAD_name),
                                  young_hic_interactions %>% filter(Young_binY_class=="Promoter") %>% 
                                    dplyr::select(c(all_of(cols), "Young_binY", "Young_binY_TAD_name")) %>% 
                                    mutate(markAge="Young") %>% dplyr::rename(promoter_hic_bin=Young_binY, promoter_TAD=Young_binY_TAD_name)) %>%
  arrange(ID)

aged_hic_TAD_df = left_join(aged_hic_df, aged_hic_promoter_df, by=c(colnames(aged_hic_interactions)[1:4], "ID", "Group", "promoterAge", "HiCAge", "markAge"))
young_hic_TAD_df = left_join(young_hic_df, young_hic_promoter_df, by=c(colnames(young_hic_interactions)[1:4], "ID", "Group", "promoterAge", "HiCAge", "markAge"))

interaction_df = rbind(aged_hic_TAD_df, young_hic_TAD_df)
interaction_df = interaction_df %>% mutate(TAD_class = ifelse(non_promoter_TAD=="." | promoter_TAD==".", "Asymmetric-TAD",
                                                              ifelse(non_promoter_TAD==promoter_TAD, "Intra-TAD", "Inter-TAD")))
interaction_df = left_join(interaction_df, avg.tpm.data, by="gene")
# saveRDS(interaction_df, file.path(projdir,"promoter_interactions","interaction_DEgenes_df.RDS"))
# interaction_df = readRDS(file.path(projdir,"promoter_interactions","interaction_DEgenes_df.RDS"))
with(interaction_df, table(TAD_class, Group, markAge))
head(interaction_df)

interaction_df %>% 
  group_by(promoterAge, HiCAge, markAge) %>% 
  summarise(across(c("H3K4me3","ATAC","tss_1kb"), ~sum(.x>0)/n())) %>%
  filter((grepl("aged",promoterAge,fixed=T) & markAge=="Aged") | (grepl("young",promoterAge,fixed=T) & markAge=="Young")) %>%
  pivot_longer(cols=c("H3K4me3","ATAC","tss_1kb"),names_to="Mark",values_to="Frac") %>%
  ggplot(aes(x=promoterAge,y=Frac)) +
  facet_grid(cols=vars(HiCAge), labeller=labeller(HiCAge=c("agedHiC"="Aged", "youngHiC"="Young"))) +
  geom_col(aes(fill=Mark)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_x_discrete(labels=c("Aged","Young")) +
  theme_bw() +
  guides(fill=guide_legend(title="")) +
  labs(x="Promoter Age", y="Feature Fraction") +
  theme(axis.text=element_text(size=14, color="black", face="bold"),
        axis.title=element_text(size=14, color="black", face="bold"),
        strip.text=element_text(size=14, color="black", face="bold"),
        legend.text=element_text(size=12),
        legend.position="top")


# Interaction classification plots ----------------------------------------

interaction_df_TAD_test = interaction_df %>%
  filter((markAge=="Aged" & HiCAge=="agedHiC") | (markAge=="Young" & HiCAge=="youngHiC")) %>%
  mutate(log10counts = log10(Counts+1)) %>%
  group_by(promoterAge, HiCAge) %>%
  t_test(log10counts ~ TAD_class)
interaction_df_TAD_test = interaction_df_TAD_test %>% add_xy_position(x="HiCAge",group="TAD_class")
interaction_df_Age_test = interaction_df %>%
  filter((markAge=="Aged" & HiCAge=="agedHiC") | (markAge=="Young" & HiCAge=="youngHiC")) %>%
  mutate(log10counts = log10(Counts+1)) %>%
  group_by(promoterAge, TAD_class) %>%
  t_test(log10counts ~ HiCAge)
interaction_df_Age_test = interaction_df_Age_test %>% add_xy_position(x="HiCAge", group="TAD_class", dodge=0.9)

interaction_df %>%
  filter((markAge=="Aged" & HiCAge=="agedHiC") | (markAge=="Young" & HiCAge=="youngHiC")) %>%
  ggplot(aes(x=HiCAge, y=log10(Counts+1))) +
  facet_grid(cols=vars(promoterAge)) +
  geom_boxplot(aes(fill=TAD_class), outlier.size=0.1, position=position_dodge(0.9), color="black") +
  # geom_bracket(xmin=with(interaction_df_TAD_test, x[Group=="AtoB" & Mark=="ATAC"]),
  #              xmax=with(interaction_df_TAD_test, x[Group=="BtoA" & Mark=="ATAC"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="AtoB" & Mark=="ATAC"])+0.5,
  #              label=7.63e-8, tip.length = 0.01) +
  # stat_pvalue_manual(interaction_df_TAD_test, tip.length=0.01) +
  # stat_pvalue_manual(interaction_df_Age_test, tip.length=0.01, bracket.nudge.y = 0.5, step.group.by = "HiCAage") +
  theme_bw() +
  labs(x=NULL) +
  guides(fill=guide_legend(nrow=2)) +
  theme(axis.text=element_text(size=14, color="black", face="bold"),
        axis.title=element_text(size=14, color="black", face="bold"),
        strip.text=element_text(size=14, color="black", face="bold"),
        legend.text=element_text(size=12),
        legend.position="top")
ggsave(file.path(projdir,"promoter_interactions","interaction_counts_TADs_boxplot.png"),dpi=300,width=5,height=4)


# Plots as a function of distance -----------------------------------------

interaction_df %>% 
  mutate(distance=abs(non_promoter_hic_bin-promoter_hic_bin)) %>% 
  filter(markAge=="Aged" & HiCAge=="agedHiC") %>%
  #filter((markAge=="Aged" & HiCAge=="agedHiC") | (markAge=="Young" & HiCAge=="youngHiC")) %>%
  ggplot(aes(x=distance, y=log10(Counts+1))) +
  facet_grid(cols=vars(promoterAge)) +
  geom_point(aes(color=TAD_class), alpha=0.6, size=0.5) +
  #geom_boxplot(aes(fill=TAD_class), outlier.size=0.1, position=position_dodge(0.9), color="black") +
  theme_bw()

# Plot TAD score per group ------------------------------------------------

interaction_df %>% filter(promoterAge=="youngTop2000") %>%
  ggplot(aes(x=HiCAge,y=as.numeric(TAD_score))) +
  facet_grid(cols=vars(TAD_class), rows=vars(markAge)) +
  geom_boxplot() +
  stat_compare_means(method="t.test")
ggsave(file.path(projdir,"promoter_interactions","young_tss_TAD_score.png"),dpi=300,width=5,height=5)


# enumerate interaction classifications -----------------------------------

library(ComplexHeatmap)
library(Hmisc)
library(scales)

interaction_df %>% 
  group_by(promoterAge, HiCAge, markAge, TAD_class) %>%
  summarise(count=n())

combMat_list = list()
for(i in unique(interaction(interaction_df$promoterAge, interaction_df$HiCAge, interaction_df$markAge, interaction_df$TAD_class))) {
  pa = unlist(strsplit(i,split='.',fixed=T))[1]
  pa_name = capitalize(gsub('DETop2000','',pa))
  ha = unlist(strsplit(i,split='.',fixed=T))[2]
  ma = unlist(strsplit(i,split='.',fixed=T))[3]
  tc = unlist(strsplit(i,split='.',fixed=T))[4]
  if(pa_name != ma){
    # only look at interactions where promoter and marks are from the same Age
    next
  }
  tmp = interaction_df %>% 
    filter(promoterAge==pa & HiCAge==ha & markAge==ma & TAD_class==tc) %>%
    #mutate(IntraTAD = ifelse(TAD_class=="Inter-TAD",0,1)) %>%
    select(c("H3K4me3","ATAC","tss_1kb")) %>%
    mutate(across(everything(),~ as.numeric(.x>0))) %>% make_comb_mat()
  set_name(tmp)[set_name(tmp)=='tss_1kb'] = "Promoter"
  combMat_list[[i]] = tmp
}
saveRDS(combMat_list,file.path(projdir,"promoter_interactions","combMat_list.RDS"))

# plot individual UpSet plots for each group
for(n in names(combMat_list)) {
  cs = comb_size(combMat_list[[n]])
  max_cs = max(unlist(lapply(combMat_list, function(x) max(comb_size(x)))))
  ss = set_size(combMat_list[[n]])
  max_ss = max(unlist(lapply(combMat_list, function(x) max(set_size(x)))))
  breaks = seq(0,6,1)
  
  png(file.path(projdir,"promoter_interactions","UpSet","IndividualPlots",paste0(n,"_DEgenes_upset.png")),res=300,units="in",width=5,height=4)
  ht = UpSet(combMat_list[[n]],
             set_order = c("Promoter", "ATAC", "H3K4me3"),
             top_annotation = HeatmapAnnotation(
               "Feature Intersections" = anno_barplot(cs,
                                                       ylim = c(0, max_cs*1.1),
                                                       border = FALSE, 
                                                       gp = gpar(fill = "black"), 
                                                       height = unit(4, "cm")
               ),
               annotation_name_side = "left",
               annotation_name_rot = 90),
             right_annotation = rowAnnotation(
               "Feature Size" = anno_barplot(ss,
                                             ylim = c(0, max_ss*1.1),
                                             border = FALSE, 
                                             gp = gpar(fill = "black"), 
                                             height = unit(3, "cm"),
                                             width = unit(2, "cm")
               )))
  ht = draw(ht)
  od = column_order(ht)
  decorate_annotation("Feature Intersections", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
              default.units = "native", just = c("left", "bottom"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
  })
  dev.off()
}

# Plot stacked UpSet plots for TADs
for (i in unique(interaction(interaction_df$promoterAge, interaction_df$markAge))) {
  pa = unlist(strsplit(i,split='.',fixed=T))[1]
  pa_name = capitalize(gsub('DETop2000','',pa))
  ma = unlist(strsplit(i,split='.',fixed=T))[2]
  if(pa_name != ma){
    # only look at interactions where promoter and marks are from the same Age
    next
  }
  young_interTAD=combMat_list[[paste(pa,'youngHiC',ma,"Inter-TAD",sep=".")]]
  young_intraTAD=combMat_list[[paste(pa,'youngHiC',ma,"Intra-TAD",sep=".")]]
  young_asymmTAD=combMat_list[[paste(pa,'youngHiC',ma,"Asymmetric-TAD",sep=".")]]
  aged_interTAD=combMat_list[[paste(pa,'agedHiC',ma,"Inter-TAD",sep=".")]]
  aged_intraTAD=combMat_list[[paste(pa,'agedHiC',ma,"Intra-TAD",sep=".")]]
  aged_asymmTAD=combMat_list[[paste(pa,'agedHiC',ma,"Asymmetric-TAD",sep=".")]]
  
  sum_cs_young=sum(comb_size(young_interTAD)+comb_size(young_intraTAD)+comb_size(young_asymmTAD))
  sum_cs_aged=sum(comb_size(aged_interTAD)+comb_size(aged_asymmTAD)+comb_size(aged_asymmTAD))

  young_interTAD_comb_sets = lapply(comb_name(young_interTAD[comb_degree(young_interTAD)>0]), function(nm) extract_comb(unlist(young_interTAD), nm))
  young_interTAD_counts = lapply(young_interTAD_comb_sets, 
                                 function(x) mean(subset(interaction_df, promoterAge==pa & HiCAge=="youngHiC" & markAge==ma & TAD_class=="Inter-TAD", "Counts")[x,]))
  young_intraTAD_comb_sets = lapply(comb_name(young_intraTAD[comb_degree(young_intraTAD)>0]), function(nm) extract_comb(unlist(young_intraTAD), nm))
  young_intraTAD_counts = lapply(young_intraTAD_comb_sets, 
                              function(x) mean(subset(interaction_df, promoterAge==pa & HiCAge=="youngHiC" & markAge==ma & TAD_class=="Intra-TAD", "Counts")[x,]))
  young_asymmTAD_comb_sets = lapply(comb_name(young_asymmTAD[comb_degree(young_asymmTAD)>0]), function(nm) extract_comb(unlist(young_asymmTAD), nm))
  young_asymmTAD_counts = lapply(young_asymmTAD_comb_sets, 
                                 function(x) mean(subset(interaction_df, promoterAge==pa & HiCAge=="youngHiC" & markAge==ma & TAD_class=="Asymmetric-TAD", "Counts")[x,]))
  aged_interTAD_comb_sets = lapply(comb_name(aged_interTAD[comb_degree(aged_interTAD)>0]), function(nm) extract_comb(unlist(aged_interTAD), nm))
  aged_interTAD_counts = lapply(aged_interTAD_comb_sets, 
                             function(x) mean(subset(interaction_df, promoterAge==pa & HiCAge=="agedHiC" & markAge==ma & TAD_class=="Inter-TAD", "Counts")[x,]))
  aged_intraTAD_comb_sets = lapply(comb_name(aged_intraTAD[comb_degree(aged_intraTAD)>0]), function(nm) extract_comb(unlist(aged_intraTAD), nm))
  aged_intraTAD_counts = lapply(aged_intraTAD_comb_sets, 
                             function(x) mean(subset(interaction_df, promoterAge==pa & HiCAge=="agedHiC" & markAge==ma & TAD_class=="Intra-TAD", "Counts")[x,]))
  aged_asymmTAD_comb_sets = lapply(comb_name(aged_asymmTAD[comb_degree(aged_asymmTAD)>0]), function(nm) extract_comb(unlist(aged_asymmTAD), nm))
  aged_asymmTAD_counts = lapply(aged_asymmTAD_comb_sets, 
                                 function(x) mean(subset(interaction_df, promoterAge==pa & HiCAge=="agedHiC" & markAge==ma & TAD_class=="Asymmetric-TAD", "Counts")[x,]))
  
  
  # comb_sets = lapply(comb_sets, function(gr) {
  #   gr$TAD_score = as.numeric(interaction_df[gr,]$TAD_score)
  #   gr$Aged_TPM = as.numeric(interaction_df[gr,]$A)
  #   gr$Young_TPM = as.numeric(interaction_df[gr,]$Y)
  #   gr$Counts = as.numeric(interaction_df[gr,]$Counts)
  #   gr
  # })
  
  ht = UpSet(combMat_list[[1]][comb_degree(combMat_list[[1]])>0], # exclude group with no overlaps
             set_order = c("Promoter", "ATAC", "H3K4me3"),
             height = unit(3,"cm"),
             top_annotation=HeatmapAnnotation(
               "Young Fraction\nof Interactions" = anno_barplot(t(data.frame(rbind(comb_size(young_asymmTAD)[-8]/sum_cs_young,
                                                                             comb_size(young_interTAD)[-8]/sum_cs_young, 
                                                                                   comb_size(young_intraTAD)[-8]/sum_cs_young), 
                                                                             row.names = c("Asymmetric-TAD","Inter-TAD","Intra-TAD"))),
                                                                axis_param = list(at = seq(0,0.15,0.025)),
                                                                height = unit(2.75, "cm"),
                                                                ylim = c(0, 0.15*1.1),
                                                                gp = gpar(fill = c("white","grey","black"))),
               "Aged Fraction\nof Interactions" = anno_barplot(t(data.frame(rbind(comb_size(aged_asymmTAD)[-8]/sum_cs_aged,
                                                                                  comb_size(aged_interTAD)[-8]/sum_cs_aged, 
                                                                                  comb_size(aged_intraTAD)[-8]/sum_cs_aged), 
                                                                            row.names = c("Asymmetric-TAD","Inter-TAD","Intra-TAD"))),
                                                               axis_param = list(at = seq(0,0.15,0.025)),
                                                               height = unit(2.75, "cm"),
                                                               ylim = c(0, 0.15*1.1),
                                                               gp = gpar(fill = c("white","grey","black"))),
               #annotation_label=expression(paste("lo","g"["10"],"(Intersections)")),
               gap = unit(2, "mm"),
               annotation_name_side = "left",
               annotation_name_rot = 90),
             right_annotation = rowAnnotation(
               "Young Set" = anno_barplot(t(rbind(set_size(young_asymmTAD), set_size(young_interTAD), set_size(young_intraTAD))),
                                          axis_param = list(at = seq(0,2.5e5,5e4), labels = c(0,"5e+04","1e+05","1.5e+05","2e+05","2.5e+05","3e+05"), side="top"),
                                          height = unit(2, "cm"),
                                          width = unit(2, "cm"),
                                          ylim = c(0, 2.5e5*1.1),
                                          gp = gpar(fill = c("white","grey","black"))),
               "Aged Set" = anno_barplot(t(rbind(set_size(aged_asymmTAD), set_size(aged_interTAD), set_size(aged_intraTAD))),
                                         axis_param = list(at = seq(0,2.5e5,5e4), labels = c(0,"5e+04","1e+05","1.5e+05","2e+05","2.5e+05","3e+05"), side="top"),
                                         height = unit(2, "cm"),
                                         width = unit(2, "cm"),
                                         ylim = c(0, 2.5e5*1.1),
                                         gp = gpar(fill = c("white","grey","black"))),
               annotation_name_side = "top",
               gap = unit(2, "mm")),
             bottom_annotation = HeatmapAnnotation(
               "Young Counts\n(Obs/Exp)" = anno_barplot(t(data.frame(rbind(unlist(young_asymmTAD_counts), unlist(young_interTAD_counts), unlist(young_intraTAD_counts)))),
                                                        axis_param = list(at = seq(0,340,40)),
                                                        ylim = c(0, 300*1.1),
                                                        gp = gpar(fill = c("white","grey","black")),
                                                        height = unit(2, "cm")),
               "Aged Counts\n(Obs/Exp)" = anno_barplot(t(data.frame(rbind(unlist(aged_asymmTAD_counts), unlist(aged_interTAD_counts), unlist(aged_intraTAD_counts)))),
                                                       axis_param = list(at = seq(0,340,40)),
                                                       ylim = c(0, 300*1.1),
                                                       gp = gpar(fill = c("white","grey","black")),
                                                       height = unit(2, "cm")),
               gap = unit(2, "mm")
             )
  )
  png(file.path(projdir,"promoter_interactions","UpSet",paste0(i,"_stackedTAD_stackedAge_DEgenes_upset.png")),res=300,units="in",width=5,height=6)
  draw(ht, annotation_legend_side = "left")
  draw(Legend(labels = c("Asymmetric-TAD","Inter-TAD","Intra-TAD"), legend_gp = gpar(fill = c("white","grey","black"), fontsize=14), title = "TAD Group", border="black"),
       x = unit(0.95, "npc"), y = unit(0.85, "npc"), just = c("right", "top"))
  dev.off()
}


# Make virtual 4C plot ----------------------------------------------------

gene_list_4C = c("Ddit3","Tcf7l2","Sox2","Myod1","Prdm16")

young_hic_4C = interaction_df %>% 
  filter(gene %in% gene_list_4C, HiCAge=="youngHiC", markAge=="Young", promoterAge=="youngTop2000") %>%
  select(-ID) %>%
  distinct()
aged_hic_4C = interaction_df %>% 
  filter(gene %in% gene_list_4C, HiCAge=="agedHiC", markAge=="Aged", promoterAge=="agedTop2000") %>%
  select(-ID) %>%
  distinct()

bin_range = range(c(young_hic_4C$non_promoter_hic_bin, aged_hic_4C$non_promoter_hic_bin))
hic_4C_df = data.frame(bins=seq(bin_range[1],bin_range[2],by=5000))

for (g in gene_list_4C) {
  hic_4C_df[paste("young",g,sep="_")] = young_hic_4C[young_hic_4C$gene==g, "Counts"][match(hic_4C_df$bins,young_hic_4C[young_hic_4C$gene==g,]$non_promoter_hic_bin)]
  hic_4C_df[paste("aged",g,sep="_")] = aged_hic_4C[aged_hic_4C$gene==g, "Counts"][match(hic_4C_df$bins,aged_hic_4C[aged_hic_4C$gene==g,]$non_promoter_hic_bin)]
}

vline_df = unique(aged_hic_4C[aged_hic_4C$gene %in% gene_list_4C,c("gene","start","end")])

hic_4C_df %>%
  pivot_longer(cols=!contains("bin"),names_to="AgeGene",values_to="counts") %>%
  separate(AgeGene, into=c("Age","gene"), sep="_") %>%
  ggplot(aes(x=bins, y=counts)) +
  facet_grid(rows=vars(gene)) +
  geom_point(aes(color=Age), alpha=0.6, size=0.5) +
  #geom_vline(data=vline_df, aes(xintercept=start)) +
  #geom_vline(data=vline_df, aes(xintercept=end))+
  scale_color_manual(values=pal_nejm()(2)) +
  theme_bw() +
  labs(x="Genomic Bins (5kb)", y="Obs/Exp Counts") +
  theme(legend.position="top",
        axis.text=element_text(face="bold",size=10,color="black"),
        axis.title=element_text(face="bold",size=10))
ggsave(file.path(projdir,"promoter_interactions","virtual_4C_plots.png"),dpi=300,width=4,height=5)  

# Export bedpe files for individual genes ---------------------------------

interaction_df %>% 
  filter(gene=="Ddit3", HiCAge=="youngHiC", markAge=="Young", promoterAge=="youngDETop2000") %>%
  mutate(non_promoter_hic_bin_end = non_promoter_hic_bin+5000, chr2=chr) %>%
  select(chr,start,end,chr2,non_promoter_hic_bin,non_promoter_hic_bin_end,ID,Counts) %>%
  distinct() %>%
  write.table(file.path(projdir,"promoter_interactions","Ddit3_young_interactions.bedpe"), sep="\t", row.names=F, col.names=F, quote=F)
  

# H3K4me3
aged_hic_interactions %>%
  pivot_longer(cols=contains("H3K4me3"),names_to="H3K4me3_group",values_to="H3K4me3_count") %>%
  group_by(H3K4me3_group) %>%
  summarise(num = sum(H3K4me3_count>0),
            frac = num/n())
young_hic_interactions %>%
  pivot_longer(cols=contains("H3K4me3"),names_to="H3K4me3_group",values_to="H3K4me3_count") %>%
  group_by(H3K4me3_group) %>%
  summarise(num = sum(H3K4me3_count>0),
            frac = num/n())
# ATAC
aged_hic_interactions %>%
  pivot_longer(cols=contains("ATAC"),names_to="ATAC_group",values_to="ATAC_count") %>%
  group_by(ATAC_group) %>%
  summarise(num = sum(ATAC_count>0),
            frac = num/n())
young_hic_interactions %>%
  pivot_longer(cols=contains("ATAC"),names_to="ATAC_group",values_to="ATAC_count") %>%
  group_by(ATAC_group) %>%
  summarise(num = sum(ATAC_count>0),
            frac = num/n())
# Promoters
aged_hic_interactions %>%
  pivot_longer(cols=contains("tss_1kb"),names_to="tss_1kb_group",values_to="tss_1kb_count") %>%
  group_by(tss_1kb_group) %>%
  summarise(num = sum(tss_1kb_count>0),
            frac = num/n())
young_hic_interactions %>%
  pivot_longer(cols=contains("tss_1kb"),names_to="tss_1kb_group",values_to="tss_1kb_count") %>%
  group_by(tss_1kb_group) %>%
  summarise(num = sum(tss_1kb_count>0),
            frac = num/n())

young_hic_interactions

ggplot(interaction_stat_df, aes(x=Young_counts, y=Aged_counts)) +
  facet_grid(cols=vars(Age)) +
  geom_point(alpha=0.5) +
  geom_abline(slope=1,intercept=0,lty=2) +
  coord_fixed() +
  #scale_x_continuous(limits=c(0,9)) +
  #scale_y_continuous(limits=c(0,9)) +
  scale_color_jama() +
  theme_bw() +
  theme(legend.position="top")
ggsave(file.path(projdir,"sum_counts_per_contact_scatterplot.png"),dpi=300,width=6,height=4)

interaction_list_test = interaction_list_df %>%
  pivot_longer(cols=c(Young_mean,Aged_mean), names_to="mean_age", values_to="mean") %>%
  group_by(Group)  %>%
  wilcox_test(mean~mean_age)
interaction_list_test = interaction_list_test %>% add_xy_position(x="Group",group="mean_age",dodge=0.9,fun="mean_se")

interaction_list_df %>%
  pivot_longer(cols=c(Young_mean,Aged_mean), names_to="mean_age", values_to="mean") %>%
  mutate(mean_age = factor(mean_age, levels=c("Young_mean","Aged_mean"))) %>%
  ggplot(aes(x=Group, y=mean)) +
  stat_summary(aes(group=mean_age),fun.data="mean_se",geom="errorbar",width=0.5,position=position_dodge(0.9)) +
  stat_summary(aes(fill=mean_age),fun.data="mean_se",geom="bar",position=position_dodge(0.9)) + 
  scale_y_continuous(expand=expansion(mult=c(0,0.1)), breaks=seq(0,1.5,0.25)) +
  scale_x_discrete(labels=c("H3K4me3-","H3K4me3+\nATAC-","H3K4me3+\nATAC+")) +
  stat_pvalue_manual(interaction_list_test, tip.length=0.01, bracket.nudge.y = 0.1, label.size=5) +
  theme_bw() +
  scale_fill_manual(values=rev(pal_nejm()(2)), labels=c("Y","A")) +
  guides(fill=guide_legend(title="Age")) +
  #geom_boxplot(aes(group=interaction(Group,sum_age)),position=position_dodge(0.9)) +
  #geom_violin(aes(fill=sum_age), position=position_dodge(0.9)) +
  #geom_boxplot(aes(group=interaction(Group,sum_age)), position=position_dodge(0.9), width=0.1, outlier.shape = NA) +
  labs(x=NULL,
       y="Mean Normalized Counts per Contact") +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,color="black",face="bold"),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.position = "top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
ggsave(file.path(projdir,"mean_count_per_contact.png"),dpi=300,width=4.5,height=5)


interaction_list_test = interaction_list_df %>%
  pivot_longer(cols=c(Young_sum,Aged_sum), names_to="sum_age", values_to="sum") %>%
  group_by(Group)  %>%
  wilcox_test(sum~sum_age)
interaction_list_test = interaction_list_test %>% add_xy_position(x="Group",group="sum_age",dodge=0.9,fun="mean_se")

interaction_list_df %>%
  pivot_longer(cols=c(Young_sum,Aged_sum), names_to="sum_age", values_to="sum") %>%
  mutate(sum_age = factor(sum_age, levels=c("Young_sum","Aged_sum"))) %>%
  ggplot(aes(x=Group, y=sum)) +
  stat_summary(aes(group=sum_age),fun.data="mean_se",geom="errorbar",width=0.5,position=position_dodge(0.9)) +
  stat_summary(aes(fill=sum_age),fun.data="mean_se",geom="bar",position=position_dodge(0.9)) + 
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_x_discrete(labels=c("H3K4me3-","H3K4me3+\nATAC-","H3K4me3+\nATAC+")) +
  stat_pvalue_manual(interaction_list_test, tip.length=0.01, bracket.nudge.y = 0.1, label.size=5) +
  theme_bw() +
  scale_fill_manual(values=rev(pal_nejm()(2)), labels=c("Y","A")) +
  guides(fill=guide_legend(title="Age")) +
  #geom_boxplot(aes(group=interaction(Group,sum_age)),position=position_dodge(0.9)) +
  #geom_violin(aes(fill=sum_age), position=position_dodge(0.9)) +
  #geom_boxplot(aes(group=interaction(Group,sum_age)), position=position_dodge(0.9), width=0.1, outlier.shape = NA) +
  labs(x=NULL,
       y="sum Normalized Counts per Contact") +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,color="black",face="bold"),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.position = "top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
ggsave(file.path(projdir,"sum_count_per_contact.png"),dpi=300,width=4.5,height=5)


# Look at heatmaps --------------------------------------------------------

mat = interaction_list_df %>% filter(Group=="H3K4me3+/ATAC-" & Age=="aged") %>% remove_rownames() %>% column_to_rownames("genes") %>% select(Young_sum, Aged_sum) %>% as.matrix()
young_neg_H3K4me3.plt = pheatmap(mat, cluster_cols=F, show_rownames=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)),
                                 cellwidth=50, cellheight=0.5)
young_neg_H3K4me3.clusters = data.frame(cluster=sort(cutree(young_neg_H3K4me3.plt$tree_row,k=2)))

cluster_rowAnno = data.frame(cluster=young_neg_H3K4me3.clusters[rownames(mat),], row.names=rownames(mat))
cluster_rowAnno$cluster = sapply(cluster_rowAnno$cluster, function(x) paste("Cluster",x,sep=" "))
cluster_rowAnno$cluster = factor(cluster_rowAnno$cluster, ordered=T)

young_neg_H3K4me3.plt = pheatmap(mat, cluster_cols=F, show_rownames=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)),
                                 cellwidth=50, cellheight=0.5, annotation_row = cluster_rowAnno,
                                 filename=file.path(projdir,"aged_posH3K4me3_negATAC_sum_interactions.png"))

young_neg_H3K4me3.tpm = filter.tpm.data[rownames(young_neg_H3K4me3.clusters),]
young_neg_H3K4me3.avg.tpm = avg.tpm.data[rownames(young_neg_H3K4me3.clusters),]
young_neg_H3K4me3.avg.tpm$cluster = factor(young_neg_H3K4me3.clusters[rownames(young_neg_H3K4me3.avg.tpm),])

young_neg_H3K4me3.avg.tpm.test = young_neg_H3K4me3.avg.tpm %>%
  pivot_longer(cols=c(A,Y),names_to="Age",values_to="TPM") %>%
  mutate(log2TPM=log2(TPM+1)) %>%
  group_by(Age) %>%
  wilcox_test(log2TPM~cluster) %>%
  add_xy_position(x="Age",group="cluster",dodge=0.9)

young_neg_H3K4me3.avg.tpm %>%
  pivot_longer(cols=c(A,Y),names_to="Age",values_to="TPM") %>%
  ggplot(aes(x=Age,y=log2(TPM+1))) + 
  geom_boxplot(aes(fill=cluster),color="black",position=position_dodge(0.9)) +
  theme_bw() +
  labs(x=NULL) +
  stat_pvalue_manual(young_neg_H3K4me3.avg.tpm.test, tip.length = 0.01, bracket.nudge.y = 0.5) +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,color="black",face="bold"),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.position = "top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
ggsave(file.path(projdir,"aged_posH3K4me3_negATAC_sum_interactions_expression_boxplot.png"),dpi=300,width=3,height=5)

pheatmap(apply(young_neg_H3K4me3.avg.tpm[,c("A","Y")],2,function(x) log2(x+1)), cluster_cols=F, annotation_row = rowAnno,
         show_rownames=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)),
         cellwidth=50, cellheight=0.5)

rowAnno = data.frame(cluster=young_neg_H3K4me3.clusters[rownames(young_neg_H3K4me3.tpm),], row.names=rownames(young_neg_H3K4me3.tpm))
rowAnno$cluster = sapply(rowAnno$cluster, function(x) paste("Cluster",x,sep=" "))
rowAnno$cluster = factor(rowAnno$cluster, ordered=T)

pheatmap(young_neg_H3K4me3.tpm, cluster_cols=F, annotation_row = rowAnno,
         show_rownames=F, scale="row", color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)),
         cellwidth=50, cellheight=0.5,
         filename=file.path(projdir,"aged_posH3K4me3_negATAC_sum_interactions_expression.png"))


# Merge interactions with expression --------------------------------------

avg.tpm.data$genes = rownames(avg.tpm.data)
interaction_list_expression_df = left_join(interaction_list_df, avg.tpm.data, by="genes") %>%
  pivot_longer(cols=contains("sum"), names_to="sum_age", values_to="sum") %>%
  pivot_longer(cols=contains("mean"), names_to="mean_age", values_to="mean") %>%
  pivot_longer(cols=c("A","Y"), names_to="TPM_age", values_to="TPM")
head(interaction_list_expression_df)

# all promoter classifications
interaction_list_expression_df %>%
  ggplot(aes(x=mean, y=log2(TPM+1))) +
  facet_grid(cols=vars(TPM_age),rows=vars(mean_age),labeller=labeller(TPM_age=c(A="Aged TPM",Y="Young TPM"))) +
  geom_point(aes(color=Group),size=1,alpha=0.25) +
  theme_bw() +
  labs(x="Contact Counts Mean") +
  theme(legend.position="top")
ggsave(file.path(projdir,"expression_vs_mean_by_group.png"),dpi=300,width=6,height=6)

# only H3K4me3+/ATAC+
interaction_list_expression_df %>%
  filter(Group=="H3K4me3+/ATAC+") %>%
  ggplot(aes(x=sum, y=log2(TPM+1))) +
  facet_grid(cols=vars(TPM_age),rows=vars(sum_age),labeller=labeller(TPM_age=c(A="Aged TPM",Y="Young TPM"))) +
  geom_point(size=1,alpha=0.25) +
  theme_bw() +
  labs(x="Contact Counts Sum") +
  theme(legend.position="top")
ggsave(file.path(projdir,"expression_vs_sum_posH3K4me3_posATAC.png"),dpi=300,width=6,height=6)

pheatmap(as.matrix(gene_list_df_expressed[gene_list_df_expressed$Age=="young" & gene_list_df_expressed$Group=="tss.1kb.negH3K4me3.bed", c("log2Y","log2A")]), 
         cluster_cols=F, show_rownames=F)
pheatmap(as.matrix(gene_list_df_expressed[gene_list_df_expressed$Age=="young" & gene_list_df_expressed$Group=="tss.1kb.posH3K4me3.posATAC.bed", c("log2Y","log2A")]), 
         cluster_cols=F, show_rownames=F)
pheatmap(as.matrix(gene_list_df_expressed[gene_list_df_expressed$Age=="young" & gene_list_df_expressed$Group=="tss.1kb.posH3K4me3.negATAC.bed", c("log2Y","log2A")]), 
         cluster_cols=F, show_rownames=F)
