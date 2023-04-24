library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(reshape2)

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS\\APA"


# 5kb APA -----------------------------------------------------------------

aged = read.table(file.path(projdir, "aged_merged_10kb_gw_normedAPA.txt"))
aged = apply(aged,2,function(x) gsub(",","",x))
aged = apply(aged,2,function(x) gsub("[","",x,fixed=T))
aged = apply(aged,2,function(x) gsub("]","",x,fixed=T))
aged = apply(aged,2,as.numeric)
# breaks = seq(0,7000,1000), legend_breaks = seq(0,7000,1000),
pheatmap(aged, color = colorRampPalette(brewer.pal(9, "Reds"))(500), legend_breaks = seq(0,6,1),
         cluster_rows=F, cluster_cols=F, show_rownames = F, show_colnames = F, border_color = NA)

# APA: 3.049
aged.melt = melt(as.matrix(aged), varnames=c("row","col")) %>% 
  mutate(col=factor(col),
         row=factor(row,levels=rev(unique(row))),
         Age="Aged")

aged_plt = ggplot(aged.melt, aes(x=col, y=row))+ 
  geom_tile(aes(fill=value), color=NA) +
  geom_text(aes(5,20),label="APA=4.07",size=5) +
  geom_text(aes(17,20),label="max=4.95",size=5) +
  geom_text(aes(23,21/2), label="Aged", hjust=0.5, vjust=0.5, angle=270, size=7) +
  scale_fill_gradient(low="white", high="red", limits=c(0,max(aged)), na.value="black") +
  scale_x_discrete(labels=c("-100kb",rep("",9),"0",rep("",9),"+100kb")) +
  scale_y_discrete(labels=c("-100kb",rep("",9),"0",rep("",9),"+100kb")) +
  coord_fixed(clip='off') +
  theme_void() +
  labs(x=NULL,y="Aged") +
  theme(legend.position="none",
        axis.text=element_text(color="black"))

# APA: 3.116
young = read.table(file.path(projdir, "young_merged_10kb_gw_normedAPA.txt"))
young = apply(young,2,function(x) gsub(",","",x))
young = apply(young,2,function(x) gsub("[","",x,fixed=T))
young = apply(young,2,function(x) gsub("]","",x,fixed=T))
young = apply(young,2,as.numeric)

pheatmap(young, color = colorRampPalette(brewer.pal(9, "Reds"))(500), 
         cluster_rows=F, cluster_cols=F, show_rownames = F, show_colnames = F, border_color = NA)

young.melt = melt(as.matrix(young), varnames=c("row","col")) %>% 
  mutate(col=factor(col),
         row=factor(row,levels=rev(unique(row))),
         Age="Young")

young_plt = ggplot(young.melt, aes(x=col, y=row))+ 
  geom_tile(aes(fill=value), color=NA) +
  geom_text(aes(5,20),label="APA=3.95",size=5) +
  geom_text(aes(17,20),label="max=4.76",size=5) +
  scale_fill_gradient(low="white", high="red", limits=c(0,max(young)), na.value="black") +
  scale_x_discrete(labels=c("-100kb",rep("",9),"0",rep("",9),"+100kb")) +
  scale_y_discrete(labels=c("-100kb",rep("",9),"0",rep("",9),"+100kb")) +
  geom_text(aes(23,21/2), label="Young", hjust=0.5, vjust=0.5, angle=270, size=7) +
  coord_fixed(clip='off') +
  theme_void() +
  labs(x=NULL,y=NULL) +
  theme(legend.position="none",
        axis.text=element_text(color="black"))

young_plt / aged_plt
ggsave(file.path(projdir,"10kb_combined_APA.png"),dpi=300,width=4,height=5)

# combined = rbind(young.melt,aged.melt)
# combined$Age = factor(combined$Age, levels=c("Young","Aged"))
# label.df=data.frame(x=c(5,5,17,17), y=20, 
#                     lab=c("APA=3.52","APA=3.58","max=5.18","max=5.47"), 
#                     Age=factor(c("Young","Aged","Young","Aged"), levels=c("Young","Aged")))
# 
# 
# ggplot(combined, aes(x=col, y=row)) +
#   facet_grid(rows=vars(Age)) +
#   geom_tile(aes(fill=value), color=NA) +
#   geom_text(data=label.df, 
#             aes(x=x,y=y,label=lab),
#             size=4) +
#   scale_fill_gradient(low="white", high="red", limits=c(0,max(combined$value)), na.value="black", name=NULL) +
#   scale_x_discrete(labels=c("-100kb",rep("",9),"0",rep("",9),"+100kb")) +
#   scale_y_discrete(labels=c("-100kb",rep("",9),"0",rep("",9),"+100kb")) +
#   coord_fixed() +
#   theme_minimal() +
#   labs(x=NULL,y=NULL) +
#   theme(strip.text=element_text(face="bold",size=14),
#         axis.text=element_text(color="black"),
#         legend.position="top")
# ggsave(file.path(projdir,"combined_APA.png"),dpi=300,width=3,height=5)

# 25kb normed APA ---------------------------------------------------------

aged = read.table(file.path(projdir, "aged_merged_25kb_gw_normedAPA.txt"))
aged = apply(aged,2,function(x) gsub(",","",x))
aged = apply(aged,2,function(x) gsub("[","",x,fixed=T))
aged = apply(aged,2,function(x) gsub("]","",x,fixed=T))
aged = apply(aged,2,as.numeric)
# breaks = seq(0,7000,1000), legend_breaks = seq(0,7000,1000),
pheatmap(aged, color = colorRampPalette(brewer.pal(9, "Reds"))(500), legend_breaks = seq(0,6,1),
         cluster_rows=F, cluster_cols=F, show_rownames = F, show_colnames = F, border_color = NA)

# APA: 3.58439235269049
aged.melt = melt(as.matrix(aged), varnames=c("row","col")) %>% 
  mutate(col=factor(col),
         row=factor(row,levels=rev(unique(row))),
         Age="Aged")

ggplot(aged.melt, aes(x=col, y=row))+ 
  geom_tile(aes(fill=value), color=NA) +
  scale_fill_gradient(low="white", high="red", limits=c(0,5.5), na.value="black") +
  coord_fixed() +
  theme_void() +
  labs(x=NULL,y=NULL) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank())

# APA: 3.5230796050025144
young = read.table(file.path(projdir, "young_merged_25kb_gw_normedAPA.txt"))
young = apply(young,2,function(x) gsub(",","",x))
young = apply(young,2,function(x) gsub("[","",x,fixed=T))
young = apply(young,2,function(x) gsub("]","",x,fixed=T))
young = apply(young,2,as.numeric)

pheatmap(young, color = colorRampPalette(brewer.pal(9, "Reds"))(500), 
         cluster_rows=F, cluster_cols=F, show_rownames = F, show_colnames = F, border_color = NA)

young.melt = melt(as.matrix(young), varnames=c("row","col")) %>% 
  mutate(col=factor(col),
         row=factor(row,levels=rev(unique(row))),
         Age="Young")

ggplot(young.melt, aes(x=col, y=row))+ 
  geom_tile(aes(fill=value), color=NA) +
  geom_text(aes(x=3,20),label="3.58",size=7) +
  scale_fill_gradient(low="white", high="red", limits=c(0,5.5), na.value="black") +
  coord_fixed() +
  theme_void() +
  labs(x=NULL,y=NULL) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank())

combined = rbind(young.melt,aged.melt)
combined$Age = factor(combined$Age, levels=c("Young","Aged"))
label.df=data.frame(x=c(5,5,17,17), y=20, 
                    lab=c("APA=3.52","APA=3.58","max=5.18","max=5.47"), 
                    Age=factor(c("Young","Aged","Young","Aged"), levels=c("Young","Aged")))


ggplot(combined, aes(x=col, y=row)) +
  facet_grid(rows=vars(Age)) +
  geom_tile(aes(fill=value), color=NA) +
  geom_text(data=label.df, 
            aes(x=x,y=y,label=lab),
            size=4) +
  scale_fill_gradient(low="white", high="red", limits=c(0,5.5), na.value="black", name=NULL) +
  scale_x_discrete(labels=c("-100kb",rep("",9),"0",rep("",9),"+100kb")) +
  scale_y_discrete(labels=c("-100kb",rep("",9),"0",rep("",9),"+100kb")) +
  coord_fixed() +
  theme_minimal() +
  labs(x=NULL,y=NULL) +
  theme(strip.text=element_text(face="bold",size=14),
        axis.text=element_text(color="black"),
        legend.position="top")
ggsave(file.path(projdir,"combined_APA.png"),dpi=300,width=3,height=5)
