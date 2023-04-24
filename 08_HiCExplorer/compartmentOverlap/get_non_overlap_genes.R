library(dplyr)
library(tidyr)

# tss=read.table("/nas/homes/benyang/HiC/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed")

tasks=grep(pattern="non",list.files(pattern="promoter.knownGenes.bed",path="promoterOverlap"),fixed=T,invert=T,value=T)
for(taskname in tasks) {
    fname=gsub(pattern=".promoter.knownGenes.bed","",taskname)
    age_comp=paste0(unlist(strsplit(fname,split="\\."))[1:2],collapse=".")
    comp=unlist(strsplit(fname,split="\\."))[2]
    TAD_group=unlist(strsplit(fname,split="\\."))[3]

    tss = if(comp %in% c("A","B")) {
        read.table(file.path("/nas/homes/benyang/HiC/04_FANC/compartmentExpression", paste0(age_comp,".1kb.tss.bed")))
    } else {
        read.table(file.path("/nas/homes/benyang/HiC/04_FANC/compartmentExpression", paste0(comp,".1kb.tss.bed")))
    }

    data=read.table(file.path("promoterOverlap",taskname))
    data_genes=data$V4

    write.table(tss[!(tss$V4 %in% data_genes),], 
    file.path("promoterOverlap",paste0(age_comp,".non_",TAD_group,".promoter.knownGenes.bed")), 
    quote=F, row.names=F, col.names=F)
}
