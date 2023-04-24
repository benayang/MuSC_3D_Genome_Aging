library(dplyr)
library(tidyr)

tasks=grep(pattern=".tab",list.files(pattern="ID","get_flank_coverage"),value=T)

for(taskname in tasks){
    fname=gsub(pattern=".tab",replacement="",taskname)
    print(fname)
    tmp=read.delim2(file.path("get_flank_coverage",taskname), col.names=c("chr","start","end","H4K20me1","ATAC","H3K27me3","H3K4me3","ID"))
    meanCoverage=tmp %>% pivot_longer(cols=c(H4K20me1,ATAC,H3K27me3,H3K4me3), names_to="Marks", values_to="Coverage") %>% 
    mutate(Coverage=as.numeric(Coverage)) %>% group_by(ID,Marks) %>% summarise(meanCoverage=mean(Coverage)) %>% pivot_wider(names_from=Marks,values_from=meanCoverage)
    write.table(meanCoverage,file.path("get_averaged_flank",paste0(fname,"_avg.tsv")),quote=F,row.names=F,col.names=T,sep="\t")
}