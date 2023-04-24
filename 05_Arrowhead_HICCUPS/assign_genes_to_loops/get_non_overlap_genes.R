library(dplyr)
library(tidyr)

young_domain=read.table("young.merged.loop_domain.knownGenes.bed")
aged_domain=read.table("aged.merged.loop_domain.knownGenes.bed")
young_diff_domain=read.table("young.merged.diff_loop_domain.knownGenes.bed")
aged_diff_domain=read.table("aged.merged.diff_loop_domain.knownGenes.bed")
tss=read.table("/nas/homes/benyang/HiC/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed")

write.table(anti_join(tss,young_domain), "young.merged.non_loop_domain.knownGenes.bed", quote=F, row.names=F, col.names=F)
write.table(anti_join(tss,aged_domain), "aged.merged.non_loop_domain.knownGenes.bed", quote=F, row.names=F, col.names=F)
write.table(anti_join(tss,young_diff_domain), "young.merged.diff_non_loop_domain.knownGenes.bed", quote=F, row.names=F, col.names=F)
write.table(anti_join(tss,aged_diff_domain), "aged.merged.diff_non_loop_domain.knownGenes.bed", quote=F, row.names=F, col.names=F)