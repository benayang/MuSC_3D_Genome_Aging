library(dplyr)

projdir = '/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS'

young_merged = read.table(file.path(projdir, 'young_merged_hiccups', 'merged_loops_domains.bed'))
young_diff = read.table(file.path(projdir, 'hiccups_diff', 'young_diff_loop_domains.bed'))
aged_merged = read.table(file.path(projdir, 'aged_merged_hiccups', 'merged_loops_domains.bed'))
aged_diff = read.table(file.path(projdir, 'hiccups_diff', 'aged_diff_loop_domains.bed'))

young_non_diff = anti_join(young_merged, young_diff, by=c("V1","V2","V3"))
aged_non_diff = anti_join(aged_merged, aged_diff, by=c("V1","V2","V3"))

write.table(young_non_diff, file.path(projdir, 'hiccups_diff', 'young_non_diff_loop_domains.bed'), sep='\t', quote=F, row.names=F, col.names=F)
write.table(aged_non_diff, file.path(projdir, 'hiccups_diff', 'aged_non_diff_loop_domains.bed'), sep='\t', quote=F, row.names=F, col.names=F)