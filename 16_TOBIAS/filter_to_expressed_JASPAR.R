library(dplyr)
library(tidyr)
library(biomaRt)
library(Hmisc)

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C"


# Get gene data -----------------------------------------------------------

tpm.data = read.table(file.path(projdir, "04_FANC", "compartmentExpression", "RNA_transformed_tpm_minus_sva_contribs.txt"))
gene.sym = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
gene.sym.tbl = data.frame(table(gene.sym)) # find duplicate gene symbols
duplicate.genes = gene.sym.tbl[gene.sym.tbl$Freq>1,"gene.sym"]

tpm.data = tpm.data[-which(gene.sym %in% duplicate.genes), 1:8] # isolate day0 data for young and aged and remove gene duplicates
rownames(tpm.data) = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
tpm.data = tpm.data[complete.cases(tpm.data), ] # remove NA rows
tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) all(x>=0))), ] # keep only genes with TPM>=0 in either young or aged


# Collect motif data ------------------------------------------------------

processFile = function(filepath) {
  con = file(filepath, "r")
  
  matrix.id = c()
  matrix.name = c()
  uniprot.id = c()
  
  while (TRUE) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    } else if (grepl("^AC MA", line)) {
      matrix.id = c(matrix.id, substr(line, 4, nchar(line)))
    } else if (grepl("^ID", line)) {
      matrix.name = c(matrix.name, gsub("ID ", "", line, fixed=T))
    } else if (grepl("^CC uniprot_ids", line)) {
      uniprot.id = c(uniprot.id, gsub("CC uniprot_ids:", "", line, fixed=T))
    }
  }
  
  close(con)
  return(data.frame(id = matrix.id, name = matrix.name, uniprot = uniprot.id))
}

motif.file = file.path(projdir, "17_TOBIAS", "JASPAR2022_CORE_vertebrates_non-redundant_pfms_transfac.txt")

# make data frame of JASPAR motif name and gene name
motif.data = processFile(motif.file)
motif.data$symbol = capitalize(tolower(motif.data$name))

# keep only genes with TPM>=1 in young or aged
tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) any(x>=1))), ] 

# excluded motifs likely came from Human or are for alternative cell fates (bone, neural)
# 1st pass -- get expressed motifs
motif.data.expressed = unique(motif.data$symbol[motif.data$symbol %in% rownames(tpm.data)])
# 2nd pass -- get expressed heterodimer motifs (include if at least one of the partners is expressed)
motif.data.heterodimer = grep("::", motif.data$symbol, value=T)
heterodimer.to.keep = sapply(strsplit(motif.data.heterodimer, split="::", fixed = T), function(x) all(capitalize(tolower(x)) %in% rownames(tpm.data)))
motif.data.heterodimer.expressed = motif.data.heterodimer[heterodimer.to.keep]

motif.data.expressed = c(motif.data.expressed, motif.data.heterodimer.expressed)
motif.id.expressed = motif.data$id[motif.data$symbol %in% motif.data.expressed]

length(motif.data.expressed)
length(motif.id.expressed)

# motif.uniprot = trimws(unlist(strsplit(motif.data$uniprot, split=';', fixed=T)))
# 
# ensembl <- useEnsembl(biomart = "genes")
# grep('GRCm38', listDatasets(ensembl)$version)
# searchDatasets(ensembl, "uniprot")
# ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
# 
# mapping = getBM(attributes = c('external_gene_name', 'mgi_symbol','uniprot_gn_id','uniprotswissprot'),
#                 filters = 'uniprotswissprot',
#                 values = unique(motif.uniprot), 
#                 mart = ensembl)
# 
# mapping = getBM(attributes = c('external_gene_name', 'mgi_symbol','uniprot_gn_id','uniprotswissprot'),
#                 filters = 'uniprot_gn_id',
#                 values = 'P05554',
#                 mart = ensembl)
mapping = getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                filters = 'mgi_symbol',
                values = 'Pax7',
                mart = ensembl)

library(EnsDb.Mmusculus.v79)
edb = EnsDb.Mmusculus.v79

ensembl_mapping = genes(edb, filter=SymbolFilter(unique(motif.data$symbol)), return.type="data.frame")
motif.data$gene_id = ensembl_mapping$gene_id[match(motif.data$symbol, ensembl_mapping$gene_name)]
motif.data$map_id = paste(motif.data$name, motif.data$id, sep="_")

motif2gene = motif.data[,c("map_id", "gene_id")]
write.table(motif2gene, file.path(projdir, "17_TOBIAS", "JASPAR_to_Ensembl.txt"), sep='\t', col.names=F, row.names=F, quote=F)

# Prepare JASPAR motifs ---------------------------------------------------

processJASPAR = function(filepath, outfilepath, motif.id.expressed) {
  con = file(filepath, "r")
  outfile = file(outfilepath, "w")
  
  while (TRUE) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    } else if (grepl("^>", line)) {
      id = trimws(unlist(strsplit(line, '\t', fixed=T))[1], which = "left", whitespace = ">")
      if (id %in% motif.id.expressed) {
        writeLines(line, outfile)
        writeLines(readLines(con, n=4), outfile)
      }
    } 
  }
  
  close(con)
  close(outfile)
}

file.create(file.path(projdir,'17_TOBIAS','expressed_JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt'))

processJASPAR(file.path(projdir,'17_TOBIAS','JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt'),
              file.path(projdir,'17_TOBIAS','expressed_JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt'),
              motif.id.expressed)
