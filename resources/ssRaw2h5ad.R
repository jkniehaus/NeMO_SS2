library(optparse)
option_list = list(
  make_option(c("-n", "--ncore"), type="integer", default=NULL,
              help="number cores used", metavar="character"),
  make_option(c("-m", "--mem"), type="integer", default=NULL,
              help="memory use in gb", metavar="character"),
  make_option(c("-o", "--prefix"), type="character", default=NULL,
              help="output file name prefix", metavar="character"),
  make_option(c("-d", "--prefix2"), type="character", default=NULL,
              help="output file name prefix", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
rawdir <- paste0(opt$prefix,'/',opt$prefix,'GeneFull_ExonOverIntron/raw/')
rawdir2 <- paste0(opt$prefix2,'/',opt$prefix,'GeneFull_ExonOverIntron/raw/')

print('Loading Libraries')
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(findPC)
library(future)
library(dplyr)
library(HGNChelper)
library(openxlsx)
library('Matrix')
library('reshape2')
library('sctransform')
library('knitr')
library('ggrepel')
library('patchwork')
library(data.table)
library(rhdf5)
library(reticulate)
library(rtracklayer)
library(anndata)
use_python('/nas/longleaf/rhel9/apps/python/3.12.2/bin/python')
py_config()

list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "kableExtra"
  )

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
      )
    )
}

cl<- makePSOCKcluster(opt$ncore-1)
clusterSetRNGStream(cl)
registerDoParallel(cl,cores=opt$ncore-1)


knit_hooks$set(optipng = hook_optipng)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = FALSE,
  warning = FALSE,
  digits = 2,
  tidy = TRUE,
  tidy.opts = list(width.cutoff=80),
  optipng = '-o 5 -strip all -quiet',
  fig.width=6.5, fig.height=2.5, dpi=100, out.width = '80%'
)
options(dplyr.summarise.inform = FALSE)
options(tibble.width = Inf)
old_theme <- theme_set(theme_bw(base_size=10))
set.seed(6646428)
tic <- proc.time()


#compress files if necessary
if (length(grep('*.gz',list.files(rawdir)))==0) {
  use <- paste0('gzip ',rawdir,'*')
  system(use)
}

s=Read10X(rawdir)
if (dir.exists(paste0(opt$prefix2))) {
    if (length(grep('*.gz',list.files(rawdir2)))==0) {
        use <- paste0('gzip ',rawdir2,'*')
        system(use)
    }
    s2 <- Read10X(rawdir2)
    s2 <- s2[,colnames(s2) %in% colnames(s)]
    s <- s[,colnames(s) %in% colnames(s2)]
}


meta=read.csv('/proj/gs25/projects/Allen_SmartSeq/U19.SS.metadat.allcol.csv',header=T,row.names=1)
cellnames=colnames(s)
cells=substr(cellnames,1,nchar(cellnames)-1)
metacells=rownames(meta) #cell names of AIBS cells
badcells=cells[which(cells %in% metacells == FALSE)]
metapaths=meta$fastq_path_list
extracells=cells #for safe keeping
print('parallel renaming')
system.time(
cells <- foreach(i = 1:length(extracells), .combine=c, .inorder = TRUE) %dopar%{
        j=extracells[i]
        if (j %in% metacells == FALSE && length(grep(j,metapaths))>0) {
                metacells[grep(j,metapaths)]
        } else {
        extracells[i]
        }
}
)
dups=cells[duplicated(cells)]
dupnum=which(duplicated(cells))
if (length(dups)>0) {
print(dim(s))
m3=s[grep(dups[1],cells),] #check if duplicates = same cell
print(head(m3,20))
print('after dup cell collapse gene matrix')
s <- s[,-c(dupnum)]
if (exists('s2')) {
    s2 <- s2[,-c(dupnum)]
}
print(dim(s))
}
uniqcells=unique(cells)
colnames(s)=uniqcells
metasub=meta[colnames(s),]
s=CreateSeuratObject(s,meta.data=metasub)
s@meta.data$nCount=colSums(s[['RNA']]$counts)
s@meta.data$nGeneFeatures=colSums(s[['RNA']]$counts > 0)
nfeatcell=s@meta.data$nGeneFeatures
s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^mt-")
s=subset(s,subset=percent.mt<15)
if (exists('s2')) {
   colnames(s2) <- uniqcells
   s2 <- CreateAssayObject(s2)
   s[['exons']] <- s2
}
dir.create('rds')
saveRDS(s,paste0('rds/',opt$prefix,'.rds'))

gtf.file = "/proj/gs25/users/Jesse/references/annotations/Mus_musculus.GRCm39.106.ensembl.filtered.gtf"
gtf.gr = rtracklayer::import(gtf.file) # creates a GRanges object
gtf.df = as.data.frame(gtf.gr)
genes = unique(gtf.df[ ,c("gene_id","gene_name")])
#rename empty gene names = gene_id
for (i in 1:length(genes$gene_name)){
  if (genes$gene_name[i] == '' | is.na(genes$gene_name[i]) == TRUE) {
    genes$gene_name[i] <- genes$gene_id[i]
  }
}

mtx <- as(s[['RNA']]$counts, Class='dgCMatrix')
counts  <- t(mtx)
counts <- counts[,colnames(counts) %in% genes$gene_name]
feats   <- colnames(counts)
genes <- genes[genes$gene_name %in% feats,]
genes.ordered <- genes[match(feats,genes$gene_name),]
var.use <- genes.ordered$gene_id
samples <- rownames(counts)
sparse_counts <- as(counts, "dgCMatrix")  # Convert to sparse matrix, if not already
print('generating anndata')
countAD <- AnnData(X   = sparse_counts,   # Create the anndata object
                   var = data.frame(genes=var.use,row.names=var.use),
                   obs = data.frame(samples=samples,row.names=samples))
print('writing h5ad')
dir.create('h5ads')
write_h5ad(countAD, paste0('h5ads/',opt$prefix,'.h5ad')) # Write it out as h5ad
