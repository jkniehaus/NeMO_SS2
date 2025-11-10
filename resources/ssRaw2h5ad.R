library(optparse)
option_list = list(
  make_option(c("-n", "--ncore"), type="integer", default=NULL,
              help="number cores used", metavar="character"),
  make_option(c("-m", "--mem"), type="integer", default=NULL,
              help="memory use in gb", metavar="character"),
  make_option(c("-o", "--prefix"), type="character", default=NULL,
              help="output file name prefix", metavar="character"),
  make_option(c("-d", "--prefix2"), type="character", default=NULL,
              help="output file name prefix", metavar="character"),
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
rawdir <- paste0(opt$prefix,'/',opt$prefix,'GeneFull_ExonOverIntron/raw/')
rawdir2 <- paste0(opt$prefix2,'/',opt$prefix,'Gene/raw/')

print('Loading Libraries')
library(optparse)
library(Seurat)
library(Matrix)
library(rtracklayer)
library(reticulate)
library(anndata)
library(foreach)
library(doParallel)

cl<- makePSOCKcluster(opt$ncore-1)
clusterSetRNGStream(cl)
registerDoParallel(cl,cores=opt$ncore-1)

set.seed(6646428)
tic <- proc.time()
#compress files if necessary; assumes unix system with gzip installed
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

meta=read.csv('resources/U19.SS.metadat.allcol.csv',header=T,row.names=1)
cellnames=colnames(s)
cells=substr(cellnames,1,nchar(cellnames)-1)
metacells=rownames(meta) #cell names of AIBS cells
badcells=cells[which(cells %in% metacells == FALSE)]
metapaths=meta$fastq_path_list
extracells=cells #for safe keeping
print('parallel renaming')
if (opt$ncore > 2) {
  cl <- makePSOCKcluster(opt$ncore - 1)
  clusterSetRNGStream(cl)
  registerDoParallel(cl, cores = opt$ncore - 1)

  cells <- foreach(i = seq_along(extracells), .combine = c, .inorder = TRUE) %dopar% {
    j <- extracells[i]
    if (!(j %in% metacells) && any(grepl(j, metapaths, fixed = TRUE))) {
      metacells[grep(j, metapaths, fixed = TRUE)][1]
    } else {
      j
    }
  }

  message("Loop finished, number of results: ", length(cells))
  stopCluster(cl)

} else {
  message("Running serial version, skipping cluster setup")
  cells <- vapply(extracells, function(j) {
    if (!(j %in% metacells) && any(grepl(j, metapaths, fixed = TRUE))) {
      metacells[grep(j, metapaths, fixed = TRUE)][1]
    } else {
      j
    }
  }, FUN.VALUE = character(1))
  message("Loop finished, number of results: ", length(cells))
}
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
if (exists('s2')) {
  colnames(s2) <- uniqcells
  s2_assay <- CreateAssayObject(s2)
  s[['exons']] <- s2_assay
}
s <- subset(s,subset=percent.mt<15)
dir.create('rds')
saveRDS(s,paste0('rds/',opt$prefix,'.rds'))

gtf.file = "resources/Mus_musculus.GRCm39.106.ensembl.filtered.gtf"
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
write_h5ad(countAD, paste0('h5ads/',opt$prefix,'.h5ad')) # Write it out as h5ad
