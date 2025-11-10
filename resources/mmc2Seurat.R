library(optparse)
option_list = list(
  make_option(c("-m", "--mem"), type="integer", default=1,
              help="memory use in gb", metavar="character"),
  make_option(c("-p", "--prefix"), type = 'character', default = '',
              help='prefix for rds file', metavar='character')
);
#library(remotes)
#remotes::install_github('jkniehaus/seuratHelpR')
#library(seuratHelpR)
source('/proj/gs25/users/Jesse/scripts/DSSeurat_functions_v5working.R')
library(DoubletFinder)
library(Seurat)
library(future)
opt_parser = OptionParser(option_list=option_list)
opt=parse_args(opt_parser)
options(future.globals.maxSize = opt$mem * 1024 ^ 3)
options(warn=1)
s <- readRDS(paste0('rds/',opt$prefix,'.rds'))
h <- paste0('rds/',opt$prefix,'.csv')
comment_lines <- sum(readLines(h, n = 10) %in% grep ("^#", readLines(h), value = TRUE))
mp <- read.csv(h, skip = comment_lines, header=TRUE)
mp <- mp[mp$cell_id %in% colnames(s),]
s$Class <- mp$class_name
s$ClassScore <- mp$class_bootstrapping_probability
s$Subclass <- mp$subclass_name
s$SubclassScore <- mp$subclass_bootstrapping_probability
s$Supertype <- mp$supertype_name
s$SupertypeScore <- mp$supertype_bootstrapping_probability
s$Cluster <- mp$cluster_name
s$ClusterScore <- mp$cluster_bootstrapping_probability
Idents(s) <- s$Class

neuro <- c('Glut','GABA','Dopa','Sero','Chol','Hist','Nora','Glyc')
basic <- sapply(mp$class_name, function(x) {
# Extract the last four characters of each string and see if it's key 'neuro'
last_four <- substr(x, nchar(x) - 3, nchar(x))
ifelse(length(grep(last_four, neuro)) > 0, 'Neuron', 'Glia')
})
names(basic) <- NULL
s$BasicCellType <- basic
n <- subset(s,subset=BasicCellType=='Neuron')
s$Neurotransmitter <- substr(mp$class_name,nchar(mp$class_name)-3,nchar(mp$class_name))
saveRDS(s,paste0('rds/',opt$prefix,'.rds'))
