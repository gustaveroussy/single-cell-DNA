set.seed(1)
library("dplyr")
library("dsb")
library("Seurat")
library("here")
library("ggplot2")
library("optparse")
library("Matrix")
library("textshape")


# Get parameters
option_list <- list(
  ### Project
  make_option("--input.file", help="Reads Counts proteine in tsv (DSB read tsv file)"),
  make_option("--output.file", help="Output file path"),
  make_option("--isotype.control",help="list of isotype control")
)

parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
input.file <- args$options$input.file
output.file <- args$options$output.file
isotype_control_name_vec <- args$option$isotype.control

isotype_control_name_vec <- as.vector(unlist(strsplit(isotype_control_name_vec,",")))

counts = read.delim(file = input.file, sep = "\t",header = T)

# make dataframe
counts = as.data.frame(counts)
counts[is.na(counts)] =  0 

prot = counts %>% 
  textshape::column_to_rownames("read_counts") %>% t() %>% as.data.frame()
  
# calculate library size of droplets to make rough thresholds for cell containing and ambient droplets 
prot_size = log10(Matrix::colSums(prot, na.rm = TRUE))
md = as.data.frame(prot_size)

md$bc = colnames(prot)
#png(hist.path)
#hist(md$prot_size, breaks = 100)
#dev.off()

background_drops = md[md$prot_size < quantile(md$prot_size,0.1)& md$prot_size > quantile(md$prot_size,0.01), ]$bc
negative_mtx_rawprot = prot[ , background_drops] %>%  as.matrix()

# define a vector of cell-containing droplet barcodes based on protein library size 
positive_cells = md[md$prot_size > quantile(md$prot_size,0.1)+0.01, ]$bc
cells_mtx_rawprot = prot[ , positive_cells] %>% as.matrix()


dsb_norm_prot = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot,
  empty_drop_matrix = negative_mtx_rawprot,
  denoise.counts = TRUE,
  use.isotype.control = TRUE,
  isotype.control.name.vec = isotype_control_name_vec
)

df_to_save <- t(dsb_norm_prot)
write.csv(dsb_norm_prot,output.file)