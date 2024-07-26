#################################################
### GWAS summary stats to PRS-CS weights      ###
#################################################  

# clone repository and read requirements: https://github.com/lisamhobson/PRScs/

##### packages #####
list_of_packages <- c("tidyverse", "glue")
for (i in list_of_packages){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}
#####

# read in parameter file, containing "data_path", "gwas_path", "output_path"
source("")

# gwas info
id = "" # gwas file name (without .file extension)
n = 0 # n of gwas

# prs-cs parameters 
phi = 1e-2 # shrinkage, 1e-2 for highly polygenic
seed = 1 

###############################################################################################
# set paths
chr <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # chromosomes to run (run as array)
p2_prscs <- paste0(data_path, "PRScs/PRScs.py") # PRScs package from github repo
p2_gwas <- paste0(gwas_path)  # gwas path
p2_ref <- paste0(data_path, "ldblk_ukbb_eur") # From the PRScs github repo - specific to ancestry
p2_bim <- paste0(data_path, "") # bim file 
p2_res <- paste0(output_path) # weights output path

# create directory for output
p2_res_id <- paste0(p2_res, id, "/")
if (!file.exists(p2_res_id)) {
  dir.create(p2_res_id, showWarnings = T)
  next
}

# create output file per chromosome
paste0(p2_res_id, "_pst_eff_a1_b0.5_phi1e-2_chr", chr, ".txt")

# read in gwas file 
p2_file <- paste0(p2_gwas, id, ".tsv")
gwas <- read_tsv(p2_file, col_names = T, show_col_types = F)

# check column names of file, edit as appropriate
colnames(gwas)

# rename columns for prs-cs
colnames(gwas)[which(names(gwas) == "variant_id")] <- "SNP"
colnames(gwas)[which(names(gwas) == "p_value")] <- "P"
colnames(gwas)[which(names(gwas) == "chromosome")] <- "CHR"
colnames(gwas)[which(names(gwas) == "base_pair_location")] <- "LOC"
colnames(gwas)[which(names(gwas) == "effect_allele")] <- "A1"
colnames(gwas)[which(names(gwas) == "other_allele")] <- "A2"
colnames(gwas)[which(names(gwas) == "odds_ratio")] <- "OR"
colnames(gwas)[which(names(gwas) == "beta")] <- "BETA"
colnames(gwas)[which(names(gwas) == "standard_error")] <- "SE"

# select column to use to generate weights (BETA + SE / OR + SE / BETA + P / OR + P)
gwas2 <- select(gwas, c("SNP", "A1","A2", "OR", "SE"))

# temp file
temp <- paste0(tempfile(), ".txt")
write.table(gwas2, temp, sep = "\t", row.names = F, quote = F)

# set prs-cs parameters using above variables
cmd <- glue::glue("
  python {p2_prscs} --ref_dir={p2_ref} --bim_prefix={p2_bim} --sst_file={temp} --n_gwas={n} --phi={phi} --seed={seed} --out_dir={p2_res_id} --chrom={chr}
  ")

# run prs-cs
system(cmd)


