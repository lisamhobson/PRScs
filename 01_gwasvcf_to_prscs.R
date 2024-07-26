##############################################
### GWAS vcf sum stats to PRS-CS weights   ###
##############################################  

list_of_packages <- c("tidyverse", "glue")
for (i in list_of_packages){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

source("") # parameter file containing paths for file locations

set_plink(plink_path) 
set_bcftools()

# set paths
chr <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # chromosomes to run
p2_prscs <- paste0(git_prs, "PRScs2/PRScs.py") # PRScs package from github repo
p2_gwas <- paste0(openGWAS)  # gwas path
p2_ref <- paste0(data_path, "ldblk_ukbb_eur") # From the PRScs github repo - specific to ancestry
p2_bim <- paste0(data_path, "UKBB_10K") # bim file 
p2_res <- paste0(output_path) # weights output path

id <- "" # gwas id
n <-  # gwas n
phi=1e-2 # shrinkage, 1e-2 for highly polygenic
seed=1 

p2_res_id <- paste0(p2_res, id, "/")
dir.create(p2_res_id, showWarnings = F)

if (file.exists(paste0(p2_res_id, "_pst_eff_a1_b0.5_phi1e-02_chr", chr, ".txt"))) {
  next  
}

p2_file <- paste0(p2_gwas, id, "/", id, ".vcf.gz")
if (!file.exists(p2_file)) {
  next
}

vcf <- VariantAnnotation::readVcf(p2_file) %>%
  gwasvcf::vcf_to_tibble()

if (nrow(vcf) < 1) {
  stop("Number of rows < 1")
}

vcf$LP <- 10^-vcf$LP
# Cut down and save temp file
dat <- vcf[c("rsid", "ALT", "REF", "ES", "LP")]
names(dat) <- c("SNP", "A1", "A2", "BETA", "P")

temp <- paste0(tempfile(), ".txt")
write.table(dat, temp, sep = "\t", row.names = F, quote = F)

cmd <- glue::glue("
  python {p2_prscs} --ref_dir={p2_ref} --bim_prefix={p2_bim} --sst_file={temp} --n_gwas={n} --phi={phi} --seed={seed} --out_dir={p2_res_id} --chrom={chr}
  ")

system(cmd)

rm(vcf)
