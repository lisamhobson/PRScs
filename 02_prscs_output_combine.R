################################################################
### combine output of PRS-CS weight when run per chromosome  ###
################################################################  

# read in parameter file, containing "output_path"
source("")

# from 01_gwascatalog_to_prscs.R script 
id <- "" # gwas file name (without .file extension)
p2_res <- paste0(output_path) # weights output path

# import each chromosome weight file 
dat <- lapply(Sys.glob(paste0(p2_res, "/", id, "/_pst_eff_a1_b0.5_phi1e-02_chr*.txt")), read.table)

# compile into dataframe
dat_total <- dat[[1]]
for(i in 2:length(dat)) 
  dat_total <- rbind.data.frame(dat_total, dat[[i]])

# write to file
write.table(dat_total, file=paste0(p2_res, id, "/", id, "_chr1-22.txt"), quote=F, row.names = F, col.names = F)
