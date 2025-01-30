library(here); library(tidyverse); library(bigsnpr); library(data.table)
library(rhdf5)

setwd(str_remove(here(), '/scripts'))

# read in the AH-HapMap3-HGDP dataset that is used for supervised ADMIXTURE 5
ref <- bed('ah_hapmap3_hgdp_s5_nosex.bed') # n = 9590

# distinguish between AH respondents and the references
ah_fam <- fread('addhealth/ah_unrel_s2_3a.fam')
ind_row <- which(!ref$fam$family.ID %in% ah_fam$V1 & !ref$fam$sample.ID %in% ah_fam$V2) # n = 424 = 9590 - 9166

# perform PCA that creates the first 5 PCs
obj_svd <- bed_autoSVD(ref, k = 5, ind.row = ind_row, thr.r2 = NA, ncores = nb_cores())
# Skipping clumping.
# Discarding 80 variants with MAC < 10.
# 
# Iteration 1:
#   Computing SVD..
# 0 outlier variant detected..
# 
# Converged!

# identify and remove any outlier(s)
prob <- bigutilsr::prob_dist(obj_svd$u, ncores = nb_cores())
S <- prob$dist.self / sqrt(prob$dist.nn)
ind_row2 <- which(S < 0.5) # n = 423
(outlier <- ind_row[-ind_row2]) # index: 56
(outlier_id <- ref$fam[outlier,]$sample.ID) # HGDP00621

# rerun PCA
ind_col <- attr(obj_svd, 'subset') # n = 210613 SNPs
obj_svd2 <- bed_autoSVD(ref, ind.row = ind_row[ind_row2], ind.col = ind_col, thr.r2 = NA,
                        k = 5, ncores = nb_cores())
# Skipping clumping.
# Discarding 1 variants with MAC < 10.
# 
# Iteration 1:
#   Computing SVD..
# 0 outlier variants detected..
# 
# Converged!

# project AH respondents onto the reference PC space
ah <- bed('addhealth/ahpruned_s2_3a.bed')
ah_project <- bed_projectSelfPCA(obj_svd2, ah, ind.col = attr(obj_svd2, 'subset'),
                                 ind.row = rows_along(ah), ncores = nb_cores())
ah_PCs <- as.data.frame(ah_project$OADP_proj) # n = 9974
colnames(ah_PCs) <- paste0('PC', 1:5)
ah_PCs$FID <- ah$fam$family.ID
ah_PCs$IID <- ah$fam$sample.ID

# save
saveRDS(obj_svd2, 'ref_5PC_bigsnpr_aligned_no_outliers.rds')
saveRDS(ah_project, 'ah_PCprojection.rds')
saveRDS(ah_PCs, 'ah_projected_PC5.rds')

