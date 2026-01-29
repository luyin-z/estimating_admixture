library(tidyverse); library(bigsnpr); library(data.table); library(rhdf5)

# read in the AH-Ref data for ADMIXTURE supervised 4
ref <- bed('supervised4.bed') # N = 11581

# identify AH respondents and references
ah_fam <- fread('addhealth/ahgwas.fam')
ind_row <- which(!ref$fam$family.ID %in% ah_fam$V1 & !ref$fam$sample.ID %in% ah_fam$V2) # N = 2254


# perform PCA on the reference set
obj_svd <- bed_autoSVD(ref, k = 4, ind.row = ind_row, thr.r2 = NA, ncores = nb_cores())

# identify and remove outliers
prob <- bigutilsr::prob_dist(obj_svd$u, ncores = nb_cores())
S <- prob$dist.self / sqrt(prob$dist.nn)

ind_row2 <- which(S < 0.5) # N = 2254
(outlier <- ind_row[-ind_row2]) # none
(outlier_id <- ref$fam[outlier,]$sample.ID) # none
ind_col <- attr(obj_svd, 'subset') # N = 275766

# rerun SVD
obj_svd2 <- bed_autoSVD(ref, ind.row = ind_row[ind_row2], ind.col = ind_col, thr.r2 = NA,
                        k = 4, ncores = nb_cores())


# project AH onto the PCA space
system('plink1.9 --bfile addhealth/ahgwas --extract ahref_all.bim --a1-allele addhealth/ahgwas.bim 5 2 --make-bed --out temp')
ah <- bed('temp.bed')
system('rm temp.*')

ah_project <- bed_projectSelfPCA(obj_svd2, ah, ind.col = attr(obj_svd2, 'subset'),
                                 ind.row = rows_along(ah), ncores = nb_cores())

ah_PCs <- as.data.frame(ah_project$OADP_proj) # N = 9974
colnames(ah_PCs) <- paste0('PC', 1:4)
ah_PCs$FID <- ah$fam$family.ID
ah_PCs$IID <- ah$fam$sample.ID

saveRDS(ah_PCs, 'ah_projected_PC4.rds')

