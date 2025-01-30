library(data.table); library(magrittr); library(purrr)

ah <- fread("addhealth/ah_common_3a.bim") %>%
  setnames(., colnames(.), c("CHR", "SNP", "V3", "BP", "A1", "A2"))
hapmap <- fread("hapmap3/hapmap_common_3a.bim") %>%
  setnames(., colnames(.), c("CHR", "SNP", "V3", "BP", "B.A1", "B.A2"))

joint <- merge(ah, hapmap, all = F)
joint[, `:=`(C.A1 = fcase(B.A1 == "A", "T", B.A1 == "T", "A", B.A1 == "C", "G", B.A1 == "G", "C"),
             C.A2 = fcase(B.A2 == "A", "T", B.A2 == "T", "A", B.A2 == "C", "G", B.A2 == "G", "C"))]

# compatible SNPs
compat <- joint[(A1==B.A1 & A2==B.A2) | (A1==C.A1 & A2==C.A2) | (A1==B.A2 & A2==B.A1) | (A1==C.A2 & A2==C.A1)]

# SNPs to flip
flip <- compat[!(A1==B.A1 & A2==B.A2) & !(A1==B.A2 & A2==B.A1)]

fwrite(compat[, .(SNP, A1)], "hapmap3/comp_snp_3a.a1", col.names = F, sep = "\t")
fwrite(flip[, .(SNP)], "hapmap3/flip_snp_3a.txt", col.names = F, sep = "\t")

