library(readr); library(tidyr); library(dplyr); library(purrr)

rel_hapmap <- read_table("hapmap3/relation_4a.txt", 
                    col_names = c("FID","IID","F","M","S","P","A"),
                    col_types = "ccccnnc") %>% select(-P)
hgdp_3a <- read_table("hgdp/hgdp_3a.psam") %>% select(IID = `#IID`, A = region)
ah_hapmap <- read_table("ah_hapmap3_s3_nosex.fam", col_names = c("FID","IID","F","M","S","P"),
                        col_types = "ccccnn")
ah_hapmap_hgdp <- list()
for (k in paste0("s", 4:6)) {
  ah_hapmap_hgdp[[k]] <- read_table(paste0("ah_hapmap3_hgdp_",k,"_nosex.fam"), 
                                    col_names = c("FID","IID","F","M","S","P"),
                                    col_types = "ccccnn") %>% select(-P) %>%
    left_join(hgdp_3a)
}

merge <- ah_hapmap_hgdp %>% map(~{
  .x <- .x %>%
    left_join(rel_hapmap %>% rename(Ah = A)) %>%
    mutate(A = case_when(
      A=="AMERICA" ~ "NA",
      A=="MIDDLE_EAST" ~ "ME",
      A=="OCEANIA" ~ "OC",
      !is.na(Ah) ~ Ah,
      T ~ "-"
    )) %>%
    select(-Ah)
})

merge2 <- ah_hapmap %>%
  left_join(rel_hapmap) %>%
  mutate(A = ifelse(is.na(A), "-", A))
  
lapply(merge, function(d) {
  table(d$A, useNA = "ifany")
  #    -  CEU  CHB  NA  YRI
  # 9166   36  137  34   83
  #    -  CEU  CHB   ME   NA  YRI
  # 9166   36  137  134   34   83
  #    -  CEU  CHB   ME   NA   OC  YRI
  # 9166   36  137  134   34   28   83
})
  
table(merge2$A, useNA = "ifany")
#    -  CEU  CHB  YRI 
# 9166   36  137   83

lapply(paste0("s", 4:6), function(k) {
  write.table(merge[[k]]$A, paste0("ah_hapmap3_hgdp_", k, "_nosex.pop"), 
              quote = F, row.names = F, col.names = F)
})
write.table(merge2[merge2$A!="CHB",]$A, "ah_hapmap3_s2_nosex.pop", quote = F, row.names = F, col.names = F)
write.table(merge2$A, "ah_hapmap3_s3_nosex.pop", quote = F, row.names = F, col.names = F)

