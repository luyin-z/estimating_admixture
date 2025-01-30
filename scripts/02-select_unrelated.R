library(readr); library(tidyr); library(purrr)

dat <- list(read_table('addhealth/ahgwas.fam', col_names = c('FID', 'IID', 'F', 'M', 'S', 'P'))) # N = 9974
dat <- append(dat, read_table('hapmap3/relation_4a.txt', col_names = c('FID', 'IID', 'F', 'M', 'S', 'P', 'Pop')) %>%
  group_by(Pop) %>%
  group_split())

set.seed(218)

unrelated <- dat %>%
  map(~{
    .x <- .x %>%
      group_by(FID) %>%
      slice_sample(n = 1)
  })

write.table(unrelated[[1]], 'addhealth/ahgwas_unrelated.fam', row.names = F, col.names = F)
write.table(unrelated[[c(2,3,5)]], 'hapmap3/hapmap_unrelated_3a.txt', quote = F, row.names = F, col.names = F)
# write.table(unrelated[[2:5]], 'hapmap3/hapmap_unrelated_4a.txt', quote = F, row.names = F, col.names = F)

