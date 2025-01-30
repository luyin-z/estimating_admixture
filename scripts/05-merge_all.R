# 0. Init ----------------------------------------------------------------------
library(here); library(readr); library(haven); library(tidyverse)

setwd(str_remove(here(), '/scripts'))



# 1. Import data ---------------------------------------------------------------
# 1.1 ADMIXTURE results
# genetic similarity proportions
us <- list()
for (k in 2:7) {
  unrel <- paste0('ah_us_nosex.',k,'.Q')
  rel <- paste0('ah_rest_nosex.',k,'.Q')
  us[[paste0('us',k)]] <- read_table(unrel, col_names = FALSE) %>%
    rbind(read_table(rel, col_names = FALSE))
}

s <- list()
for (k in 2:6) {
  if (k <= 3) {
    unrel <- paste0('ah_hapmap3_s',k,'_nosex.',k,'.Q')
  } else {
    unrel <- paste0('ah_hapmap3_hgdp_s',k,'_nosex.',k,'.Q')
  }
  rel <- paste0('ah_rest_s',k,'_nosex.',k,'.Q')
  s[[paste0('s',k)]] <- read_table(unrel, col_names = FALSE) %>%
    rbind(read_table(rel, col_names = FALSE))
}

# SEs
se <- list()
for (k in 2:6) {
  unrel <- paste0('ah_us_nosex.',k,'.Q_se')
  rel <- paste0('ah_rest_nosex.',k,'.Q_se')
  se[[paste0('us',k)]] <- read_table(unrel, col_names = FALSE) %>%
    rbind(read_table(rel, col_names = FALSE))
}
for (k in 2:6) {
  if (k <= 3) {
    unrel <- paste0('ah_hapmap3_s',k,'_nosex.',k,'.Q_se')
  } else {
    unrel <- paste0('ah_hapmap3_hgdp_s',k,'_nosex.',k,'.Q_se')
  }
  rel <- paste0('ah_rest_s',k,'_nosex.',k,'.Q_se')
  se[[paste0('s',k)]] <- read_table(unrel, col_names = FALSE) %>%
    rbind(read_table(rel, col_names = FALSE))
}


# 1.2 .fam files
us_fam <- read_table('addhealth/ah_us_nosex.fam', col_types = 'ccnnnn', 
                     col_names = c('FID','IID','F','M','S','P')) %>%
  rbind(read_table('addhealth/ah_rest_nosex.fam', col_types = 'ccnnnn',
                   col_names = c('FID','IID','F','M','S','P')))

s_fam <- list()
for (k in 2:6) {
  if (k<=3) {
    unrel <- paste0('ah_hapmap3_s',k,'_nosex.fam')
  } else {
    unrel <- paste0('ah_hapmap3_hgdp_s',k,'_nosex.fam')
  }
  rel <- paste0('addhealth/ah_rest_s',k,'_nosex.fam')
  s_fam[[paste0('s',k)]] <- read_table(unrel, col_types = 'ccnnnn',
                                       col_names = c('FID','IID','F','M','S','P')) %>%
    rbind(read_table(rel, col_types = 'ccnnnn', col_names = c('FID','IID','F','M','S','P')))
}


# 1.3 linkage data
link <- read_dta('addhealth/GID_link.dta')


# 1.4 survey data
ah <- list()
for (i in 1:5) {
  ah[[paste0('w',i)]] <- read_dta(paste0('addhealth/wave', i, '.dta'))
}
ah$wt4 <- read_dta('addhealth/weights4.dta')
ah$ses <- read_dta('addhealth/conses.dta')


# 1.5 geographic ancestry
fam_origin <- readRDS('addhealth/family_origin_dum.rds')



# 2. Prepare data --------------------------------------------------------------
# 2.1 combine ADMIXTURE results with corresponding .fam files
us_dat <- us %>%
  map(~{
    .x <- .x %>%
      cbind(us_fam)
  }) %>%
  set_names(paste0('us',2:6))
s_dat <- map2(s, s_fam, function(x,y) cbind(x,y))
dat <- c(us_dat, s_dat)


# 2.2 key variables from survey
ppheno <- ah$w4 %>% # W4 variables
  transmute(
    aid,
    # demographics
    sex = bio_sex4 %>% factor(labels = c('Male', 'Female')),
    byear = h4od1y,
    bmonth = h4od1m,
    age4 = (iyear4-byear)*12 + imonth4-bmonth,
    edu = case_when(
      h4ed2 %in% c(1,2) ~ 'Less than high school',
      h4ed2==3 ~ 'High school',
      h4ed2 %in% 4:6 ~ 'Some college',
      h4ed2>=7 & h4ed2<=13 ~ "Bachelor's degree or above",
      T ~ NA
    ) %>% factor(levels = c('Less than high school', 'High school', 'Some college',
                            "Bachelor's degree or above")),
    edu2 = case_when(
      h4ed2<=6 ~ edu,
      h4ed2 %in% c(7:8,12) ~ 'College degree',
      h4ed2 %in% c(9:11,13) ~ 'Graduate degree',
      T ~ NA
    ) %>% factor(levels = c('Less than high school', 'High school', 'Some college',
                            'College degree', 'Graduate degree')),
    # attractiveness & racial classification
    h4ir1, h4ir4,
    # interviewer ID
    intid4
  ) %>%
  # W3 variables
  full_join(
    ah$w3 %>% 
      mutate(age3 = (iyear3-h3od1y)*12 + imonth3-h3od1m,
             across(c('h3od2', starts_with('h3od4')), function(x) x = ifelse(x>1, NA, x)),
             racec = case_when(
               h3od2==1 ~ 'Hispanic (of any race)',
               h3od2==0 & h3od4a+h3od4b+h3od4c+h3od4d>1 ~ 'Non-Hispanic Multiple Races',
               h3od2==0 & h3od4a==1 & h3od4b+h3od4c+h3od4d==0 ~ 'Non-Hispanic White', 
               h3od2==0 & h3od4b==1 & h3od4a+h3od4c+h3od4d==0 ~ 'Non-Hispanic Black', 
               h3od2==0 & h3od4c==1 & h3od4a+h3od4b+h3od4d==0 ~ 'Non-Hispanic Native American', 
               h3od2==0 & h3od4d==1 & h3od4a+h3od4b+h3od4c==0 ~ 'Non-Hispanic Asian & Pacific Islander', 
               T ~ NA
             ) %>% factor(levels = c('Non-Hispanic White', 'Non-Hispanic Black', 
                                     'Non-Hispanic Native American',
                                     'Non-Hispanic Asian & Pacific Islander', 
                                     'Hispanic (of any race)', 'Non-Hispanic Multiple Races')),
             skin = h3ir17 %>% factor(labels = c('Black', 'Dark Brown', 'Medium Brown', 'Light Brown', 'White'))) %>%
      select(aid, intid3, age3, racec, skin, h3od2, starts_with('h3od4'), h3od6, h3ir1, h3ir4)) %>%
  # W1 and W2 variables
  full_join(ah$w1 %>% select(aid, intid1 = intid, iyear, imonth, h1gi1y, h1gi1m, 
                             h1ir1, h1ir4, h1gi4, h1gi8, h1gi9, starts_with('h1gi6'))) %>%
  full_join(ah$w2 %>% select(aid, intid2, iyear2, imonth2, h2gi1y, h2gi1m, h2ir1)) %>%
  mutate(intid1 = as.numeric(intid1), intid2 = as.numeric(intid2), 
         across(c(h1gi1y, h1gi1m), ~ ifelse(.>=96, NA, .)),
         across(c(h2gi1y, h2gi1m), ~ ifelse(.>=98, NA, .)),
         age1 = (iyear-h1gi1y)*12 + imonth-h1gi1m,
         age2 = (iyear2-h2gi1y)*12 + imonth2-h2gi1m,
         across(matches('h[1-4]ir1'), ~ ifelse(.>5, NA, .))) %>%
  mutate(
    # interviewer's racial classification of respondents
    irace = factor(ifelse(h1gi9>=6, NA, h1gi9), 
                   labels = c('White', 'Black', 'Native American', 'Asian & Pacific Islander', 'Other')),
    irace3 = factor(h3ir4, labels = c('White', 'Black', 'Native American', 'Asian & Pacific Islander')),
    irace4 = factor(h4ir4, labels = c('White', 'Black', 'Native American', 'Asian & Pacific Islander')),
    # average attractiveness across waves
    attract = rowMeans(select(., ends_with('ir1')), na.rm = T)) %>%
  # sampling weights, cluster variable: school ID, and strata variable: region
  full_join(ah$wt4 %>% select(aid, wt4 = gswgt4_2, clusterid = psuscid, strataid = region)) %>%
  mutate(strataid = factor(strataid, labels = c('West', 'Midwest', 'South', 'Northeast'))) %>%
  # W5 race information
  full_join(ah$w5 %>% select(aid, starts_with('h5od4'), h5od8)) %>%
  mutate(
    # best race
    w3nhw = factor(h3od4a, labels = c('NA','NHW')),
    w3nhb = factor(h3od4b, labels = c('NA','NHB')),
    w3aina = factor(h3od4c, labels = c('NA','AINA')),
    w3aapi = factor(h3od4d, labels = c('NA','AAPI')),
    w3all = paste(w3nhw, w3nhb, w3aina, w3aapi, sep = ', '),
    w3all = case_when(w3all=='NA, NA, NA, NA' ~ NA, T ~ str_remove_all(w3all, '^(NA, )*|(, NA)*')),
    h3od6 = factor(ifelse(h3od6<=4, h3od6, NA), labels = c('NHW', 'NHB', 'AINA', 'AAPI')),
    best_race = case_when(
      h3od2==1 ~ 'Hispanic',
      h3od2==0 & !is.na(w3all) & !grepl(',', w3all) ~ w3all,
      h3od2==0 & !is.na(w3all) & str_count(w3all, ',')==1 & grepl('NHW', w3all) ~ 
        str_remove(w3all, '^(NHW, )*|(, NHW)*'),
      h3od2==0 & !is.na(w3all) & str_count(w3all, ',')==1 & !grepl('NHW', w3all) ~ h3od6,
      h3od2==0 & !is.na(w3all) & str_count(w3all, ',')>1 & !is.na(h3od6) ~ h3od6,
      T ~ NA
    ),
    w1nhw = factor(h1gi6a==1, labels = c('NA','NHW')),
    w1nhb = factor(h1gi6b==1, labels = c('NA','NHB')),
    w1aina = factor(h1gi6c==1, labels = c('NA','AINA')),
    w1aapi = factor(h1gi6d==1, labels = c('NA','AAPI')),
    w1all = paste(w1nhw, w1nhb, w1aina, w1aapi, sep = ', '),
    w1all = case_when(w1all=='NA, NA, NA, NA' ~ NA, T ~ str_remove_all(w1all, '^(NA, )*|(, NA)*')),
    h1gi8 = factor(ifelse(h1gi8<=4, h1gi8, NA), labels = c('NHW', 'NHB', 'AINA', 'AAPI')),
    best_race = case_when(
      is.na(best_race) & h1gi4==1 ~ 'Hispanic',
      is.na(best_race) & h1gi4==0 & !is.na(w1all) & !grepl(',', w1all) ~ w1all,
      is.na(best_race) & h1gi4==0 & !is.na(w1all) & str_count(w1all, ',')==1 & grepl('NHW', w1all) ~ 
        str_remove(w1all, '^(NHW, )*|(, NHW)*'),
      is.na(best_race) & h1gi4==0 & !is.na(w1all) & str_count(w1all, ',')==1 & !grepl('NHW', w1all) ~ h1gi8,
      is.na(best_race) & h1gi4==0 & !is.na(w1all) & str_count(w1all, ',')>1 & !is.na(h1gi8) ~ h1gi8,
      T ~ best_race
    ),
    w5nhw = factor(h5od4a==1, labels = c('NA','NHW')),
    w5nhb = factor(h5od4b==1, labels = c('NA','NHB')),
    w5aina = factor(h5od4f==1, labels = c('NA','AINA')),
    w5aapi = factor(h5od4d==1|h5od4e==1, labels = c('NA','AAPI')),
    w5all = paste(w5nhw, w5nhb, w5aina, w5aapi, sep = ', '),
    w5all = case_when(w5all=='NA, NA, NA, NA' ~ NA, T ~ str_remove_all(w5all, '^(NA, )*|(, NA)*')),
    h5od8 = case_when(
      h5od8==1 ~ 'NHW',
      h5od8==2 ~ 'NHB',
      h5od8==21 ~ 'AINA',
      h5od8 %in% c(9:15,17:18,20) ~ 'AAPI',
      h5od8 %in% c(3:7) ~ 'Hispanic'
    ),
    best_race = case_when(
      is.na(best_race) & h5od4c==1 ~ 'Hispanic',
      is.na(best_race) & h5od4c==0 & !is.na(w5all) & !grepl(',', w5all) ~ w5all,
      is.na(best_race) & h5od4c==0 & !is.na(w5all) & str_count(w5all, ',')==1 & grepl('NHW', w5all) ~ 
        str_remove(w5all, '^(NHW, )*|(, NHW)*'),
      is.na(best_race) & h5od4c==0 & !is.na(w5all) & str_count(w5all, ',')==1 & !grepl('NHW', w5all) ~ h5od8,
      is.na(best_race) & h5od4c==0 & !is.na(w5all) & str_count(w5all, ',')>1 & !is.na(h5od8) ~ h5od8,
      T ~ best_race
    ) %>% factor(levels = c('NHW', 'NHB', 'AINA', 'AAPI', 'Hispanic'))
  ) %>%
  # W1 constructed SES variables
  full_join(ah$ses %>% select(aid, nhood_ses = nhood1_d, fam_ses = sespc_al)) %>%
  select(aid, sex, byear, bmonth, paste0('age',1:4), edu, edu2, racec, skin, attract,
         attract1 = h1ir1, attract2 = h2ir1, attract3 = h3ir1, attract4 = h4ir1,
         matches('intid[1-4]'), wt4, clusterid, strataid, irace, irace3, irace4, 
         best_race, nhood_ses, fam_ses)


# 2.3 create inverse probability weights
ppheno <- ppheno %>%
  mutate(dna = as.numeric(aid %in% link$aid),
         dna = ifelse(!aid %in% ah$w1$aid, NA, dna))
write_dta(ppheno, 'interm_pheno.dta')

# run in STATA
# use interm_pheno.dta, clear
# logit dna i.racec i.sex byear i.edu i.strataid [pweight = wt4], or
# scalar loglik = e(ll)
# scalar pr2 = e(r2_p)
# predict prob, pr
# gen ipwgt = wt4 / prob
# margins, dydx(*) post
# estadd scalar loglik = loglik
# estadd scalar pr2 = pr2
# est store m
# esttab m using "tables/appendix_table8.csv", replace nogaps compress label ///
#   star(* 0.05 ** 0.01 *** 0.001) b(%20.3f) se(%20.3f) obslast eqlabel(none) ///
#   scalars("pr2 Pseudo R2") coeflabels(byear "Birth year (wave 4)") nonumbers mtitles("AME")
# save all_pheno.dta, replace

ppheno <- ppheno %>%
  left_join(read_dta('all_pheno.dta') %>% select(aid, ipwgt))

pheno <- link %>% 
  transmute(aid, FID = as.character(gfam), IID = as.character(gid)) %>%
  left_join(ppheno)
summary(pheno)


# 2.4 combine genetic & survey data
all <- dat %>%
  map(~{
    .x <- .x %>%
      left_join(pheno) %>%
      # drop HapMap3/HGDP reference samples from supervised ADMIXTURE
      filter(!is.na(aid)) %>%
      left_join(fam_origin)
  })

# X1: Sub-Saharan African
# X2: European
# X3: East Asian
# X4: Indigenous American
# X5: Middle Eastern
# X6: Oceania
colnames(all$us2)[1:2] <- c('X2', 'X1')
colnames(all$us3)[1:3] <- c('X3', 'X2', 'X1')
colnames(all$us4)[1:4] <- c('X2', 'X4', 'X1', 'X3')
colnames(all$us5)[1:5] <- c('X1', 'X5', 'X4', 'X2', 'X3')
colnames(all$us6)[1:6] <- c('X2', 'X5', 'X3', 'X6', 'X4', 'X1')
colnames(all$s2)[1:2] <- c('X2', 'X1')
colnames(all$s3)[1:3] <- c('X2', 'X3', 'X1')
colnames(all$s4)[1:4] <- c('X4', 'X2', 'X3', 'X1')
colnames(all$s5)[1:5] <- c('X5', 'X4', 'X2', 'X3', 'X1')
colnames(all$s6)[1:6] <- c('X6', 'X5', 'X4', 'X2', 'X3', 'X1')


# 2.5 create ADMIXTURE SE data
colnames(se$us2)[1:2] <- c('X2', 'X1')
colnames(se$us3)[1:3] <- c('X3', 'X2', 'X1')
colnames(se$us4)[1:4] <- c('X2', 'X4', 'X1', 'X3')
colnames(se$us5)[1:5] <- c('X1', 'X5', 'X4', 'X2', 'X3')
colnames(se$us6)[1:6] <- c('X2', 'X5', 'X3', 'X6', 'X4', 'X1')
colnames(se$s2)[1:2] <- c('X2', 'X1')
colnames(se$s3)[1:3] <- c('X2', 'X3', 'X1')
colnames(se$s4)[1:4] <- c('X4', 'X2', 'X3', 'X1')
colnames(se$s5)[1:5] <- c('X5', 'X4', 'X2', 'X3', 'X1')
colnames(se$s6)[1:6] <- c('X6', 'X5', 'X4', 'X2', 'X3', 'X1')

se <- se %>%
  map(~{
    .x <- .x %>%
      rename_at(vars(starts_with('X')), function(x) x = paste0(x, '_se') )
  })

b <- all %>%
  map(~{
    .x <- .x %>% 
      rename_at(vars(starts_with('X')), function(x) x = paste0(x, '_b') )
  })

se_fam <- list()
for (i in 2:6) {
  se_fam[[paste0('us',i)]] <- us_fam
}
se_fam <- c(se_fam, s_fam)
se_dat <- map2(se, se_fam, function(x,y) cbind(x,y))

all_stats <- map2(b, se_dat, function(x,y) left_join(x,y)) %>%
  bind_rows(.id = 'admixture') %>%
  select(aid, admixture, starts_with('X'))



# 3. Export data ---------------------------------------------------------------
saveRDS(all, 'all_data.rds')
saveRDS(all_stats, 'all_stats.rds')
saveRDS(ppheno, 'all_pheno.rds')

