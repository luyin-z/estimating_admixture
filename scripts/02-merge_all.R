# 0. Init ----------------------------------------------------------------------
library(here); library(readr); library(haven)
library(tidyverse)

setwd(str_remove(here(), '/scripts'))



# 1. Import data ---------------------------------------------------------------
# 1.1 ADMIXTURE point estimates
us <- readRDS('unsupervised_output.rds')
s <- readRDS('supervised_output.rds')
dat <- c(us, s)


# 1.2 linkage data
link <- read_xpt('GID_link.xpt') %>% rename_all(tolower) # n = 9974


# 1.3 survey data
ah <- list()
ah$w1 <- read_xpt('allwave1.xpt') %>% rename_all(tolower)
for (i in 2:5) {
  ah[[paste0('w',i)]] <- read_xpt(paste0('wave', i, '.xpt')) %>% rename_all(tolower)
}
ah$wt4 <- read_xpt('weights4.xpt') %>% rename_all(tolower)
ah$ses <- read_xpt('conses.xpt') %>% rename_all(tolower)


# 1.4 geographic ancestry
fam_origin <- readRDS('family_origin_final.rds')



# 2. Prepare data --------------------------------------------------------------
# 2.1 key variables from survey
ppheno <- ah$w4 %>%
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
    # racial classification
    h4ir4,
    # interviewer ID
    intid4
  ) %>%
  # W3 variables
  full_join(
    ah$w3 %>% 
      mutate(age3 = (iyear3-h3od1y)*12 + imonth3-h3od1m,
             across(c('h3od2', starts_with('h3od4')), function(x) x = ifelse(x>1, NA, x)),
             w3nhw = factor(h3od4a, labels = c('NA','NHW')),
             w3nhb = factor(h3od4b, labels = c('NA','NHB')),
             w3aina = factor(h3od4c, labels = c('NA','AINA')),
             w3aapi = factor(h3od4d, labels = c('NA','AAPI')),
             w3all = paste(w3nhw, w3nhb, w3aina, w3aapi, sep = ', '),
             w3all = case_when(w3all=='NA, NA, NA, NA' ~ NA, T ~ str_remove_all(w3all, '^(NA, )*|(, NA)*')),
             race = case_when(
               !grepl(',', w3all) ~ w3all,
               is.na(w3all) ~ NA,
               T ~ 'Multiracial') %>% factor(levels = c('NHW', 'NHB', 'AINA', 'AAPI', 'Multiracial')), 
             hispanic = factor(h3od2, labels = c('Not', 'Hispanic')),
             racec = case_when(hispanic=='Hispanic' ~ 'Hispanics', hispanic=='Not' ~ race, T ~ NA) %>%
               factor(levels = c('NHW', 'NHB', 'AINA', 'AAPI', 'Hispanics', 'Multiracial'),
                      labels = c('Non-Hispanic White', 'Non-Hispanic Black', 'Non-Hispanic Native American',
                                 'Non-Hispanic Asian & Pacific Islander', 'Hispanic (of any race)',
                                 'Non-Hispanic Multiple Races')),
             skin = h3ir17 %>% factor(labels = c('Black', 'Dark Brown', 'Medium Brown', 'Light Brown', 'White'))) %>%
      select(aid, intid3, age3, racec, skin, h3od2, starts_with('h3od4'), h3od6, h3ir1, h3ir4)) %>%
  # W1 and W2 variables
  full_join(ah$w1 %>% select(aid, intid1 = intid, iyear, imonth, h1gi1y, h1gi1m, h1gi9)) %>%
  full_join(ah$w2 %>% select(aid, intid2, iyear2, imonth2, h2gi1y, h2gi1m, h2ir1)) %>%
  mutate(intid1 = as.numeric(intid1), intid2 = as.numeric(intid2), 
         across(c(h1gi1y, h1gi1m), function(x) x = ifelse(x>=96, NA, x)),
         across(c(h2gi1y, h2gi1m), function(x) x = ifelse(x>=98, NA, x)),
         age1 = (iyear-h1gi1y)*12 + imonth-h1gi1m,
         age2 = (iyear2-h2gi1y)*12 + imonth2-h2gi1m) %>%
  mutate(
    # interviewer's racial classification of respondents
    irace = ifelse(h1gi9>=6, NA, h1gi9),
    irace = factor(irace, labels = c('White', 'Black', 'Native American', 'Asian & Pacific Islander', 'Other')),
    irace3 = factor(h3ir4, labels = c('White', 'Black', 'Native American', 'Asian & Pacific Islander')),
    irace4 = factor(h4ir4, labels = c('White', 'Black', 'Native American', 'Asian & Pacific Islander'))) %>%
  full_join(ah$wt4 %>% select(aid, wt4 = gswgt4_2, clusterid = psuscid, strataid = region)) %>% # cluster variable: school ID; strata variable: region
  mutate(strataid = factor(strataid, labels = c('West', 'Midwest', 'South', 'Northeast'))) %>%
  # merge in W1 constructed SES variables
  full_join(ah$ses %>% transmute(aid = as.character(aid), nhood_ses = nhood1_d, fam_ses = sespc_al)) %>%
  select(aid, sex, byear, bmonth, paste0('age',1:4), edu, racec, skin,
         matches('intid[1-4]'), wt4, clusterid, strataid, irace, irace3, irace4,
         nhood_ses, fam_ses)


# 2.2 create inverse probability weights
ppheno <- ppheno %>%
  mutate(dna = as.numeric(aid %in% link$aid),
         dna = ifelse(!aid %in% ah$w1$aid, NA, dna))
write_dta(ppheno, 'interm_pheno.dta')

# run the STATA script
system('"C:/Program Files/StataNow19/StataSE-64.exe" -b do "P:/AddHealth/Contract/27110801-Conley/Work/Luyin Zhang/admixture/scripts/03-create_ipw.do"')

ppheno <- ppheno %>%
  left_join(read_dta('all_pheno.dta') %>% select(aid, prob, ipwgt))

pheno <- link %>% 
  transmute(aid, FID = as.character(gfam), IID = as.character(gid)) %>%
  left_join(ppheno)
summary(pheno)


# 2.3 combine genetic & survey data
all <- dat %>%
  map(~{
    .x <- .x %>%
      mutate(across(c(FID,IID), as.character)) %>%
      left_join(pheno) %>%
      filter(!is.na(aid)) %>% # drop 1kg/hgdp ref sample for supervised admixture outputs
      left_join(fam_origin)
  })



# 3. Export data ---------------------------------------------------------------
saveRDS(all, 'all_data.rds')
saveRDS(ppheno, 'all_pheno.rds')

