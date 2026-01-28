# 0. Init & import data --------------------------------------------------------
library(here); library(haven); library(readxl)
library(tidyverse)
library(caret)
library(kableExtra); library(modelsummary); library(gtsummary); library(htmlTable)

setwd(str_remove(here(), '/scripts'))

subregion_label <- c('No other country', 'North Africa', 'Sub-Saharan Africa', 'British Isles',
                  'East Europe', 'North Europe', 'South Europe', 'West Europe',
                  'East Asia', 'Oceania', 'Philippines', 'South Asia', 'Southeast Asia', 
                  'West Asia', 'America', 'Caribbean', 'Central America', 'Canada',
                  'Mexico', 'Native America', 'South America',  'Other')
region_label <- c('Sub-Saharan Africa', 'Europe', 'Eastern Asia & Oceania', 
                       'Middle East & Southern Asia', 'The Americas', 'USA & Canada')

# read in Add Health survey and linkage data
wave1 <- read_xpt('allwave1.xpt') %>% rename_all(tolower)
wave3 <- read_xpt('wave3.xpt') %>% rename_all(tolower)
wave5 <- read_xpt('wave5.xpt') %>% rename_all(tolower)
link <- read_xpt('GID_link.xpt') %>% rename_all(tolower)

# read in prepared codebook and metadata
recode <- read_excel('geoanc_codebook.xlsx', na = '.') %>%
  transmute(value = value, country = factor(harmonized), 
            subregion = factor(subregion, levels = subregion_label))
meta <- read_excel('multilevel_regions.xlsx') %>%
  mutate(subregion = factor(subregion, levels = subregion_label), 
         region = factor(region, levels = region_label))

country_label <- levels(recode$country)



# 1. prepare geographic ancestry --------------------------------------------------------
# 1.1 label the self-reported geographic ancestries at country and subregion levels
merged <- link %>% 
  left_join(wave3 %>% select(aid, matches('h3od7'), h3od8), by = 'aid') %>%
  mutate_at(vars(matches('h3od7'), 'h3od8'), ~ ifelse(.>=990|.==203, NA, .))

v <- c('a', 'b', 'c', 'd')
for (i in (1:4)) {
  t <- v[[i]]
  d <- recode %>%
    rename_with(~ paste0(.,'_',t))
  var1 <- paste0('h3od7', t)
  var2 <- paste0('value', '_', t)
  merged <- merge(merged, d, by.x = var1, by.y = var2, all.x = T)
}
merged <- merge(merged, recode, by.x = 'h3od8', by.y = 'value', all.x = T)


# 1.2 basic processing and recategorization
all <- merged %>%
  mutate(
    # single ancestral country reported
    single_anc = case_when(
      h3od7b==0 ~ 1,
      is.na(h3od7b) | is.na(h3od7a) ~ NA, # same with or without is.na(h3od7a)
      T ~ 0
    ),
    # best geographic ancestry
    subregion = case_when(
      is.na(subregion) & !is.na(subregion_a) & single_anc==1 ~ subregion_a,
      T ~ subregion
    ), # 2728 reported America, 34 reported Canada
    h3od8 = case_when(
      is.na(h3od8) & !is.na(h3od7a) & single_anc==1 ~ h3od7a,
      T ~ h3od8
    ),
    country = case_when(
      is.na(country) & !is.na(country_a) & single_anc==1 ~ country_a,
      T ~ country
    ),
    across(starts_with('subregion_'), 
           ~ case_when(. %in% c('No other country', 'America', 'Canada', 'Other') ~ NA, T ~ .)),
    across(starts_with('country_'), 
           ~ case_when(. %in% c('No other country', 'America', 'Canada', 'Other') ~ NA, T ~ .)),
    # reported countries in SSA as one of geographic ancestries
    include_ssa = ifelse(subregion_a=='Sub-Saharan Africa' | subregion_b=='Sub-Saharan Africa' |
                           subregion_c=='Sub-Saharan Africa' | subregion_d=='Sub-Saharan Africa', 1, NA),
    # reported countries in The Americas as one of geographic ancestries
    include_america = case_when(
      subregion_a %in% c('Caribbean', 'Central America', 'Mexico', 'Native America', 'South America') ~ 1,
      subregion_b %in% c('Caribbean', 'Central America', 'Mexico', 'Native America', 'South America') ~ 1,
      subregion_c %in% c('Caribbean', 'Central America', 'Mexico', 'Native America', 'South America') ~ 1,
      subregion_d %in% c('Caribbean', 'Central America', 'Mexico', 'Native America', 'South America') ~ 1,
      T ~ NA
    ),
    # reported countries in Europe as one of geographic ancestries
    include_europe = case_when(
      subregion_a %in% c('British Isles', 'East Europe', 'North Europe', 'South Europe', 'West Europe') ~ 1,
      subregion_b %in% c('British Isles', 'East Europe', 'North Europe', 'South Europe', 'West Europe') ~ 1,
      subregion_c %in% c('British Isles', 'East Europe', 'North Europe', 'South Europe', 'West Europe') ~ 1,
      subregion_d %in% c('British Isles', 'East Europe', 'North Europe', 'South Europe', 'West Europe') ~ 1,
      T ~ NA
    ),
    # reported countries in Middle East & Southern Asia as one of geographic ancestries
    include_mesa = case_when(
      subregion_a %in% c('Middle East', 'South Asia') ~ 1,
      subregion_b %in% c('Middle East', 'South Asia') ~ 1,
      subregion_c %in% c('Middle East', 'South Asia') ~ 1,
      subregion_d %in% c('Middle East', 'South Asia') ~ 1,
      T ~ NA
    ),
    # reported countries in Eastern Asia & Oceania as one of geographic ancestries
    include_eao = case_when(
      subregion_a %in% c('East Asia', 'Southeast Asia', 'Philippines', 'Oceania') ~ 1,
      subregion_b %in% c('East Asia', 'Southeast Asia', 'Philippines', 'Oceania') ~ 1,
      subregion_c %in% c('East Asia', 'Southeast Asia', 'Philippines', 'Oceania') ~ 1,
      subregion_d %in% c('East Asia', 'Southeast Asia', 'Philippines', 'Oceania') ~ 1,
      T ~ NA
    )
  ) %>%
  # recategorize people who reported America/Canada as best ancestry and only one non-America/Canada ancestry
  mutate(
    subregion2 = case_when(
      subregion %in% c('America', 'Canada') & !is.na(subregion_a) & 
        rowSums(!is.na(.[paste0('subregion_',c('b','c','d'))]))==0 ~ subregion_a,
      subregion %in% c('America', 'Canada') & !is.na(subregion_b) & 
        rowSums(!is.na(.[paste0('subregion_',c('a','c','d'))]))==0 ~ subregion_b,
      subregion %in% c('America', 'Canada') & !is.na(subregion_c) & 
        rowSums(!is.na(.[paste0('subregion_',c('a','b','d'))]))==0 ~ subregion_c,
      subregion %in% c('America', 'Canada') & !is.na(subregion_d) & 
        rowSums(!is.na(.[paste0('subregion_',c('a','b','c'))]))==0 ~ subregion_d,
      T ~ subregion
    ),
    country2 = case_when(
      country %in% c('America', 'Canada') & !is.na(country_a) & 
        rowSums(!is.na(.[paste0('country_',c('b','c','d'))]))==0 ~ country_a,
      country %in% c('America', 'Canada') & !is.na(country_b) & 
        rowSums(!is.na(.[paste0('country_',c('a','c','d'))]))==0 ~ country_b,
      country %in% c('America', 'Canada') & !is.na(country_c) & 
        rowSums(!is.na(.[paste0('country_',c('a','b','d'))]))==0 ~ country_c,
      country %in% c('America', 'Canada') & !is.na(country_d) & 
        rowSums(!is.na(.[paste0('country_',c('a','b','c'))]))==0 ~ country_d,
      T ~ country
    )
  )


# 1.3 create a fractional version of geographic ancestry
# define a function for assigning a person's observation weight to different ancestries
ancestry_recode_USCA <- function(var) {
  # for respondents whose self-reported best geographic ancestry is US/Canada
  us_canada <- all %>%
    # n = 2762 (2728 + 34)
    filter(subregion %in% c('America', 'Canada')) %>%
    # reshape to dummy variables
    pivot_longer(cols = starts_with(paste0(var,'_')), names_to = 'names', values_to = 'category') %>%
    group_by(aid, category) %>%
    summarise(count = n()) %>%
    pivot_wider(names_from = category, values_from = count, values_fill = 0) %>%
    ungroup() %>%
    # whether reported US/Canada only or US/Canada & Other
    mutate(Other = 0, no_other_anc = `NA`==4) %>%
    # create weights for reported ancestry groups, excluding "Other"
    mutate(sum = rowSums(.[3:(length(.)-1)])) %>%
    mutate(across(3:(length(.)-2), ~ ifelse(no_other_anc, 0, ./sum))) %>%
    left_join(all) %>%
    # for those only reported US/Canada or US/Canada & Other
    mutate(America = ifelse(no_other_anc & .data[[var]]=='America', 1, 0),
           Canada = ifelse(no_other_anc & .data[[var]]=='Canada', 1, 0)) %>%
    select(-c(`NA`, sum))
 
  # for respondents whose self-reported best geographic ancestry is not US/Canada or missing
  non_us_canada <- all %>%
    cbind(predict(dummyVars(paste0('~ ', var), all, sep = NULL), newdata = all) %>%
            as.data.frame() %>%
            rename_with(~ str_remove(., var))) %>%
    filter(!(.data[[var]] %in% c('America', 'Canada'))) %>%
    # n = 7212
    mutate(no_other_anc = NA) %>% 
    select(-`No other country`)
  
  if (length(us_canada)!=length(non_us_canada)) {
    to_add <- names(non_us_canada)[!names(non_us_canada) %in% names(us_canada)]
    us_canada[to_add] <- 0
  }
  
  all_dum <- rbind(us_canada, non_us_canada)
  
  return(all_dum)
}

# create the fractional version of ancestry
all_dum <- ancestry_recode_USCA('subregion') %>%
  mutate(
    `REGION Sub-Saharan Africa` = `Sub-Saharan Africa`,
    `REGION Europe` = `British Isles` + `East Europe` + `North Europe` + `South Europe` + `West Europe`,
    `REGION Eastern Asia & Oceania` = `East Asia` + `Oceania` + `Philippines` + `Southeast Asia`,
    `REGION Middle East & Southern Asia` = `North Africa` + `West Asia` + `South Asia`,
    `REGION The Americas` = `Caribbean` + `Central America` + `Mexico` + `Native America` + `South America`,
  ) %>%
  rename_with(~ paste('SUBREG', .), all_of(subregion_label[-1])) %>%
  left_join(
    ancestry_recode_USCA('country') %>%
      rename_with(~ paste('CNTR', .), all_of(country_label[!country_label %in% 'No other country']))
  )

# some respondents reported multiple ancestral countries that belong to a single subregion
for (r in subregion_label[-1]) {
  all_dum[all_dum[[paste0('SUBREG ',r)]]==1 & all_dum$subregion2 %in% c('Canada','America','Other'),]$subregion2 <- r
}


# 1.4 impute ancestry using race information
# for respondents who were missing on the ancestry questions or who reported America/Canada only
all_na <- all_dum %>%
  filter(is.na(subregion) | `SUBREG America`+`SUBREG Canada`==1) %>% # n = 1783 + 1680
  # merge in racial identity and detailed Asian/Hispanic/Pacific Islander questions
  left_join(wave3 %>% select(aid, matches('^h3od[2-5]'))) %>%
  left_join(wave5 %>% select(aid, matches('^h5od[4-7]'))) %>%
  left_join(wave1 %>% select(aid, matches('^h1gi[4-7]'))) %>%
  # first impute the fractional version using W3-W5-W1
  # prioritize W3
  mutate(
    across(matches('^(h3od[2-5]|h5od[4-7]|h1gi[4-7])'), ~ ifelse(.>1, NA, .)),
    # Chinese, Japanese, & Korean
    `SUBREG East Asia` = case_when(h3od5a==1 | h3od5c==1 | h3od5e==1 ~ 1, T ~ `SUBREG East Asia`),
    `CNTR China` = case_when(h3od5a==1 ~ 1, T ~ `CNTR China`),
    `CNTR Japan` = case_when(h3od5c==1 ~ 1, T ~ `CNTR Japan`),
    `CNTR North and South Korea` = case_when(h3od5e==1 ~ 1, T ~ `CNTR North and South Korea`),
    # Philippines
    `SUBREG Philippines` = case_when(h3od5b==1 ~ 1, T ~ `SUBREG Philippines`),
    # Indian
    `SUBREG South Asia` = case_when(h3od5d==1 ~ 1, T ~ `SUBREG South Asia`),
    # Vietnames
    `SUBREG Southeast Asia` = case_when(h3od5f==1 ~ 1, T ~ `SUBREG Southeast Asia`),
    # Mexican/Mexican American, Chicano/Chicana
    `SUBREG Mexico` = case_when(h3od3a==1 | h3od3b==1 ~ 1, T ~ `SUBREG Mexico`),
    # Cuban/Cuban American, Puerto Rican
    `SUBREG Caribbean` = case_when(h3od3c==1 | h3od3d==1 ~ 1, T ~ `SUBREG Caribbean`),
    `CNTR Cuba` = case_when(h3od3c==1 ~ 1, T ~ `CNTR Cuba`),
    `CNTR Puerto Rico` = case_when(h3od3d==1 ~ 1, T ~ `CNTR Puerto Rico`),
    # Central American
    `SUBREG Central America` = case_when(h3od3e==1 ~ 1, T ~ `SUBREG Central America`),
    # Sub-Saharan Africa
    `SUBREG Sub-Saharan Africa` = case_when(h3od4b==1 ~ 1, T ~ `SUBREG Sub-Saharan Africa`),
  ) %>%
  # whether successfully imputed using W3
  mutate(imputed = rowSums(.[2:20], na.rm = T) > 0, imputed3 = imputed) %>%
  # if not, use W5
  mutate(
    # Chinese, Japanese, & Korean
    `SUBREG East Asia` = case_when(!imputed & (h5od6b==1 | h5od6d==1 | h5od6e==1) ~ 1, T ~ `SUBREG East Asia`),
    `CNTR China` = case_when(!imputed & h5od6b==1 ~ 1, T ~ `CNTR China`),
    `CNTR Japan` = case_when(!imputed & h5od6d==1 ~ 1, T ~ `CNTR Japan`),
    `CNTR North and South Korea` = case_when(!imputed & h5od6e==1 ~ 1, T ~ `CNTR North and South Korea`),
    # Philippines
    `SUBREG Philippines` = case_when(!imputed & h5od6c==1 ~ 1, T ~ `SUBREG Philippines`),
    # Indian
    `SUBREG South Asia` = case_when(!imputed & h5od6a==1 ~ 1, T ~ `SUBREG South Asia`),
    # Vietnames
    `SUBREG Southeast Asia` = case_when(!imputed & h5od6f==1 ~ 1, T ~ `SUBREG Southeast Asia`),
    # Mexican/Mexican American, Chicano/Chicana
    `SUBREG Mexico` = case_when(!imputed & h5od5a==1 ~ 1, T ~ `SUBREG Mexico`),
    # Cuban/Cuban American, Puerto Rican
    `SUBREG Caribbean` = case_when(!imputed & (h5od5c==1 | h5od5b==1) ~ 1, T ~ `SUBREG Caribbean`),
    `CNTR Cuba` = case_when(!imputed & h5od5c==1 ~ 1, T ~ `CNTR Cuba`),
    `CNTR Puerto Rico` = case_when(!imputed & h5od5b==1 ~ 1, T ~ `CNTR Puerto Rico`),
    # Central American
    `SUBREG Central America` = case_when(!imputed & h5od5d==1 ~ 1, T ~ `SUBREG Central America`),
    # South American
    `SUBREG South America` = case_when(!imputed & h5od5e==1 ~ 1, T ~ `SUBREG South America`),
    `CNTR Unknown South America` = `SUBREG South America`,
    # Native Hawaiian, Samoan, Guamanian or Chamorro, Other Pacific Islander
    `SUBREG Oceania` = case_when(!imputed & (h5od7a==1 | h5od7b==1 | h5od7c==1 | h5od7d==1) ~ 1, T ~ `SUBREG Oceania`),
    `CNTR Pacific Islands, Australia, and Atlantic Islands` = `SUBREG Oceania`,
    # Sub-Saharan Africa
    `SUBREG Sub-Saharan Africa` = case_when(!imputed & h5od4b==1 ~ 1, T ~ `SUBREG Sub-Saharan Africa`),
  ) %>%
  # whether successfully imputed using W5
  mutate(imputed = rowSums(.[2:20], na.rm = T) > 0, imputed5 = !imputed3 & imputed) %>%
  # if not, use W1
  mutate(
    # Chinese, Japanese, & Korean
    `SUBREG East Asia` = case_when(!imputed & (h1gi7a==1 | h1gi7c==1 | h1gi7e==1) ~ 1, T ~ `SUBREG East Asia`),
    `CNTR China` = case_when(!imputed & h1gi7a==1 ~ 1, T ~ `CNTR China`),
    `CNTR Japan` = case_when(!imputed & h1gi7c==1 ~ 1, T ~ `CNTR Japan`),
    `CNTR North and South Korea` = case_when(!imputed & h1gi7e==1 ~ 1, T ~ `CNTR North and South Korea`),
    # Philippines
    `SUBREG Philippines` = case_when(!imputed & h1gi7b==1 ~ 1, T ~ `SUBREG Philippines`),
    `CNTR Philippines` = `SUBREG Philippines`,
    # Indian
    `SUBREG South Asia` = case_when(!imputed & h1gi7d==1 ~ 1, T ~ `SUBREG South Asia`),
    `CNTR India` = `SUBREG South Asia`,
    # Vietnames
    `SUBREG Southeast Asia` = case_when(!imputed & h1gi7f==1 ~ 1, T ~ `SUBREG Southeast Asia`),
    `CNTR Vietnam` = `SUBREG Southeast Asia`,
    # Mexican/Mexican American, Chicano/Chicana
    `SUBREG Mexico` = case_when(!imputed & (h1gi5a==1 | h1gi5b==1) ~ 1, T ~ `SUBREG Mexico`),
    `CNTR Mexico` = `SUBREG Mexico`,
    # Cuban/Cuban American, Puerto Rican
    `SUBREG Caribbean` = case_when(!imputed & (h1gi5c==1 | h1gi5d==1) ~ 1, T ~ `SUBREG Caribbean`),
    `CNTR Cuba` = case_when(!imputed & h1gi5c==1 ~ 1, T ~ `CNTR Cuba`),
    `CNTR Puerto Rico` = case_when(!imputed & h1gi5d==1 ~ 1, T ~ `CNTR Puerto Rico`),
    # Central American
    `SUBREG Central America` = case_when(!imputed & h1gi5e==1 ~ 1, T ~ `SUBREG Central America`),
    `CNTR Unknown Central America` = `SUBREG Central America`,
    # Sub-Saharan Africa
    `SUBREG Sub-Saharan Africa` = case_when(!imputed & h1gi6b==1 ~ 1, T ~ `SUBREG Sub-Saharan Africa`),
    `CNTR Unknown Sub-Saharan Africa` = `SUBREG Sub-Saharan Africa`
  ) %>%
  # recode NA to 0 for imputed individuals
  mutate(
    imputed = rowSums(.[2:20], na.rm = T) > 0,
    imputed1 = !imputed3 & !imputed5 & imputed,
    across(c(paste0('SUBREG ',c('America','Canada')), paste0('CNTR ',c('America','Canada'))), 
           ~ case_when(.==1 & imputed ~ 0, T ~ .)),
    across(c(starts_with('SUBREG '), starts_with('CNTR ')), 
           ~ case_when(is.na(.) & imputed ~ 0, T ~ .))
  ) %>%
  # create fractional weights for imputed ancestries
  mutate(
    # create subregion-level weights
    sum1 = rowSums(select(., starts_with('SUBREG '))),
    across(starts_with('SUBREG '), ~ case_when(imputed ~ ./sum1, T ~ .)),
    # create country-level weights
    sum2 = rowSums(select(., starts_with('CNTR '))),
    across(starts_with('CNTR '), ~ case_when(imputed ~ ./sum2, T ~ .)),
    # create region-level weights
    `REGION Sub-Saharan Africa` = `SUBREG Sub-Saharan Africa`,
    `REGION Europe` = `SUBREG British Isles` + `SUBREG East Europe` + `SUBREG North Europe` + `SUBREG South Europe` + `SUBREG West Europe`,
    `REGION Eastern Asia & Oceania` = `SUBREG East Asia` + `SUBREG Oceania` + `SUBREG Philippines` + `SUBREG Southeast Asia`,
    `REGION Middle East & Southern Asia` = `SUBREG North Africa` + `SUBREG West Asia` + `SUBREG South Asia`,
    `REGION The Americas` = `SUBREG Caribbean` + `SUBREG Central America` + `SUBREG Mexico` + `SUBREG Native America` + `SUBREG South America`,
  ) %>%
  # impute best ancestry for respondents with only one imputed subregion/country
  mutate(
    subregion2 = case_when(
      # Chinese, Japanese, & Korean
      `SUBREG East Asia`==1 ~ 'East Asia',
      # Philippines
      `SUBREG Philippines`==1 ~ 'Philippines',
      # Indian
      `SUBREG South Asia`==1 ~ 'South Asia',
      # Vietnames
      `SUBREG Southeast Asia`==1 ~ 'Southeast Asia',
      # Mexican/Mexican American, Chicano/Chicana
      `SUBREG Mexico`==1 ~ 'Mexico',
      # Cuban/Cuban American, Puerto Rican
      `SUBREG Caribbean`==1 ~ 'Caribbean',
      # Central American
      `SUBREG Central America`==1 ~ 'Central America',
      # South American
      `SUBREG South America`==1 ~ 'South America',
      # Native Hawaiian, Samoan, Guamanian or Chamorro, Other Pacific Islander
      `SUBREG Oceania`==1 ~ 'Oceania',
      # Sub-Saharan Africa
      `SUBREG Sub-Saharan Africa`==1 ~ 'Sub-Saharan Africa',
      T ~ subregion2
    ),
    country2 = case_when(
      # Chinese, Japanese, & Korean
      `CNTR China`==1 ~ 'China',
      `CNTR Japan`==1 ~ 'Japan',
      `CNTR North and South Korea`==1 ~ 'North and South Korea',
      # Philippines
      `CNTR Philippines`==1 ~ 'Philippines',
      # Indian
      `CNTR India`==1 ~ 'India',
      # Vietnames
      `CNTR Vietnam`==1 ~ 'Vietnam',
      # Mexican/Mexican American, Chicano/Chicana
      `CNTR Mexico`==1 ~ 'Mexico',
      # Cuban/Cuban American, Puerto Rican
      `CNTR Cuba`==1 ~ 'Cuba',
      `CNTR Puerto Rico`==1 ~ 'Puerto Rico',
      # Central American
      `CNTR Unknown Central America`==1 ~ 'Unknown Central America',
      # South American
      `CNTR Unknown South America`==1 ~ 'Unknown South America',
      # Native Hawaiian, Samoan, Guamanian or Chamorro, Other Pacific Islander
      `CNTR Pacific Islands, Australia, and Atlantic Islands`==1 ~ 'Pacific Islands, Australia, and Atlantic Islands',
      # Sub-Saharan Africa
      `CNTR Unknown Sub-Saharan Africa`==1 ~ 'Unknown Sub-Saharan Africa',
      T ~ country2
    )
  )

all_dum2 <- all_dum %>%
  filter(!is.na(subregion), `SUBREG America`+`SUBREG Canada`==0) %>%
  mutate(`CNTR Unknown Central America` = 0) %>%
  rbind(all_na %>% select(all_of(names(all_dum)), 'CNTR Unknown Central America')) %>%
  # merge in region-level variables
  left_join(meta) %>%
  left_join(meta %>% rename(region2 = region, subregion2 = subregion))

# some respondents have multiple ancestral subregions that belong to a single best region
for (r in region_label[-6]) {
  condition1 <- all_dum2[[paste0('REGION ',r)]]==1
  condition1 <- ifelse(is.na(condition1), F, condition1)
  condition2 <- all_dum2$region2!=r
  condition2 <- ifelse(is.na(condition2), T, condition2)
  all_dum2[condition1 & condition2,]$region2 <- r
}

# further impute the best ancestral region for Hispanics and non-Hispanic Black
all_reg <- all_dum2 %>% 
  filter(is.na(region2) | region2=='USA & Canada') %>% # n = 2084
  left_join(all_na %>% select(aid, paste0('imputed',c(1,3,5)))) %>%
  left_join(wave3 %>% select(aid, h3od2, h3od4b)) %>%
  left_join(wave5 %>% select(aid, h5od4c, h5od4b)) %>%
  left_join(wave1 %>% select(aid, h1gi4, h1gi6b)) %>%
  mutate(
    across(c(h3od2, h3od4b, h5od4c, h5od4b, h1gi4, h1gi6b), ~ ifelse(.>1, NA, .)),
    region3 = case_when(
      imputed3 & h3od2==1 ~ 'The Americas',
      imputed3 & h3od2==0 & h3od4b==1 ~ 'Sub-Saharan Africa',
      T ~ region2
    )
  ) %>%
  mutate(
    region3 = case_when(
      imputed5 & h5od4c==1 ~ 'The Americas',
      imputed5 & h5od4c==0 & h5od4b==1 ~ 'Sub-Saharan Africa',
      T ~ region3
    )
  ) %>%
  mutate(
    region3 = case_when(
      imputed1 & h1gi4==1 ~ 'The Americas',
      imputed1 & h1gi4==0 & h1gi6b==1 ~ 'Sub-Saharan Africa',
      T ~ region3
    )
  )


# 1.5 create the final geographic ancestry data
all_final <- all_dum2 %>%
  left_join(all_reg %>% select(aid, region3)) %>%
  mutate(region3 = case_when(!is.na(region3) ~ region3,
                                  T ~ region2) %>% factor(levels = levels(all_dum2$region2)),
         subregion2 = case_when(region3=='Sub-Saharan Africa' ~ 'Sub-Saharan Africa', T ~ subregion2),
         country2 = case_when(region3=='Sub-Saharan Africa' ~ 'Unknown Sub-Saharan Africa', T ~ country2),
         include_ssa = ifelse(`REGION Sub-Saharan Africa` > 0 |
                                region3=='Sub-Saharan Africa', 1, include_ssa),
         include_america = ifelse(`REGION The Americas` > 0 |
                                    region3=='The Americas', 1, include_america),
         include_europe = ifelse(`REGION Europe` > 0 |
                                   region3=='Europe', 1, include_europe),
         include_mesa = ifelse(`REGION Middle East & Southern Asia` > 0 |
                                 region3=='Middle East & Southern Asia', 1, include_mesa),
         include_eao = ifelse(`REGION Eastern Asia & Oceania` > 0 |
                                region3=='Eastern Asia & Oceania', 1, include_eao),
         # in a very rare case: a respondent reported Hispanic origin but was missing
         # on the detailed background question
         # make the fractional and best versions consistent
         across(matches('^SUBREG|CNTR '), 
                ~ case_when(`REGION The Americas`==0 & region3=='The Americas' ~ NA, T ~ .)),
         across(paste0('REGION ',region_label[1:4]), 
                ~ case_when(`REGION The Americas`==0 & region3=='The Americas' ~ 0, T ~ .)),
         `REGION The Americas` = case_when(`REGION The Americas`==0 & region3=='The Americas' ~ 1,
                                             T ~ `REGION The Americas`)) %>%
  # combine West Asia and North Africa and rename to Middle East
  mutate(subregion = case_when(subregion %in% c('North Africa', 'West Asia') ~ 'Middle East', T ~ subregion),
         subregion2 = case_when(subregion2 %in% c('North Africa', 'West Asia') ~ 'Middle East', T ~ subregion2),
         `SUBREG Middle East` = `SUBREG West Asia` + `SUBREG North Africa`) %>%
  select(-c(`SUBREG West Asia`, `SUBREG North Africa`))

all_final %>% select(starts_with('SUBREG '), -paste0('SUBREG ',c('America','Canada','Other'))) %>%
  drop_na() %>% sum() # 8055
all_final %>% select(starts_with('REGION ')) %>% drop_na() %>% sum() # 8056



# 2. Export --------------------------------------------------------------------
saveRDS(all_final %>% select(aid, region, subregion, country, region2, subregion2,
                             country2, region3, single_anc, starts_with('REGION '),
                             starts_with('SUBREG '), starts_with('CNTR '), starts_with('include_')),
        'family_origin_final.rds')

# export region meta data
meta %>%
  left_join(recode %>% 
              select(-value) %>%
              rbind(data.frame(subregion = 'Central America', country = 'Unknown Central America')) %>%
              unique()) %>%
  mutate(subregion = case_when(subregion %in% c('North Africa', 'West Asia') ~ 'Middle East', T ~ subregion)) %>%
  saveRDS('region_meta_final.rds')

