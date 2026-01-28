# 0. Init ----------------------------------------------------------------------
library(here); library(haven)
library(tidyverse); library(survey); library(scales)
library(estimatr); library(caret); library(fixest); library(lme4)
library(lmtest); library(sandwich); library(marginaleffects)
library(modelsummary); library(gtsummary); library(flextable); library(riverplot)
library(grDevices); library(ftExtra); library(ggh4x); library(ggtern)

setwd(str_remove(here(), '/scripts'))

region_label <- c('Europe', 'The Americas', 'Sub-Saharan Africa', 'Eastern Asia & Oceania',
                       'Middle East & South Asia')
ancestry_label <- c('Sub−Saharan African', 'European', 'East Asian', 'Indigenous American')
ancestry_label2 <- c('AFR', 'EUR', 'EAS', 'IAM')

source('scripts/00-functions.R')

set_gtsummary_theme(list('style_number-arg:big.mark' = ''))
fig_format <- '.pdf'



# 1. Import data ---------------------------------------------------------------
all <- readRDS('all_data.rds')
se <- readRDS('supervised_k45_se_output.rds')
pheno <- readRDS('all_pheno.rds')
fam_origin <- readRDS('family_origin_final.rds')
meta <- readRDS('region_meta_final.rds')
best_race <- readRDS('best_race.rds') # from another project



# 2. Tables --------------------------------------------------------------------
# 2.1 Table 1 (main text): descriptive statistics
# left panel of table 1: unweighted
d <- all$s4 %>%
  mutate(has_reg = 1-as.numeric(is.na(region3)),
         has_subreg = 1-as.numeric(is.na(subregion2)),
         has_reg = case_when(region3=='USA & Canada' ~ 0, T ~ has_reg),
         has_subreg = case_when(subregion2 %in% c('America','Canada','Other') ~ 0, T ~ has_subreg))

left <- d %>%
  select(sex, byear, edu, racec, strataid, has_reg, has_subreg, region, skin, fam_ses, nhood_ses) %>%
  tbl_summary(  
    statistic = list(all_continuous() ~ '{mean}-{sd}',
                     all_categorical() ~ '{p}%- '),
    digits = everything() ~ 2,
    # type = list(all_categorical() ~ 'categorical'),
    label = list(sex ~ 'Gender',
                 byear ~ 'Birth year',
                 edu ~ 'Educational attainment',
                 racec ~ 'Race & enthnicity',
                 strataid ~ 'High school census region',
                 has_reg ~ 'Has ancestral region data',
                 has_subreg ~ 'Has ancestral sub-region data',
                 region ~ 'Best ancestry',
                 skin ~ 'Skin tone',
                 fam_ses ~ 'W1 social origins score',
                 nhood_ses ~ 'W1 neighborhood disadvantage')                            
  ) %>%
  # do not show missingness rows
  remove_row_type(type = 'missing') %>%
  # add sample sizes
  add_n() %>%
  # remove footnote
  modify_footnote(c(all_stat_cols()) ~ NA) %>%
  # modify header
  modify_header(stat_0 ~ '**Mean (SD) / %**')

# add statistics for region/subregion data availability
n1 <- d %>% select(starts_with('reg ')) %>% sum(na.rm = T)
n2 <- d %>% select(starts_with('SUBREG '), -paste0('SUBREG ',c('America','Canada','Other'))) %>% sum(na.rm = T)
left_new <- left$table_body %>%
  filter(variable %in% c('has_reg', 'has_subreg')) %>%
  mutate(
    stat_0 = ifelse(variable=='has_reg', n1, n2),
    stat_0 = paste0(sprintf('%.2f', stat_0/as.numeric(n)*100), '%- ')
  )

# add statistics for the fractional version of geographic ancestry
left_add <- d %>%
  select(starts_with('REGION ')) %>%
  drop_na() %>%
  mutate(n = sum(rowSums(.)),
         across(starts_with('REGION '), function(x) x = sum(x)/sum(rowSums(.)))) %>%
  unique() %>%
  gather('region', 'prop', starts_with('REGION ')) %>%
  transmute(variable = 'region2', 
            var_type = 'categorical',
            var_label = 'Best ancestry - adjusted',
            row_type = 'level',
            label = str_remove(region, 'REGION '),
            stat_0 = paste0(sprintf('%.2f', prop*100), '%- '),
            n) %>%
  add_row(variable = 'region2', var_type = 'categorical', var_label = 'Best ancestry - adjusted',
          row_type = 'label', label = var_label, stat_0 = NA, n = unique(.$n)) %>%
  mutate(n = ifelse(row_type=='label', n, NA))

# complete left panel
left_all <- left %>%
  modify_table_body(
    ~.x %>% 
      filter(!(variable %in% c('skin','fam_ses','nhood_ses','has_reg','has_subreg'))) %>%
      rbind(left_add[c(6,1:5),]) %>%
      rbind(left_new) %>%
      rbind(.x %>% filter(variable %in% c('skin','fam_ses','nhood_ses')))
  )

# right panel of table 1: weighted
right <- svydesign(data = d %>% filter(!is.na(clusterid), !is.na(strataid), !is.na(ipwgt)),
                   id = ~ clusterid, strata = ~ strataid, weights = ~ ipwgt) %>%
  tbl_svysummary(     
    statistic = list(all_continuous() ~ '{mean}-{sd}',       
                     all_categorical() ~ '{p}%- '),
    digits = everything() ~ 2,
    type = all_categorical() ~ 'categorical',
    label = list(sex ~ 'Gender',
                 byear ~ 'Birth year',
                 edu ~ 'Educational attainment',
                 racec ~ 'Race & enthnicity',
                 strataid ~ 'High school census region',
                 has_reg ~ 'Has ancestral region data',
                 has_subreg ~ 'Has ancestral sub-region data',
                 region ~ 'Best ancestry',
                 skin ~ 'Skin tone',
                 fam_ses ~ 'W1 social origins score',
                 nhood_ses ~ 'W1 neighborhood disadvantage'),
    include = c(sex, byear, edu, racec, strataid, has_reg, has_subreg, region, skin, fam_ses, nhood_ses)
  ) %>%
  # do not show missingness rows
  remove_row_type(type = 'missing') %>%
  # add sample sizes
  add_n() %>%
  # remove footnote
  modify_footnote(c(all_stat_cols()) ~ NA) %>%
  # modify header
  modify_header(stat_0 ~ '**Mean (SD) / %**')

# add statistics for region/subregion data availability
n1 <- d %>% 
  mutate(across(starts_with('REGION '), function(x) x = ipwgt*x)) %>%
  select(starts_with('REGION ')) %>% 
  sum(na.rm = T)
n2 <- d %>% 
  mutate(across(starts_with('SUBREG '), function(x) x = ipwgt*x)) %>% 
  select(starts_with('SUBREG '), -paste0('SUBREG ',c('America','Canada','Other'))) %>% 
  sum(na.rm = T)
right_new <- right$table_body %>%
  filter(variable %in% c('has_reg', 'has_subreg'), row_type=='label') %>%
  mutate(
    stat_0 = ifelse(variable=='has_reg', n1, n2),
    stat_0 = paste0(sprintf('%.2f', stat_0/as.numeric(n)*100), '%- '),
    var_type = 'dichotomous'
  )

# add statistics for the fractional version of geographic ancestry
right_svy <- svydesign(data = d %>%
                         mutate(n = 1) %>%
                         filter(!is.na(clusterid), !is.na(strataid), !is.na(ipwgt)),
                       id = ~ clusterid, strata = ~ strataid, weights = ~ ipwgt)

result <- list()
for (reg in d %>% select(starts_with('REGION ')) %>% names()) {
  result[[reg]] <- svytotal(as.formula(paste0('~`', reg, '`')), na.rm = TRUE, design = right_svy)[[1]]
}

right_add <- result %>%
  bind_rows(.id = 'region') %>%
  mutate(n = rowSums(.)) %>%
  gather('region', 'sum', starts_with('REGION ')) %>%
  transmute(variable = 'region2', 
            var_type = 'categorical',
            var_label = 'Best ancestry - adjusted',
            row_type = 'level',
            label = str_remove(region, 'REGION '),
            stat_0 = paste0(sprintf('%.2f', sum/n*100), '%- '),
            n) %>%
  add_row(variable = 'region2', var_type = 'categorical', var_label = 'Best ancestry - adjusted',
          row_type = 'label', label = var_label, stat_0 = NA, n = unique(.$n)) %>%
  mutate(n = ifelse(row_type=='label', n, NA))

# complete right panel
wt_n <- svytotal(~n, design = right_svy)[[1]]
## sample size with valid weights
N <-  d %>% filter(!is.na(clusterid), !is.na(strataid), !is.na(ipwgt)) %>% nrow

right_all <- right %>%
  modify_table_body(
    ~.x %>% 
      filter(!(variable %in% c('skin','fam_ses','nhood_ses','has_reg','has_subreg'))) %>%
      rbind(right_add[c(6,1:5),]) %>%
      rbind(right_new) %>%
      rbind(.x %>% filter(variable %in% c('skin','fam_ses','nhood_ses'))) %>%
      # rescale sample size
      mutate(n = sprintf('%.0f', as.numeric(n)/wt_n*N),
             n = ifelse(n=='NA', NA, n))
  )

# combine the left and right panels of table 1
tbl1 <- left_all %>%
  modify_table_body(
    ~.x %>%
      left_join(right_all$table_body %>% rename(n1 = n, stat_1 = stat_0)) %>%
      separate(stat_0, c('stat_0a','stat_0b'), '-') %>%
      separate(stat_1, c('stat_1a','stat_1b'), '-')
  ) %>% 
  modify_header(n1 = '**N**', stat_0a = '**Mean / %**', stat_0b = '**SD**',
                stat_1a = '**Mean / %**', stat_1b = '**SD**') %>%
  modify_spanning_header(c(n, stat_0a, stat_0b) ~ '**Unweighted**', c(n1, stat_1a, stat_1b) ~ '**Weighted**')

tbl1 %>% as_gt() %>% gt::gtsave(file = 'tables/table1.docx')


# 2.2 Table 2: supervised K = 4 (the fraction version)
k4 <- all[c('us4', 's4')]
colnames(k4$us4)[1:4] <- paste0('us4_', colnames(k4$us4)[1:4])
colnames(k4$s4)[1:4] <- paste0('s4_', colnames(k4$s4)[1:4]) %>%
  str_replace('X1', 'Sub−Saharan African') %>%
  str_replace('X2', 'European') %>%
  str_replace('X3', 'East Asian') %>%
  str_replace('X4', 'Indigenous American')

table_k4 <- k4$us4 %>%
  select(aid, starts_with('REGION '), starts_with('SUBREG '), starts_with('us4_')) %>%
  left_join(k4$s4 %>% select(aid, starts_with('s4_')))

# unsupervised and supervised GSP summary statistics
upper <- summarytable2(table_k4, paste0('us4_X', 1:4))
lower <- summarytable2(table_k4, paste0('s4_', ancestry_label))

# reliability
reliability <- se$s4 %>%
  mutate(stats = 'se') %>%
  rbind(all$s4 %>% 
          select(FID, IID, paste0('X',1:4)) %>%
          mutate(stats = 'b')) %>%
  gather(X, value, paste0('X',1:4)) %>%
  spread(stats, value) %>%
  group_by(X) %>%
  summarise(
    mean_b = mean(b), tss_b = sum((b-mean_b)^2), sum_var = sum(se^2),
    rel = 1-sum_var/tss_b
  ) %>%
  ungroup() %>%
  select(X, rel) %>%
  spread(X, rel) %>%
  mutate(across(X1:X4, ~ sprintf('%.3f',.)))

# correlation
cor_k4 <- table_k4 %>%
  select(starts_with('us4'), starts_with('s4')) %>%
  mutate(X1 = sprintf('%.3f', cor(us4_X1, `s4_Sub−Saharan African`)),
         X2 = sprintf('%.3f', cor(us4_X2, `s4_European`)),
         X3 = sprintf('%.3f', cor(us4_X3, `s4_East Asian`)),
         X4 = sprintf('%.3f', cor(us4_X4, `s4_Indigenous American`))) %>% 
  select(X1:X4) %>%
  unique()

# make panelA: unsupervised 4
panelA <- upper %>%
  mutate(Region = case_when(Subregion=='Total' ~ Region, T ~ ''),
         `N = 8055` = paste0(sprintf('%.0f',`N = 8055`), ' (',
                             sprintf('%.0f', prop*100), '%)')) %>%
  rename_with(function(x) x = str_replace(x, 'us4', 'Unsupervised 4'), starts_with('us4')) %>%
  select(-prop) %>%
  rename(`Sub-Region` = Subregion) %>%
  flextable() %>%
  split_header() %>%
  merge_h(part = 'header') %>%
  set_header_labels(values = c(rep('',3), '*P*^1^','*P*^2^','*P*^3^','*P*^4^')) %>%
  align(align = 'center', part = 'header') %>%
  colformat_md(part = 'header') %>%
  align(j = 4:7, align = 'center') %>%
  autofit()

# make panelB: supervised 4 + reliability
panelB <- lower %>%
  mutate(Region = case_when(Subregion=='Total' ~ Region, T ~ ''),
         `N = 8055` = paste0(sprintf('%.0f',`N = 8055`), ' (',
                             sprintf('%.0f', prop*100), '%)')) %>%
  rename_with(function(x) x = str_remove(x, 's4_'), starts_with('s4_')) %>%
  select(-prop) %>%
  rbind(reliability %>%
          transmute(Region = 'Reliability', Subregion = '', `N = 8055` = '', 
                    `Sub−Saharan African` = X1, `European` = X2, `East Asian` = X3,
                    `Indigenous American` = X4)) %>%
  rename(`Sub-Region` = Subregion) %>%
  flextable() %>%
  add_header_row(values = c('','*P*^AFR^','*P*^EUR^','*P*^EAS^','*P*^IAM^'),
                 colwidths = c(3,rep(1,4))) %>%
  add_header_row(values = c('', 'Supervised 4'), colwidths = c(3,4)) %>%
  align(align = 'center', part = 'header') %>%
  colformat_md(part = 'header') %>%
  merge_at(i = nrow(lower)+1, j = 1:2) %>%
  align(i = nrow(lower)+1, j = 1, align = 'center') %>%
  align(j = 4:7, align = 'center') %>%
  bold(i = nrow(lower)+1, j = 1) %>%
  autofit()

# make panelC: correlation
panelC <- cor_k4 %>%
  transmute(` ` = '', `  ` = '', `   ` = '', 
            `Correlation_X1 - Sub−Saharan African`  = X1, `Correlation_X2 - European` = X2,
            `Correlation_X3 - East Asian` = X3, `Correlation_X4 - Indigenous American` = X4) %>%
  flextable() %>%
  set_header_labels(values = c('','','','*P*^1^,*P*^AFR^','*P*^2^,*P*^EUR^','*P*^3^,*P*^EAS^',
                               '*P*^4^,*P*^IAM^')) %>%
  add_header_row(values = c('', 'Unsupervised-Supervised Correlation'), colwidths = c(3,4)) %>%
  align(align = 'center', part = 'header') %>%
  colformat_md(part = 'header') %>%
  align(j = 4:7, align = 'center') %>%
  autofit()

save_as_docx(panelB, panelC, path = 'tables/table2.docx')
save_as_docx(panelA, path = 'tables/appendix_unadx4_tbl.docx')


# 2.3 Table 3 (main text): regression table for skin tone
# prepare data
d <- all$s4 %>% 
  mutate(region = region3 %>% factor(),
         skin = 6 - as.numeric(skin),
         darkbrown = as.numeric(skin>=4),
         lightbrown = as.numeric(skin>=2),
         skin2 = case_when(region %in% c('The Americas', 'Europe') ~ lightbrown,
                           region=='Sub-Saharan Africa' ~ darkbrown, T ~ NA),
         female = as.numeric(sex=='Female'),
         across(paste0('X',1:4), ~ .*10),
         ga_comb = paste(include_europe, include_america, include_ssa, include_eao,
                         include_mesa, sep = ''))

# create a variable for ancestral country fixed effects
d <- d %>%
  select(starts_with('CNTR ')) %>%
  drop_na() %>%
  mutate_all(~ifelse(.>0, 1, .)) %>%
  transmute(cga_comb = do.call(paste, c(., sep = ' '))) %>%
  cbind(d %>% drop_na(starts_with('CNTR ')) %>% select(aid)) %>%
  right_join(d)

# include people who reported SSA or The Americas as one of their geographic ancestries
d <- d %>%
  filter(include_ssa==1) %>%
  mutate(region = 'Sub-Saharan Africa') %>%
  rbind(d %>%
          filter(include_america==1) %>%
          mutate(region = 'The Americas')) %>%
  rbind(d %>%
          filter(include_europe==1) %>%
          mutate(region = 'Europe')) %>%
  mutate(skin2 = case_when(region %in% c('The Americas', 'Europe') ~ lightbrown,
                           region=='Sub-Saharan Africa' ~ darkbrown, T ~ NA))

exp1 <- 'skin ~ female + age3 + I(age3^2) + X1 + X3 + X4'
exp2 <- 'skin ~ female + age3 + I(age3^2) + X1 + X3 + X4 + fam_ses + nhood_ses'
exp3 <- 'skin ~ female + age3 + I(age3^2) + X1 + X3 + X4 + fam_ses + nhood_ses | intid3'
exp4 <- 'skin ~ female + age3 + I(age3^2) + X1 + X3 + X4 + fam_ses + nhood_ses | intid3 + cga_comb'

# OLS regression
reg <- lapply(c('Sub-Saharan Africa','The Americas'), function(ga) {
  list(
    lm(exp1, d %>% filter(region==ga)), 
    lm(exp2, d %>% filter(region==ga)), 
    feols(as.formula(exp3), 
          d %>% 
            filter(region==ga) %>%
            drop_na(skin, female, age3, fam_ses, nhood_ses, intid3) %>%
            group_by(intid3) %>%
            filter(n()>1) %>% 
            ungroup()),
    feols(as.formula(exp4), 
          d %>% 
            filter(region==ga) %>%
            drop_na(skin, female, age3, fam_ses, nhood_ses, intid3, cga_comb) %>%
            group_by(intid3) %>%
            filter(n()>1) %>% 
            group_by(cga_comb) %>%
            filter(n()>1) %>% 
            group_by(intid3) %>%
            filter(n()>1) %>% 
            ungroup())
  )
}) %>% 
  set_names(c('Sub-Saharan Africa','The Americas'))

# logistic regression
reg2 <- lapply(c('Sub-Saharan Africa','The Americas'), function(ga) {
  list(
    glm(exp1 %>% str_replace('skin','skin2'), family = binomial(link='logit'), 
        d %>% filter(region==ga) %>% mutate(across(paste0('X',1:4), ~ ./10))),
    glm(exp2 %>% str_replace('skin','skin2'), family = binomial(link='logit'),
        d %>% filter(region==ga) %>% mutate(across(paste0('X',1:4), ~ ./10)))
  )
}) %>% set_names(c('Sub-Saharan Africa','The Americas'))

# create AMEs and pseudo-R2 for logistic models
ame <- list(`Sub-Saharan Africa` = list(), `The Americas` = list())
pseudo_r2 <- list()
for (ga in c('Sub-Saharan Africa','The Americas')) {
  ame[[ga]][[1]] <- avg_slopes(reg2[[ga]][[1]], 
                               d %>% 
                                 filter(region==ga) %>% 
                                 mutate(across(paste0('X',1:4), ~ ./10)) %>%
                                 drop_na(skin2, female, age3))
  ame[[ga]][[1]]$estimate <- ame[[ga]][[1]]$estimate/10
  ame[[ga]][[1]]$std.error <- ame[[ga]][[1]]$std.error/10
  ame[[ga]][[2]] <- avg_slopes(reg2[[ga]][[2]], 
                               d %>% 
                                 filter(region==ga) %>% 
                                 mutate(across(paste0('X',1:4), ~ ./10)) %>%
                                 drop_na(skin2, female, age3, fam_ses, nhood_ses))
  ame[[ga]][[2]]$estimate <- ame[[ga]][[2]]$estimate/10
  ame[[ga]][[2]]$std.error <- ame[[ga]][[2]]$std.error/10
  pseudo_r2[[ga]] <- lapply(reg2[[ga]], function(m) performance::r2(m)) %>%
    bind_cols()
}
pseudo_r2 <- pseudo_r2 %>%
  bind_cols() %>%
  mutate(V1 = 'Pseudo R2')
names(pseudo_r2)[1:4] <- paste0('V',2:5)

# average skin tone
mean_skin <- lapply(c('Sub-Saharan Africa','The Americas'), function(ga) {
  list(
    d %>% 
      filter(region==ga) %>%
      drop_na(skin, female, age3) %>%
      summarise(across(c(skin,skin2), mean)),
    d %>% 
      filter(region==ga) %>%
      drop_na(skin, female, age3, fam_ses, nhood_ses) %>%
      summarise(across(c(skin,skin2), mean)),
    d %>% 
      filter(region==ga) %>%
      drop_na(skin, female, age3, fam_ses, nhood_ses, intid3) %>%
      group_by(intid3) %>%
      filter(n()>1) %>%
      ungroup() %>%
      summarise(across(c(skin,skin2), mean)),
    d %>% 
      filter(region==ga) %>%
      drop_na(skin, female, age3, fam_ses, nhood_ses, intid3, cga_comb) %>%
      group_by(intid3) %>%
      filter(n()>1) %>%
      group_by(cga_comb) %>%
      filter(n()>1) %>%
      group_by(intid3) %>%
      filter(n()>1) %>%
      ungroup() %>%
      summarise(across(c(skin,skin2), mean))
  ) %>%
    bind_rows() %>%
    t()
}) %>% 
  do.call(cbind, .) %>%
  as.data.frame() %>%
  mutate(label = 'Mean outcome',
         across(!label, ~ sprintf('%.3f',.))) %>%
  select(label, V1:V8)
names(mean_skin) <- paste0('V',1:9)

# upper panel of table 3: OLS regression
rows <- c('Basic Controls', rep('X',8),
          'Socioeconomic Controls', '', rep('X',3), '', rep('X',3),
          'Interviewer FEs', rep(c('','','X','X'),2),
          'Ancestry FEs', rep(c('','','','X'),2)) %>%
  matrix(nrow = 4, byrow = T) %>%
  as.data.frame() %>%
  rbind(mean_skin['skin',])
attr(rows, 'position') <- c(7:11)

tbl3_reg <- modelsummary(c(reg$`Sub-Saharan Africa`, reg$`The Americas`), 
                         estimate = '{estimate}{stars}',
                         coef_omit = '^(?!X[1345])', 
                         gof_omit = 'IC|RMSE|Std.Errors|FE|Log.Lik.|F|Adj',
                         add_rows = rows,
                         output = 'flextable') %>%
  set_header_labels(values = c('', rep(paste('Model',1:4),2))) %>%
  add_header_row(values = c('', 'Skin tone\n(OLS)','Skin tone\n(OLS)'), 
                 colwidths = c(1,4,4), top = T) %>%
  add_header_row(values = c('', 'Self-Reported Ancestral Sub-Saharan Africans',
                            'Self-Reported Ancestral Americans'), 
                 colwidths = c(1,4,4), top = T) %>%
  align(align = 'center', part = 'header') %>%
  compose(i = 1, j = 1, value = as_paragraph(as_i('P'), as_sup('AFR'), ' (10 pp)')) %>%
  compose(i = 3, j = 1, value = as_paragraph(as_i('P'), as_sup('EAS'), ' (10 pp)')) %>%
  compose(i = 5, j = 1, value = as_paragraph(as_i('P'), as_sup('IAM'), ' (10 pp)')) %>%
  hline(i = 10) %>%
  vline(j = 5) %>%
  autofit()

# lower panel of table 3: AMEs from logistic regression
rows2 <- c('Basic Controls', rep('X',4),
           'Socioeconomic Controls', '', 'X', '', 'X') %>%
  matrix(nrow = 2, byrow = T) %>%
  as.data.frame() %>%
  rbind(mean_skin['skin2',paste0('V',c(1:3,6:7))] %>% rename(V4=V6,V5=V7)) %>%
  rbind(pseudo_r2 %>%
          mutate(across(V2:V5, ~sprintf('%.3f',.)))) 
names(rows2)[2:5] <- paste0('(',1:4,')')
rows2 <- rows2 %>% mutate(part = '', term = V1, statistic = '')

tbl3_ame <- modelsummary(c(ame$`Sub-Saharan Africa`, ame$`The Americas`),
                         estimate = '{estimate}{stars}',
                         coef_omit = '^(?!X[1345])', 
                         gof_omit = 'IC|F|RMSE|Log.Lik.',
                         add_rows = rows2,
                         output = 'data.frame') %>%
  mutate(part = ifelse(grepl('Pseudo',term), 'last', part)) %>%
  arrange(part!='estimates', part!='manual', part=='last') %>%
  mutate(term = ifelse(statistic=='std.error', '', term)) %>%
  select(-c(part, statistic)) %>%
  flextable() %>%
  set_header_labels(values = c('', paste('Model',rep(1:2,2)))) %>%
  add_header_row(values = c('', 'Dark brown or darker\n(Logistic)','Light brown or darker\n(Logistic)'), 
                 colwidths = c(1,2,2), top = T) %>%
  align(align = 'center', part = 'header') %>%
  compose(i = 1, j = 1, value = as_paragraph(as_i('P'), as_sup('AFR'), ' (10 pp)')) %>%
  compose(i = 3, j = 1, value = as_paragraph(as_i('P'), as_sup('EAS'), ' (10 pp)')) %>%
  compose(i = 5, j = 1, value = as_paragraph(as_i('P'), as_sup('IAM'), ' (10 pp)')) %>%
  hline(i = c(6,8)) %>%
  vline(j = 3) %>%
  autofit()

save_as_docx(tbl3_reg, tbl3_ame, path = 'tables/table3.docx')

# marginal R2 of genetic similarity proportions for continuous skin tone
skin_mr2 <- lapply(c('Sub-Saharan Africa','The Americas'), function(ga) {
  td <- d %>% filter(region==ga) %>% drop_na(skin, female, age3)
  r2 <- summary(reg[[ga]][[1]])[['r.squared']]
  r2o <- summary(lm(skin ~ female + age3 + I(age3^2), td))[['r.squared']]
  dr0 <- summary(lm(skin ~ X1 + X3 + X4, td))[['r.squared']]
  dr1 <- r2 - r2o
  combined <- data.frame(model = c('Null', 'Basic controls', 'Basic controls & others GSPs'),
                         term = 'Combined', marginal_r2 = sprintf('%.3f', c(dr0,dr1,dr1)))
  lapply(paste0('X',c(1,3:4)), function(x) {
    exp1 <- paste('skin ~ ', x)
    exp2 <- paste('skin ~ female + age3 + I(age3^2) +', x)
    exp3 <- paste('skin ~ female + age3 + I(age3^2) +', paste(setdiff(paste0('X',c(1,3:4)),x), collapse = ' + '))
    m1 <- lm(exp1, td)
    m2 <- lm(exp2, td)
    m3 <- lm(exp3, td)
    mr2_1 <- summary(m1)[['r.squared']]
    mr2_2 <- summary(m2)[['r.squared']] - r2o
    mr2_3 <- r2 - summary(m3)[['r.squared']]
    mr2 <- c(mr2_1, mr2_2, mr2_3)
    data.frame(model = c('Null', 'Basic controls', 'Basic controls & others GSPs'),
               term = x, marginal_r2 = sprintf('%.3f', mr2)) %>%
      return()
  }) %>%
    bind_rows() %>%
    rbind(combined) %>%
    mutate(term = case_when(term=='X1' ~ '*P*^AFR^', term=='X3' ~ '*P*^EAS^',
                            term=='X4' ~ '*P*^IAM^', T ~ term)) %>%
    return()
}) %>% 
  set_names(c('Sub-Saharan Africa','The Americas')) %>%
  bind_rows(.id = 'sample')


# 2.4 Table 4 (main text): regression table for racial classification
# prepare data
d <- all$s4 %>%
  rename(irace1 = irace) %>%
  pivot_longer(cols = matches('^age|^intid|^irace'),
               names_to = c('.value', 'wave'),
               names_pattern = '^(age|intid|irace)(\\d+)$',
               values_drop_na = TRUE,
               values_transform = list(age = as.numeric)) %>%
  arrange(aid) %>%
  mutate(region = region3 %>% factor(),
         age2 = age*age, 
         inw = 1 - as.numeric(irace=='White'),
         ib = as.numeric(irace=='Black'),
         across(paste0('X',1:4), ~ .*10),
         ga_comb = paste(include_ssa, include_europe, include_eao,
                         include_america, include_mesa, sep = ''))

# create a variable for ancestral country fixed effects
d <- all$s4 %>%
  select(starts_with('CNTR ')) %>%
  drop_na() %>%
  mutate_all(~ifelse(.>0, 1, .)) %>%
  transmute(cga_comb = do.call(paste, c(., sep = ' '))) %>%
  cbind(all$s5 %>% drop_na(starts_with('CNTR ')) %>% select(aid)) %>%
  right_join(d)

# include people who reported SSA or The Americas as one of their geographic ancestries
d <- d %>%
  filter(include_ssa==1) %>%
  mutate(region = 'Sub-Saharan Africa') %>%
  rbind(d %>%
          filter(include_america==1) %>%
          mutate(region = 'The Americas')) %>%
  rbind(d %>%
          filter(include_europe==1) %>%
          mutate(region = 'Europe'))

exp1 <- 'y ~ factor(sex) + age + age2 + X1 + X3 + X4'
exp2 <- 'y ~ factor(sex) + age + age2 + X1 + X3 + X4 + fam_ses + nhood_ses'
exp3 <- 'y ~ factor(sex) + age + age2 + X1 + X3 + X4 + skin + fam_ses + nhood_ses'
exp4 <- 'y ~ factor(sex) + age + age2 + X1 + X3 + X4 + skin + fam_ses + nhood_ses | intid'
exp5 <- 'y ~ factor(sex) + age + age2 + X1 + X3 + X4 + skin + fam_ses + nhood_ses | intid + cga_comb'

# linear probability models
reg <- list(d %>% filter(region=='Sub-Saharan Africa') %>% rename(y = ib), 
            d %>% filter(region=='The Americas' & wave=='1') %>% rename(y = inw)) %>%
  lapply(function(data) {
    list(
      lm_robust(as.formula(exp1), clusters = aid, data), 
      lm_robust(as.formula(exp2), clusters = aid, data), 
      lm_robust(as.formula(exp3), clusters = aid, data), 
      feols(as.formula(exp4), cluster = 'aid', 
            data %>%
              drop_na(sex, age, skin, fam_ses, nhood_ses, intid, y) %>%
              group_by(intid) %>%
              filter(n()>1) %>%
              ungroup()),
      feols(as.formula(exp5), cluster = 'aid', 
            data %>%
              drop_na(sex, age, skin, fam_ses, nhood_ses, intid, cga_comb, y) %>%
              group_by(intid) %>%
              filter(n()>1) %>%
              group_by(cga_comb) %>%
              filter(n()>1) %>%
              group_by(intid) %>%
              filter(n()>1) %>%
              ungroup())
    )
  }) %>%
  set_names(c('Sub-Saharan Africa', 'The Americas'))

# logistic regression
reg2 <- list(d %>% 
               filter(region=='Sub-Saharan Africa') %>% 
               mutate(y = ib, across(paste0('X',1:4), ~ ./10)), 
             d %>% filter(region=='The Americas' & wave=='1') %>% 
               mutate(y = inw, across(paste0('X',1:4), ~ ./10))) %>%
  lapply(function(data) {
    list(
      glm(exp1, data, family = binomial()),
      glm(exp2, data, family = binomial()),
      glm(exp3, data, family = binomial())
    )
  }) %>%
  set_names(c('Sub-Saharan Africa', 'The Americas'))

# create AMEs and pseudo-R2 for logistic models
ame <- list(`Sub-Saharan Africa` = list(), `The Americas` = list())
pseudo_r2 <- list()
for (ga in c('Sub-Saharan Africa','The Americas')) {
  if (ga=='The Americas') {
    td <- d %>% filter(wave=='1')
  } else {
    td <- d
  }
  ame[[ga]][[1]] <- avg_slopes(reg2[[ga]][[1]], 
                               td %>% 
                                 filter(region==ga) %>% 
                                 mutate(y = ifelse(grepl('Africa',region),ib,inw),
                                        across(paste0('X',1:4), ~ ./10)) %>%
                                 drop_na(y, sex, age))
  ame[[ga]][[1]]$estimate <- ame[[ga]][[1]]$estimate/10
  ame[[ga]][[1]]$std.error <- ame[[ga]][[1]]$std.error/10
  ame[[ga]][[2]] <- avg_slopes(reg2[[ga]][[2]], 
                               td %>% 
                                 filter(region==ga) %>% 
                                 mutate(y = ifelse(grepl('Africa',region),ib,inw),
                                        across(paste0('X',1:4), ~ ./10)) %>%
                                 drop_na(y, sex, age, fam_ses, nhood_ses))
  ame[[ga]][[2]]$estimate <- ame[[ga]][[2]]$estimate/10
  ame[[ga]][[2]]$std.error <- ame[[ga]][[2]]$std.error/10
  ame[[ga]][[3]] <- avg_slopes(reg2[[ga]][[3]], 
                               td %>% 
                                 filter(region==ga) %>% 
                                 mutate(y = ifelse(grepl('Africa',region),ib,inw),
                                        across(paste0('X',1:4), ~ ./10)) %>%
                                 drop_na(y, sex, age, fam_ses, nhood_ses, skin))
  ame[[ga]][[3]]$estimate <- ame[[ga]][[3]]$estimate/10
  ame[[ga]][[3]]$std.error <- ame[[ga]][[3]]$std.error/10
  pseudo_r2[[ga]] <- lapply(reg2[[ga]], function(m) performance::r2(m)) %>%
    bind_cols()
}
pseudo_r2 <- pseudo_r2 %>%
  bind_cols() %>%
  mutate(V1 = 'Pseudo R2')
names(pseudo_r2)[1:6] <- paste0('V',2:7)

# average racial classification
n_meany <- list(d %>% filter(region=='Sub-Saharan Africa') %>% rename(y = ib), 
                d %>% filter(region=='The Americas' & wave=='1') %>% rename(y = inw)) %>%
  lapply(function(data) {
    list(
      data %>% 
        drop_na(sex, age, y) %>%
        summarise(y = mean(y), n = length(unique(aid))),
      data %>% 
        drop_na(sex, age, fam_ses, nhood_ses, y) %>%
        summarise(y = mean(y), n = length(unique(aid))),
      data %>% 
        drop_na(sex, age, skin, fam_ses, nhood_ses, y) %>%
        summarise(y = mean(y), n = length(unique(aid))),
      data %>% 
        drop_na(sex, age, skin, fam_ses, nhood_ses, intid, y) %>%
        group_by(intid) %>%
        filter(n()>1) %>%
        ungroup() %>%
        summarise(y = mean(y), n = length(unique(aid))),
      data %>% 
        drop_na(sex, age, skin, fam_ses, nhood_ses, intid, cga_comb, y) %>%
        group_by(intid) %>%
        filter(n()>1) %>%
        group_by(cga_comb) %>%
        filter(n()>1) %>%
        group_by(intid) %>%
        filter(n()>1) %>%
        ungroup() %>%
        summarise(y = mean(y), n = length(unique(aid)))
    ) %>%
      bind_rows() %>%
      t()
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  rownames_to_column('label') %>%
  mutate(across(!label, ~ ifelse(label=='n', sprintf('%.0f',.), sprintf('%.3f',.))),
         label = ifelse(label=='n','Num.Inds.','Mean outcome')) %>%
  select(label, V1:V10)
names(n_meany) <- paste0('V',1:11)

# lower panel of table 4: linear probability models
rows <- c('Basic Controls', rep('X',10),
          'Socioeconomic Controls', '', rep('X',4), '', rep('X',4),
          'Skin Tone Control', rep(c('','','X','X','X'),2),
          'Interviewer FEs', rep(c('','','','X','X'),2),
          'Ancestry FEs', rep('',4), 'X', rep('',4), 'X') %>%
  matrix(nrow = 5, byrow = T) %>%
  rbind(n_meany) %>%
  as.data.frame()
attr(rows, 'position') <- c(7:13)

tbl4_reg <- modelsummary(c(reg$`Sub-Saharan Africa`, reg$`The Americas`),
                         estimate = '{estimate}{stars}',
                         coef_omit = '^(?!X[1345])', 
                         gof_omit = 'IC|RMSE|Std.Errors|FE|Log.Lik.|F|Adj',
                         add_rows = rows,
                         output = 'flextable') %>%
  set_header_labels(values = c('', paste('Model',rep(1:5,2)))) %>%
  add_header_row(values = c('', 'Classified as Black\n(LPM)', 'Classified as Non-White\n(LPM)'), 
                 colwidths = c(1,5,5), top = T) %>%
  align(align = 'center', part = 'header') %>%
  compose(i = 1, j = 1, value = as_paragraph(as_i('P'), as_sup('AFR'), ' (10 pp)')) %>%
  compose(i = 3, j = 1, value = as_paragraph(as_i('P'), as_sup('EAS'), ' (10 pp)')) %>%
  compose(i = 5, j = 1, value = as_paragraph(as_i('P'), as_sup('IAM'), ' (10 pp)')) %>%
  hline(i = 11) %>%
  vline(j = 6) %>%
  autofit()

# upper panel of table 4: AMEs from logistic regression
rows2 <- c('Basic Controls', rep('X',6),
           'Socioeconomic Controls', '', rep('X',2), '', rep('X',2),
           'Skin Tone Control', rep(' ',2), 'X', rep(' ',2), 'X') %>%
  matrix(nrow = 3, byrow = T) %>%
  as.data.frame() %>%
  rbind(n_meany[,paste0('V',c(1:4,7:9))] %>% rename(V5=V7,V6=V8,V7=V9)) %>%
  rbind(pseudo_r2 %>%
          mutate(across(V2:V7, ~sprintf('%.3f',.)))) 
names(rows2)[2:7] <- paste0('(',1:6,')')
rows2 <- rows2 %>% mutate(part = '', term = V1, statistic = '')

tbl4_ame <- modelsummary(c(ame$`Sub-Saharan Africa`, ame$`The Americas`),
                         estimate = '{estimate}{stars}',
                         coef_omit = '^(?!X[1345])', 
                         gof_omit = 'IC|F|RMSE|Log.Lik.',
                         add_rows = rows2,
                         output = 'data.frame') %>%
  mutate(part = ifelse(grepl('Pseudo',term), 'last', part)) %>%
  arrange(part!='estimates', part!='manual', part=='last') %>%
  mutate(term = ifelse(statistic=='std.error', '', term)) %>%
  select(-c(part, statistic)) %>%
  flextable() %>%
  set_header_labels(values = c('', paste('Model',rep(1:3,2)))) %>%
  add_header_row(values = c('', 'Classified as Black\n(Logistic)', 'Classified as Non-White\n(Logistic)'), 
                 colwidths = c(1,3,3), top = T) %>%
  add_header_row(values = c('', 'Self-Reported Ancestral Sub-Saharan Africans',
                            'Self-Reported Ancestral Americans'), colwidths = c(1,3,3)) %>%
  align(align = 'center', part = 'header') %>%
  compose(i = 1, j = 1, value = as_paragraph(as_i('P'), as_sup('AFR'), ' (10 pp)')) %>%
  compose(i = 3, j = 1, value = as_paragraph(as_i('P'), as_sup('EAS'), ' (10 pp)')) %>%
  compose(i = 5, j = 1, value = as_paragraph(as_i('P'), as_sup('IAM'), ' (10 pp)')) %>%
  hline(i = c(6,9)) %>%
  vline(j = 4) %>%
  autofit()

save_as_docx(tbl4_ame, tbl4_reg, path = 'tables/table4.docx')

# marginal R2 of genetic similarity proportions for racial classification
y_mr2 <- list(d %>% filter(region=='Sub-Saharan Africa') %>% rename(y = ib), 
              d %>% filter(region=='The Americas' & wave=='1') %>% rename(y = inw)) %>%
  lapply(function(data) {
    td <- data %>% drop_na(y, sex, age, skin)
    r2 <- summary(lm_robust(y ~ factor(sex) + age + age2 + X1 + X3 + X4, clusters = aid, td))[['r.squared']]
    r2o <- summary(lm_robust(y ~ factor(sex) + age + age2, clusters = aid, td))[['r.squared']]
    r2f1 <- summary(lm_robust(y ~ factor(sex) + age + age2 + skin + X1 + X3 + X4, clusters = aid, td))[['r.squared']]
    r2f2 <- summary(lm_robust(y ~ factor(sex) + age + age2 + skin, clusters = aid, td))[['r.squared']]
    dr0 <- summary(lm_robust(y ~ X1 + X3 + X4, clusters = aid, td))[['r.squared']]
    dr1 <- r2 - r2o
    dr2 <- r2f1 - r2f2
    combined <- data.frame(model = c('Null', 'Basic controls', 'Basic controls & others GSPs',
                                     'Basic controls, others GSPs, & skin tone'),
                           term = 'Combined', marginal_r2 = sprintf('%.3f', c(dr0,dr1,dr1,dr2)))
    lapply(paste0('X',c(1,3:4)), function(x) {
      exp1 <- paste('y ~ ', x)
      exp2 <- paste('y ~ factor(sex) + age + age2 +', x)
      exp3 <- paste('y ~ factor(sex) + age + age2 +', paste(setdiff(paste0('X',c(1,3:4)),x), collapse = ' + '))
      exp4 <- paste('y ~ factor(sex) + age + age2 + skin +', paste(setdiff(paste0('X',c(1,3:4)),x), collapse = ' + '))
      m1 <- lm_robust(as.formula(exp1), clusters = aid, td)
      m2 <- lm_robust(as.formula(exp2), clusters = aid, td)
      m3 <- lm_robust(as.formula(exp3), clusters = aid, td)
      m4 <- lm_robust(as.formula(exp4), clusters = aid, td)
      mr2_1 <- summary(m1)[['r.squared']]
      mr2_2 <- summary(m2)[['r.squared']] - r2o
      mr2_3 <- r2 - summary(m3)[['r.squared']]
      mr2_4 <- r2f1 - summary(m4)[['r.squared']]
      mr2 <- c(mr2_1, mr2_2, mr2_3, mr2_4)
      data.frame(model = c('Null', 'Basic controls', 'Basic controls & others GSPs',
                           'Basic controls, others GSPs, & skin tone'),
                 term = x, marginal_r2 = sprintf('%.3f', mr2)) %>%
        return()
    }) %>%
      bind_rows() %>%
      rbind(combined) %>%
      mutate(term = case_when(term=='X1' ~ '*P*^AFR^', term=='X3' ~ '*P*^EAS^',
                              term=='X4' ~ '*P*^IAM^', T ~ term)) %>%
      return()
  }) %>% 
  set_names(c('Sub-Saharan Africa','The Americas')) %>%
  bind_rows(.id = 'sample')

y_skin_mr2 <- list(d %>% filter(region=='Sub-Saharan Africa') %>% rename(y = ib), 
                   d %>% filter(region=='The Americas' & wave=='1') %>% rename(y = inw)) %>%
  lapply(function(data) {
    td <- data %>% drop_na(y, sex, age, skin)
    r2a <- summary(lm_robust(y ~ factor(sex) + age + age2 + X1 + X3 + X4, clusters = aid, td))[['r.squared']]
    r2b <- summary(lm_robust(y ~ factor(sex) + age + age2 + skin + X1 + X3 + X4, clusters = aid, td))[['r.squared']]
    dr <- r2b - r2a
    data.frame(model = 'Basic controls, others GSPs, & skin tone',
               term = 'Skin tone', marginal_r2 = sprintf('%.3f', dr)) %>%
      return()
  }) %>% 
  set_names(c('Sub-Saharan Africa','The Americas')) %>%
  bind_rows(.id = 'sample')

# Table 14 (appendix): marginal R2 of genetic similarity proportions for skin tone and racial classification
mr2_tbl <- skin_mr2 %>%
  mutate(outcome = 'Skin tone') %>%
  filter(model=='Basic controls & others GSPs') %>%
  rbind(y_mr2 %>%
          filter(model %in% c('Basic controls & others GSPs','Basic controls, others GSPs, & skin tone')) %>%
          rbind(y_skin_mr2) %>%
          mutate(outcome = 'Racial classification')) %>%
  mutate(outcome = factor(outcome, levels = c('Skin tone','Racial classification'))) %>%
  arrange(outcome, model, sample, term, marginal_r2) %>%
  select(outcome, model, sample, term, marginal_r2) %>%
  flextable() %>%
  set_header_labels(values = c('Outcome', 'Model', 'Sample', 'Term', 'Marginal R2')) %>%
  merge_v(j = c(1,3)) %>%
  merge_at(i = 1:8, j = 2) %>%
  merge_at(i = 9:16, j = 2) %>%
  merge_at(i = 17:24, j = 2) %>%
  hline(i = 8) %>%
  hline(i = 16, j = 2:5) %>%
  hline(i = c(4,12,21), j = 3:5) %>%
  colformat_md() %>%
  autofit()

save_as_docx(mr2_tbl, path = 'tables/appendix_gsp_mr2.docx')


# 2.5 Table 1 (appendix): top 5 countries that comprise each ancestry category (the fraction version)
top_countries_wt <- all$s4 %>%
  filter(`SUBREG America`==0, `SUBREG Canada`==0, `SUBREG Other`==0) %>%
  select(starts_with('CNTR '), -paste('CNTR', c('Other','America','Canada'))) %>%
  colSums(na.rm = T) %>%
  as.data.frame() %>%
  rownames_to_column()

names(top_countries_wt) <- c('country', 'n')

top_countries_wt <- top_countries_wt %>%
  mutate(country = str_remove(country, 'CNTR ')) %>%
  left_join(meta) %>%
  group_by(region, subregion) %>%
  arrange(subregion, desc(n)) %>%
  slice_head(n = 5) %>%
  transmute(region = region, subregion = subregion,
            country = paste0(country, ' (', round(n), ')'),
            ind = paste0('country', row_number())) %>%
  spread(ind, country) %>%
  ungroup() %>%
  rename(Region = region, `Sub-region` = subregion) %>%
  rename_with(capitalize, starts_with('country'))

tbl_s1 <- top_countries_wt %>%
  mutate(across(starts_with('country'),  ~ str_remove(., '\\*'))) %>%
  flextable() %>%
  add_footer_lines(value = as_paragraph(c('Other Middle East = Afghanistan, Armenia, Bahrain, Iran, Iraq, Israel, Jordan, Kuwait, Palestine, Oman, Qatar, Saudi Arabia, Syria, Turkey, United Arab Emirates, and Yemen',
                                          'Other East Europe = Albania, Bulgaria, Estonia, Latvia, Lithuania, Romania, and Yugoslavia',
                                          'Iceland and Atlantic Islands = Iceland, Azores, Bermuda, and Falkland Islands',
                                          'Pacific Islands, Australia, and Atlantic Islands = Figi, Guam, New Zealand, Tonga, Western Samoa, Australia, Azores, Bermuda, and Falkland Islands',
                                          'Other Southeast Asia = Borneo, Brunei, Cambodia, Indonesia, Laos, Malaysia, Mekong Valley, Myanmar, Singapore, and Thailand',
                                          'Other West Indies = Antigua and Barbuda, Aruba, Bahamas, Barbados, Cayman Islands, Dominica, Dominican Republic, Grenada, Guadeloupe, Montserrat, Netherlands Antilles, Saint Kitts and Nevis, Saint Lucia, Saint Vincent and the Grenadines, Trinidad and Tobago, and Virgin Islands of the United States',
                                          'Other Central America = Belize, Costa Rica, Guatemala, and Panama',
                                          'Other South America = Argentina, Bolivia, Brazil, Chile, Ecuador, French Guiana, Paraguay, Peru, Suriname, and Uruguay'))) %>%
  autofit()

save_as_docx(tbl_s1, path = 'tables/appendix_top5_countries.docx')


# 2.6 Table 2 (appendix): self-reported race vs geographic ancestry
d <- all$s4 %>%
  left_join(best_race) %>%
  mutate(best_race = factor(best_race, levels = c('NHW','Hispanic','NHB','AAPI','AINA')),
         region3 = factor(region3, levels = region_label)) # treat USA & Canada as missing

tbl2a <- table(d$best_race, d$region3, useNA = 'ifany') %>%
  as.data.frame() %>%
  filter(Var1!=' ') %>%
  group_by(Var1) %>%
  mutate(n = ifelse(is.na(Var2), NA, Freq),
         percent = n/sum(n, na.rm = T)*100,
         label = ifelse(is.na(Var2), Freq, paste0(Freq, ' \n(', sprintf('%.2f',percent),'%)'))) %>%
  select(` ` = Var1, Var2, label) %>%
  spread(Var2, label) %>%
  rename(Missing = `<NA>`)

tbl2b <- d %>%
  drop_na(best_race) %>%
  group_by(best_race) %>%
  summarise(across(starts_with('include_'), ~sum(., na.rm = T)),
            miss = sum(is.na(region3)), valid = n()-miss, miss = as.character(miss)) %>%
  mutate(across(starts_with('include_'), ~ paste0(.,'\n(', sprintf('%.2f',./valid*100),'%)'))) %>%
  select(best_race, paste0('include_',c('europe','america','ssa','eao','mesa')), miss) %>%
  arrange(best_race)
names(tbl2b) <- names(tbl2a)

tbl_s2 <- tbl2a %>%
  mutate(type = 'Single Ancestral Region') %>%
  rbind(tbl2b %>%
          mutate(type = 'All Ancestral Regions')) %>%
  mutate(type = ifelse(` `=='NHW', type, NA)) %>%
  relocate(type, .before = ` `) %>%
  flextable() %>%
  set_header_labels(type = '') %>%
  colformat_num(big.mark = '') %>%
  align(j = 3:8, align = 'center') %>%
  autofit()

save_as_docx(tbl_s2, path = 'tables/appendix_race_vs_ancestry.docx')


# 2.7 Table 3 (appendix): unsupervised K = 5 (the fraction version)
k5 <- all[c('us5', 's5')]
colnames(k5$us5)[1:5] <- paste0('us5_', colnames(k5$us5)[1:5])
colnames(k5$s5)[1:5] <- paste0('s5_', colnames(k5$s5)[1:5]) %>%
  str_replace('X1', 'Sub−Saharan African') %>%
  str_replace('X2', 'European') %>%
  str_replace('X3', 'East Asian') %>%
  str_replace('X4', 'Indigenous American') %>%
  str_replace('X5', 'Middle Eastern')

table_k5 <- k5$us5 %>%
  select(aid, starts_with('REGION '), starts_with('SUBREG '), starts_with('us5_')) %>%
  left_join(k5$s5 %>% select(aid, starts_with('s5_')))

tbl_us5 <- summarytable2(table_k5, paste0('us5_X', 1:5))
tbl_s5 <- summarytable2(table_k5, paste0('s5_', c(ancestry_label,'Middle Eastern')))

cor_k5 <- all[c('us5', 's5')] %>%
  map(~{
    .x <- .x %>%
      select(aid, paste0('X',1:5))
  }) %>%
  bind_rows(.id = 'admixture') %>%
  gather(ancestry, X, X1:X5) %>%
  spread(admixture, X) %>%
  group_by(ancestry) %>%
  summarise(cor = sprintf('%.3f', cor(s5, us5))) %>%
  ungroup() %>%
  spread(ancestry, cor)

panelA_us5 <- tbl_us5 %>%
  mutate(Region = case_when(Subregion=='Total' ~ Region, T ~ ''),
         `N = 8055` = paste0(sprintf('%.0f',`N = 8055`), ' (',
                             sprintf('%.0f', prop*100), '%)')) %>%
  rename_with(function(x) x = str_replace(x, 'us5', 'Unsupervised 5'), starts_with('us5')) %>%
  select(-prop) %>%
  rename(`Sub-Region` = Subregion) %>%
  flextable() %>%
  split_header() %>%
  merge_h(part = 'header') %>%
  set_header_labels(values = c(rep('',3), '*P*^1^','*P*^2^','*P*^3^',
                               '*P*^4^','*P*^5^')) %>%
  align(align = 'center', part = 'header') %>%
  colformat_md(part = 'header') %>%
  align(j = 4:8, align = 'center') %>%
  autofit()

panelB_s5 <- tbl_s5 %>%
  mutate(Region = case_when(Subregion=='Total' ~ Region, T ~ ''),
         `N = 8055` = paste0(sprintf('%.0f',`N = 8055`), ' (',
                             sprintf('%.0f', prop*100), '%)')) %>%
  rename_with(function(x) x = str_remove(x, 's5_'), starts_with('s5_')) %>%
  select(-prop) %>%
  rename(`Sub-Region` = Subregion) %>%
  flextable() %>%
  add_header_row(values = c('','*P*^AFR^','*P*^EUR^','*P*^EAS^','*P*^IAM^','*P*^MEA^'),
                 colwidths = c(3,rep(1,5))) %>%
  add_header_row(values = c('', 'Supervised 5'), colwidths = c(3,5)) %>%
  align(align = 'center', part = 'header') %>%
  colformat_md(part = 'header') %>%
  align(j = 4:8, align = 'center') %>%
  autofit()

panelC_k5 <- cor_k5 %>%
  transmute(` ` = '', `  ` = '', `   ` = '', X1, X2, X3, X4, X5) %>%
  flextable() %>%
  set_header_labels(values = c('','','','*P*^1^,*P*^AFR^','*P*^2^,*P*^EUR^','*P*^3^,*P*^EAS^',
                               '*P*^4^,*P*^IAM^','*P*^5^,*P*^MEA^')) %>%
  add_header_row(values = c('', 'Unsupervised-Supervised Correlation'), colwidths = c(3,5)) %>%
  align(align = 'center', part = 'header') %>%
  align(j = 4:8, align = 'center') %>%
  colformat_md(part = 'header') %>%
  autofit()

save_as_docx(panelA_us5, panelB_s5, panelC_k5, path = 'tables/appendix_unadx5_tbl.docx')


# 2.8 Table 4 (appendix): unsupervised K = 6 (the fraction version)
us6 <- all$us6
colnames(us6)[1:6] <- paste0('us6_', colnames(us6)[1:6])

tbl_us6 <- summarytable2(us6, paste0('us6_X', paste0(1:6)))

cor_k6 <- all[c('us6','s6')] %>%
  map(~{
    .x <- .x %>%
      select(aid, paste0('X',1:6))
  }) %>%
  bind_rows(.id = 'admixture') %>%
  gather(ancestry, X, X1:X6) %>%
  spread(admixture, X) %>%
  group_by(ancestry) %>%
  summarise(cor = sprintf('%.3f', cor(s6, us6))) %>%
  ungroup() %>%
  spread(ancestry, cor)

panelA_us6 <- tbl_us6 %>%
  mutate(Region = case_when(Subregion=='Total' ~ Region, T ~ ''),
         `N = 8055` = paste0(sprintf('%.0f',`N = 8055`), ' (',
                             sprintf('%.0f', prop*100), '%)')) %>%
  rename_with(function(x) x = str_replace(x, 'us6', 'Unsupervised 6'), starts_with('us6')) %>%
  select(-prop) %>%
  rename(`Sub-Region` = Subregion) %>%
  flextable() %>%
  split_header() %>%
  merge_h(part = 'header') %>%
  set_header_labels(values = c(rep('',3), '*P*^1^','*P*^2^','*P*^3^',
                               '*P*^4^','*P*^5^','*P*^6^')) %>%
  align(align = 'center', part = 'header') %>%
  colformat_md(part = 'header') %>%
  align(j = 4:9, align = 'center') %>%
  autofit()

panelB_us6 <- cor_k6 %>%
  transmute(` ` = '', `  ` = '', `   ` = '', X1, X2, X3, X4, X5, X6) %>%
  flextable() %>%
  set_header_labels(values = c('','','','*P*^1^,*P*^AFR^','*P*^2^,*P*^EUR^','*P*^3^,*P*^EAS^',
                               '*P*^4^,*P*^IAM^','*P*^5^,*P*^MEA^','*P*^6^,*P*^OCE^')) %>%
  add_header_row(values = c('', 'Unsupervised-Supervised Correlation'), colwidths = c(3,6)) %>%
  align(align = 'center', part = 'header') %>%
  align(j = 4:9, align = 'center') %>%
  colformat_md(part = 'header') %>%
  autofit()

save_as_docx(panelA_us6, panelB_us6, path = 'tables/appendix_unadx6_tbl.docx')



# 2.9 Table 6 (appendix): supervised K = 4 - percentiles (the fraction version)
tbl_s4p <- summarytable3(table_k4, paste0('s4_', ancestry_label)) %>%
  mutate(Region = case_when(Subregion=='Total' ~ Region, T ~ ''),
         `N = 8055` = paste0(sprintf('%.0f',`N = 8055`), ' (',
                             sprintf('%.0f', prop*100), '%)')) %>%
  rename_with(function(x) x = str_remove(x, 's4_'), starts_with('s4_')) %>%
  select(-prop) %>%
  rename(`Sub-Region` = Subregion) %>%
  flextable() %>%
  add_header_row(values = c('','*P*^AFR^','*P*^EUR^','*P*^EAS^','*P*^IAM^'),
                 colwidths = c(3,rep(1,4))) %>%
  add_header_row(values = c('', 'Supervised 4'), colwidths = c(3,4)) %>%
  align(align = 'center', part = 'header') %>%
  colformat_md(part = 'header') %>%
  align(j = 4:7, align = 'center') %>%
  autofit()

save_as_docx(tbl_s4p, path = 'tables/appendix_sadx4_percentiles.docx')


# 2.10 Table 8 (appendix): supervised K = 4 (the raw version)
table2_k4 <- k4$us4 %>%
  select(aid, subregion, region, starts_with('us4_')) %>%
  left_join(k4$s4 %>% select(aid, starts_with('s4_')))
by_subregion <- table2_k4 %>% rename(by_var = subregion)
by_region <- table2_k4 %>% rename(by_var = region)

# subregion-level & region-level statistics
tbl_s4_raw <- list(
  summarytable(by_subregion, paste0('s4_', ancestry_label), F, F, F, T, 'Subregion'),
  summarytable(by_region, paste0('s4_', ancestry_label), F, F, F, T, 'Region')
) %>%
  map(~{
    .x <- .x %>%
      modify_header(
        `s4_Sub−Saharan African` = 'Sub−Saharan African',
        s4_European = 'European',
        `s4_East Asian` = 'East Asian',
        `s4_Indigenous American` = 'Indigenous American'
      )
  }) %>%
  set_names(c('subregion', 'region'))

tbl_s4r <- tbl_s4_raw$subregion %>%
  modify_table_body(
    ~.x %>% 
      # remove undesired label rows
      filter(row_type!='label') %>%
      # append region rows
      rbind(tbl_s4_raw$region$table_body %>%
              filter(row_type!='label', label!='Sub-Saharan Africa')) %>%
      # create grouping rows
      left_join(meta %>% select(region, subregion) %>% unique(),
                by = c('label' = 'subregion')) %>%
      mutate(
        region = case_when(is.na(region) ~ label, T ~ region),
        groupname_col = factor(region, levels = c('USA & Canada',region_label)),
        label = ifelse(var_label=='Region'|label=='Sub-Saharan Africa', 'Total', label)
      ) %>%
      arrange(groupname_col, desc(label=='Total'), label)
  ) %>%
  as_flex_table() %>%
  set_header_labels(groupname_col = 'Subregion', label = 'Sub-Region') %>%
  add_header_row(values = c('','*P*^AFR^','*P*^EUR^','*P*^EAS^','*P*^IAM^'),
                 colwidths = c(3,rep(1,4))) %>%
  add_header_row(values = c('', 'Supervised 4'), colwidths = c(3,4)) %>%
  align(align = 'center', part = 'header') %>%
  colformat_md(part = 'header') %>%
  bold(bold = F, part = 'header') %>%
  align(j = 1:3, align = 'left') %>%
  autofit()

save_as_docx(tbl_s4r, path = 'tables/appendix_sadx4_raw.docx')


# 2.11 Table 9 (appendix): ICCs of the four GSPs by ancestry group
d <- all$s4 %>%
  mutate(region = factor(region3), subregion = factor(subregion2), country = factor(country2))
d <- d %>%
  select(starts_with('CNTR ')) %>%
  drop_na() %>%
  mutate_all(~ifelse(.>0, 1, .)) %>%
  transmute(cga_comb = do.call(paste, c(., sep = ' '))) %>%
  cbind(d %>% drop_na(starts_with('CNTR ')) %>% select(aid)) %>%
  right_join(d)

tbl_icc <- lapply(paste0('X',1:4), function(y) {
  lapply(c('region','subregion','country','cga_comb'), function(r) {
    lmer(as.formula(paste(y, '~ 1 |', r)), d) %>%
      performance::icc() %>%
      as.data.frame()
  }) %>%
    bind_rows() %>%
    mutate(geo_level = c('Region','Sub-region','Country','Country (multiple)'))
}) %>%
  set_names(paste0('X',1:4)) %>%
  bind_rows(.id = 'term') %>%
  transmute(Level = factor(geo_level, levels = c('Region','Sub-region','Country','Country (multiple)')),
            Term = case_when(term=='X1' ~ '*P*^AFR^', term=='X2' ~ '*P*^EUR^',
                             term=='X3' ~ '*P*^EAS^', term=='X4' ~ '*P*^IAM^'),
            ICC = sprintf('%.3f', ICC_unadjusted)) %>%
  arrange(Level) %>%
  flextable() %>%
  merge_v(j = 1) %>%
  hline(i = c(4,8,12)) %>%
  colformat_md() %>%
  autofit()

save_as_docx(tbl_icc, path = 'tables/appendix_icc_gsp_ancestry.docx')


# 2.12 Table 10 (appendix): racial classification by geographic ancestry
td <- all$s4 %>%
  select(aid, irace1 = irace, irace3, irace4, region3, starts_with('include_')) %>%
  pivot_longer(cols = matches('^irace'),
               names_to = c('.value', 'wave'),
               names_pattern = '^(irace)(\\d+)$') %>%
  filter(wave=='1')

tbl_rcga <- table(td$region3, td$irace, useNA = 'ifany') %>%
  as.data.frame() %>%
  drop_na(Var1) %>%
  group_by(Var1) %>%
  mutate(` ` = 'Self-Reported Geographic Ancestry',
         n = ifelse(is.na(Var2), NA, Freq),
         percent = n/sum(n, na.rm = T)*100,
         label = ifelse(is.na(Var2), Freq, paste0(Freq, ' \n(', sprintf('%.2f',percent),'%)'))) %>%
  select(` `, `  ` = Var1, Var2, label) %>%
  spread(Var2, label) %>%
  select(` `, `  `, Black, White, `Asian & Pacific Islander`, `Native American`, Other, Missing = `<NA>`) %>%
  flextable() %>%
  merge_v(j = 1) %>%
  rotate(j = 1, rotation = 'btlr') %>%
  add_header_row(values = c('','Interviewer-Classified Race'), colwidths = c(2,6)) %>%
  align(part = 'header', align = 'center') %>%
  align(j = 2:7, align = 'center') %>%
  autofit()

save_as_docx(tbl_rcga, path = 'tables/appendix_racialclass_vs_ancestry.docx')



# 3. Figures -------------------------------------------------------------------
# 3.1 Figure 2 (main text): the genetic ancestry of the US
# data preparation for river plot
# weighted using IPWs
fig2_dat <- all$s4 %>%
  select(ipwgt, paste0('X',1:4)) %>%
  filter(!is.na(ipwgt)) %>%
  mutate(across(X1:X4, function(x) x = x*ipwgt))

# average genetic similarity proportions at different Ks
fig2_prop <- fig2_dat %>%
  summarise(across(X1:X4, function(x) x = weighted.mean(x, ipwgt))) %>%
  gather(ID, prop) %>%
  mutate(prop = prop/sum(prop))

s3 <- fig2_prop %>%
  mutate(ID = ifelse(ID %in% c('X3','X4'), 'X3', ID)) %>%
  group_by(ID) %>%
  summarise(prop = sum(prop)) %>%
  ungroup()

s2 <- s3 %>%
  mutate(ID = ifelse(ID=='X1', ID, 'X2')) %>%
  group_by(ID) %>%
  summarise(prop = sum(prop)) %>%
  ungroup()

fig2_prop <- rbind(s2 %>% mutate(ID = paste0('s2_',ID)),
                   s3 %>% mutate(ID = paste0('s3_',ID)),
                   fig2_prop)

# prepare node data for the river plot
s4 <- data.frame(Value = colMeans(fig2_dat[,paste0('X',1:4)])) %>%
  rownames_to_column(var = 'ID') %>%
  mutate(x = 3,
         col = case_when(
           ID=='X1' ~ '#ffc300', ID=='X2' ~ '#219ebc', ID=='X3' ~ '#a7c957', 
           ID=='X4' ~ '#e76f51', T ~ NA
         ))
s4 <- s4 %>% cbind(t(col2rgb(s4$col)) %>% as.data.frame())
s3 <- s4 %>%
  mutate(x = 2, ID = case_when(ID %in% c('X3', 'X4') ~ 's3_X3', T ~ paste0('s3_', ID))) %>%
  group_by(ID) %>%
  mutate(prop = Value/sum(Value)) %>%
  summarise(Value = sum(Value), x = mean(x),
            red = sqrt(sum(prop*red^2)), green = sqrt(sum(prop*green^2)), blue = sqrt(sum(prop*blue^2)))%>%
  mutate(col = rgb(red, green, blue, maxColorValue = 255))
s2 <- s3 %>%
  mutate(x = 1, ID = case_when(ID %in% c('s3_X2', 's3_X3') ~ 's2_X2', T ~ ID %>% str_replace('s3','s2'))) %>%
  group_by(ID) %>%
  mutate(prop = Value/sum(Value)) %>%
  summarise(Value = sum(Value), x = mean(x),
            red = sqrt(sum(prop*red^2)), green = sqrt(sum(prop*green^2)), blue = sqrt(sum(prop*blue^2)))%>%
  mutate(col = rgb(red, green, blue, maxColorValue = 255))
s1 <- s2 %>%
  mutate(x = 0, ID = 's1') %>%
  group_by(ID) %>%
  mutate(prop = Value/sum(Value)) %>%
  summarise(Value = sum(Value), x = mean(x),
            red = sqrt(sum(prop*red^2)), green = sqrt(sum(prop*green^2)), blue = sqrt(sum(prop*blue^2)))%>%
  mutate(col = rgb(red, green, blue, maxColorValue = 255))

nodes <- rbind(s1,s2,s3,s4)
nodes$labels <- c('','Sub-Saharan African', 'European, \nEast Asian, & Indigenous American',
                  'Sub-Saharan African', 'European', 'East Asian &\nIndigenous American',
                  'Sub-Saharan African', 'European', 'East Asian', 'Indigenous American')
nodes <- nodes %>%
  left_join(fig2_prop) %>%
  mutate(prop =ifelse(is.na(prop), '', paste0('(',sprintf('%.2f',prop*100),'%)')),
         labels = paste0(labels, '\n', prop))

# prepare edge data for the river plot
edges <- data.frame(
  N1 = c('s1', 's1', 's2_X1', 's2_X2', 's2_X2', 's3_X1', 's3_X2', 's3_X3', 's3_X3'),
  N2 = c('s2_X1', 's2_X2', 's3_X1', 's3_X2', 's3_X3', 'X1', 'X2', 'X3', 'X4')) %>%
  left_join(nodes %>% select(N2 = ID, Value, col))

# combine and convert node and edge data into the riverplot format
riv_data <- makeRiver(nodes %>% select(-Value) %>% as.data.frame(), edges)

# correlation at the splitting branch
rho3 <- cor(all$s2$X2, all$s3$X2 + all$s3$X3)
rho4 <- cor(all$s3$X3, all$s4$X3 + all$s4$X4)

# left panel of figure 2: river plot
plot(riv_data, plot_area = c(0.9,0.9), srt = 0, textcex = 0.8, node_margin = 0.3)
text(0.36, -0.095, expression(italic('K')==2), xpd = NA)
text(0.64, -0.095, expression(italic('K')==3), xpd = NA)
text(0.64, -0.15, bquote(paste(rho['2,3'], ' = ', .(sprintf('%.2f', rho3)))), xpd = NA)
text(0.92, -0.095, expression(italic('K')==4), xpd = NA)
text(0.92, -0.15, bquote(paste(rho['3,4'], ' = ', .(sprintf('%.2f', rho4)))), xpd = NA)
text(1.05, 0.94, expression(paste('(', italic('P')^IAM, ')')), cex = 0.8, xpd = NA)
text(1.05, 0.8, expression(paste('(', italic('P')^EAS, ')')), cex = 0.8, xpd = NA)
text(1.05, 0.45, expression(paste('(', italic('P')^EUR, ')')), cex = 0.8, xpd = NA)
text(1.05, 0.09, expression(paste('(', italic('P')^AFR, ')')), cex = 0.8, xpd = NA)
segments(-0.05, -0.05, 1.05, -0.05, xpd = NA)
riverplot <- recordPlot()

# right panel of figure 2: cumulative probability
fig2_right <- all$s4 %>%
  select(aid, paste0('X',1:4), ipwgt) %>%
  filter(!is.na(ipwgt)) %>%
  gather('anc','X',X1:X4) %>%
  group_by(anc) %>%
  arrange(anc, X) %>%
  mutate(n = cumsum(ipwgt),
         anc = case_when(anc=='X1' ~ 'AFR',
                         anc=='X2' ~ 'EUR',
                         anc=='X3' ~ 'EAS',
                         anc=='X4' ~ 'IAM', 
                         T ~ NA) %>% factor(levels = ancestry_label2),
         frac = n/max(n)) %>%
  ggplot() +
  geom_line(aes(x = X, y = frac, group = anc, color = anc), linewidth = 1) + 
  geom_hline(yintercept = seq(0, 0.75, 0.25), color = 'gray', linetype = 'dotted', linewidth = 0.3) +
  geom_vline(xintercept = seq(0, 0.75, 0.25), color = 'gray', linetype = 'dotted', linewidth = 0.3) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.25), limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.25), limits = c(0,1)) +
  labs(x = bquote(italic(P)), y = 'Cumulative Probability', color = '') +
  scale_color_manual(values = c(AFR= '#ffc300', EUR = '#219ebc', EAS = '#a7c957',
                                IAM = '#e76f51'),
                     labels = c(AFR = bquote(italic(P)^'AFR'), EUR = bquote(italic(P)^'EUR'),
                                EAS = bquote(italic(P)^'EAS'), IAM = bquote(italic(P)^'IAM'))) +
  theme_custom() +
  theme(legend.position=c(0.9,0.15), legend.background = element_rect(fill='transparent'))

ggsave(paste0('figures/figure2',fig_format), width = 18, height = 6,
       plot_grid(as_grob(riverplot), as_grob(fig2_right), align = 'hv', nrow = 1, rel_heights = c(1,3)))


# 3.2 Figure 3 (main text): skin color scatter plots
s4 <- all$s4 %>% mutate(region = region3 %>% factor())
d <- s4 %>%
  filter(include_ssa==1) %>%
  mutate(region = 'Sub-Saharan Africa') %>%
  rbind(s4 %>%
          filter(include_america==1) %>%
          mutate(region = 'The Americas')) %>%
  rbind(s4 %>%
          filter(include_europe==1) %>%
          mutate(region = 'Europe')) %>%
  rbind(s4 %>%
          filter(include_eao==1) %>%
          mutate(region = 'Eastern Asia & Oceania')) %>%
  rbind(s4 %>%
          filter(include_mesa==1) %>%
          mutate(region = 'Middle East & South Asia')) %>%
  mutate(region = factor(region, levels = region_label[c(3,1,4,2,5)]))

prepare_skin <- function(data, X2_thresh, X0_thresh, both_thresh, discrete = T) {
  lm_dat <- data %>% filter(!is.na(region), !is.na(skin))
  # GSP-related inclusion criteria
  lm_dat <- lm_dat %>% filter(X2>=X2_thresh, X0>=X0_thresh, X2+X0>=both_thresh)
  # geographic ancestry-related inclusion criteria
  lm_dat <- lm_dat %>%
    group_by(region) %>%
    mutate(regn = max(row_number())) %>%
    filter(regn>=110) %>%
    ungroup()
  
  if (discrete) {
    dat <- lm_dat %>%
      transmute(aid, region, include_ssa, include_america, include_europe,
                sex, X0, skin)
    dat <- dat %>%
      cbind(predict(dummyVars('~ skin', dat, sep = '_'), newdata = dat) %>% as.data.frame()) %>%
      gather(ytype, y, starts_with('skin_')) %>% 
      mutate(ytype = str_remove_all(ytype, 'skin_'))
  } else {
    lm_dat <- lm_dat %>% mutate(skin = 6 - as.numeric(skin))
    # modeling & residualization
    Xs <- names(data %>% select(starts_with('X')))
    lXs <- Xs[!Xs %in% c('X0','X2')] %>% paste(collapse = ' + ')
    mx <- lm(paste('X0 ~', lXs), lm_dat)
    my <- lm(paste('skin ~', lXs), lm_dat)
    dat <- lm_dat %>% 
      transmute(aid, region, include_ssa, include_america, include_europe,
                sex, raw_x = X0, res_x = residuals(mx) + mean(lm_dat$X0), res2_x = X0,
                raw_y = skin, res_y = residuals(my) + mean(lm_dat$skin), res2_y = res_y)
  }
  
  return(dat)
}

my_skin_plot <- function(discrete = F) {
  if (discrete) {
    fig3_data1 <- prepare_skin(d %>% rename(X0 = X1), 0.01, 0.01, -1, T) %>%
      mutate(anc = 'AFR') %>%
      rbind(prepare_skin(d %>% rename(X0 = X4), 0.01, 0.01, -1, T) %>%
              mutate(anc = 'IAM')) %>%
      mutate(anc = factor(anc, levels = ancestry_label2[-2]), x = X0,
             ytype = factor(ytype, levels = levels(d$skin)),
             # whether include in the main text figure
             main = case_when(
               anc=='AFR' & region %in% c('Sub-Saharan Africa', 'The Americas') ~ 1,
               anc=='IAM' & region=='The Americas' ~ 1,
               T ~ 0
             )) %>%
      filter(ytype!='White', region %in% c('Sub-Saharan Africa', 'The Americas'))
    
    texty <- 0.9
    rugy <- 1.01
    rugyend <- 1.05
    yintercept <- seq(0,1,0.25)
    ylabel <- 'Fraction with a given skin tone (or darker)'
    ylimit <- c(0,1.05)
    ybreaks <- seq(0.25,1,0.25)
    linetype <- 'Interviewer-Rated Skin Tone'
  } else {
    fig3_data1 <- prepare_skin(d %>% filter(region=='Sub-Saharan Africa') %>% 
                                 rename(X0 = X1), 0.01, 0.01, -1, F) %>%
      mutate(anc = 'AFR') %>%
      rbind(prepare_skin(d %>% filter(region=='The Americas') %>%
                           rename(X0 = X1), 0.01, 0.01, -1, F) %>%
              mutate(anc = 'AFR')) %>%
      rbind(prepare_skin(d %>% filter(region=='The Americas') %>%
                           rename(X0 = X4), 0.01, 0.01, -1, F) %>%
              mutate(anc = 'IAM')) %>%
      pivot_longer(
        cols = ends_with(c('_x', '_y')),
        names_to = c('ytype', '.value'),
        names_pattern = '(.*)_(.*)'
      ) %>%
      mutate(main = 1, anc = factor(anc, levels = ancestry_label2[-2]),
             ytype = factor(ytype, levels = c('raw','res','res2'), 
                            labels = c('Raw','Residualized','Residualized2')))
    
    texty <- 4.5
    rugy <- 5.1
    rugyend <- 5.5
    yintercept <- seq(0,5,0.5)
    ylabel <- 'Skin tone'
    ylimit <- c(0,5.5)
    ybreaks <- seq(0.5,5,0.5)
    linetype <- 'Version'
  }
  
  # data for region-specific sample size & significance  
  fig3_dat1 <- fig3_data1 %>% 
    group_by(anc, region, main, ytype) %>%
    summarise(n = max(row_number())) %>%
    group_by(anc, ytype) %>%
    mutate(
      x = case_when(region=='Sub-Saharan Africa' & anc!='AFR' ~ 0.78, 
                    discrete ~ 0.18,
                    T ~ 0.02),
      y = texty,
      n = paste0('Self-Reported Ancestral',
                 ifelse(discrete,'\n',' '),
                 ifelse(region=='Sub-Saharan Africa','Sub-Saharan Africans','Americans'),
                 '\n(N = ', n, ')'),
      rowno = ifelse(region=='The Americas',2,1)
    ) %>%
    ungroup() %>%
    select(-ytype) %>%
    unique()
  fig3_data1 <- fig3_data1 %>% left_join(fig3_dat1 %>% select(anc, region, rowno))
  
  # data for plotting binned dots
  fig3_dat2 <- fig3_data1 %>%
    group_by(anc, region, rowno, ytype) %>%
    arrange(anc, region, rowno, ytype, x) %>%
    group_modify(~{
      n <- round(nrow(.x)/80)
      breaks <- quantile(.x$x, probs = seq(0,1,length.out=n+1))
      .x$bin <- cut(.x$x, breaks = breaks, labels = 1:n)
      .x$bin <- ifelse(is.na(.x$bin),1,.x$bin)
      return(.x)
    }) %>%
    group_by(anc, region, main, rowno, ytype, bin) %>%
    summarise(across(c(x,y), function(v) v = mean(v))) %>%
    ungroup()
  
  if (discrete) {
    fig3_dat2 <- fig3_dat2 %>%
      group_by(anc, region, rowno, bin) %>%
      mutate(y = cumsum(y))
  }
  
  # data for loess lines using unbinned points
  fig3_dat3 <- fig3_data1 %>%
    left_join(fig3_dat2 %>% group_by(anc, region) %>% summarise(x_end = max(x)))
  
  if (discrete) {
    fig3_dat3 <- fig3_dat3 %>%
      group_by(aid, anc, region, rowno) %>%
      arrange(aid, anc, region, rowno, ytype) %>%
      mutate(y = cumsum(y))
  }
  
  fig3_dat3 <- fig3_dat3 %>%
    group_by(anc, region, rowno, ytype) %>%
    group_modify(~{
      # fit loess lines
      ## degree: the degree of the polynomials to be used; span: the degree of smoothing
      m <- loess(y ~ x, data = .x, span = 0.8, degree = 1)
      # get fitted y values
      .x$hat_y <- predict(m)
      .x <- .x %>% filter(x<=x_end)
      return(.x)
    }) %>% 
    ungroup()
  
  if (discrete) {
    ggplot() +
      geom_point(data = fig3_dat2, show.legend = F,
                 aes(x = x, y = y, color = region, shape = ytype), size = 2, alpha = 0.4) +
      geom_line(data = fig3_dat3, 
                aes(x = x, y = hat_y, color = region, linetype = ytype), linewidth = 1) +
      geom_text(data = fig3_dat1, 
                aes(x = x, y = y, label = n, color = region), show.legend = F, size = 3.5) +
      geom_segment(data = fig3_data1 %>% select(-ytype) %>% unique(), 
                   aes(x = x, xend = x, y = 1.01, yend = 1.05, color = region), alpha = 0.1) + 
      geom_hline(yintercept = seq(0,1,0.25), color = 'gray', linetype = 'dotted', linewidth = 0.3) +
      geom_vline(xintercept = seq(0,1,0.25), color = 'gray', linetype = 'dotted', linewidth = 0.3) +
      scale_x_continuous(name = '', limit = c(0,1), breaks = seq(0.25,1,0.25), expand = c(0,0)) +
      scale_y_continuous(name = 'Fraction with a given skin tone (or darker)', limit = c(0,1.05), breaks = seq(0.25,1,0.25), expand = c(0,0)) +
      facet_grid(region ~ anc, switch = 'x', labeller = label_bquote(cols = bold(italic(P)^.(as.character(anc))))) +
      scale_color_manual(name = 'Self-Reported Geographic Ancestry', 
                         values = c(`The Americas`='#dd0426', `Sub-Saharan Africa`='#f9a620'),
                         labels = c(`The Americas`='Ancestral Americans', 
                                    `Sub-Saharan Africa`='Ancestral Sub-Saharan Africans')) +
      scale_shape_manual(name = 'Interviewer-Rated Skin Tone', values = c(0,1,2,5)) +
      scale_linetype_manual(name = 'Interviewer-Rated Skin Tone', values = c(1,3,2,6)) +
      theme_custom() +
      theme(legend.position = 'top', legend.box = 'vertical', legend.spacing.y = unit(0.05,'cm'),
            strip.text.y = element_blank(), strip.placement = 'outside')
  } else {
    ggplot() +
      geom_point(data = fig3_dat2 %>% filter(main==1, ytype=='Residualized'), show.legend = F,
                 aes(x = x, y = y, color = region), shape = 15, size = 2) +
      geom_point(data = fig3_data1 %>% filter(main==1, ytype=='Residualized'), show.legend = F,
                 aes(x = x, y = y, color = region), shape = 1, size = 2, alpha = 0.1) +
      geom_line(data = fig3_dat3 %>% filter(main==1, ytype=='Residualized'),
                aes(x = x, y = hat_y, color = region), linewidth = 1) +
      geom_text(data = fig3_dat1 %>% filter(main==1), 
                aes(x = x, y = y, label = n, color = region), show.legend = F, hjust = 0, size = 3.5) +
      geom_segment(data = fig3_data1 %>% filter(main==1) %>% select(-ytype) %>% unique(), 
                   aes(x = x, xend = x, y = rugy, yend = rugyend, color = region), alpha = 0.1) + 
      geom_hline(yintercept = yintercept, color = 'gray', linetype = 'dotted', linewidth = 0.3) +
      geom_vline(xintercept = seq(0,1,0.25), color = 'gray', linetype = 'dotted', linewidth = 0.3) +
      scale_x_continuous(name = '', limit = c(0,1), breaks = seq(0.25,1,0.25), expand = c(0,0)) +
      scale_y_continuous(name = ylabel, limit = ylimit, breaks = ybreaks, expand = c(0,0)) +
      facet_wrap(rowno ~ anc, scales = 'free', strip.position = 'bottom', 
                 labeller = label_bquote(cols = bold(italic(P)^.(as.character(anc))))) +
      scale_color_manual(name = 'Self-Reported Geographic Ancestry', 
                         values = c(`The Americas`='#dd0426', `Sub-Saharan Africa`='#f9a620'),
                         labels = c(`The Americas`='Ancestral Americans', 
                                    `Sub-Saharan Africa`='Ancestral Sub-Saharan Africans')) +
      theme_custom() +
      theme(legend.position = 'top', legend.box = 'vertical', legend.spacing.y = unit(0.05,'cm'),
            strip.text.y = element_blank(), strip.placement = 'outside')
  }
}

# Figure 3 (main text): skin tone as continuous
my_skin_plot(F)
ggsave(paste0('figures/figure3',fig_format), width = 12, height = 5)

# Figure 17 (appendix): skin tone as discrete
my_skin_plot(T)
ggsave(paste0('figures/appendix_ordinal_figure3',fig_format), width = 10, height = 8)


# 3.3 Figure 4 (main text): self-reported, interviewer-classified race, & skin color
d <- d %>%
  gather(var, irace, irace:irace4) %>%
  filter(region=='Sub-Saharan Africa' | (region=='The Americas' & var=='irace'), 
         !is.na(irace), !is.na(skin)) %>%
  mutate(y = ifelse(region=='Sub-Saharan Africa', as.numeric(irace=='Black'), 
                    1 - as.numeric(irace=='White')))

# racial classification by skin tone category
prepare_irace <- function(data, ga_sample) {
  if (ga_sample=='Sub-Saharan Africa') {
    y_label <- 'Fraction classified as Black'
  } else if (ga_sample=='The Americas') {
    y_label <- 'Fraction classified as non-White'
  } else {
    print('No such ancestral region category.')
  }
  
  data <- data %>% filter(region==ga_sample)
  
  if (any(duplicated(data$aid))) {
    dat <- summary(lm_robust(y ~ skin - 1, cluster = aid, data))$coefficients
  } else {
    dat <- summary(lm(y ~ skin - 1, data))$coefficients %>%
      cbind(confint(lm(y ~ skin - 1, data)))
    colnames(dat)[grepl('%', colnames(dat))] <- c('CI Lower', 'CI Upper')
  }
  
  rd <- dat %>%
    as.data.frame() %>%
    rownames_to_column('skin') %>%
    mutate(skin = str_remove(skin, 'skin')) %>%
    left_join(data %>% select(aid, skin) %>% unique() %>% group_by(skin) %>% summarise(n = max(row_number()))) %>%
    ungroup() %>%
    transmute(skin, skin_prop = n/sum(n), skin_label = paste0(skin,' (',round(skin_prop*100),'%)'),
              y = Estimate, se = `Std. Error`, lower = `CI Lower`, upper = `CI Upper`,
              y_label = y_label, sample = ga_sample, 
              sample_label = case_when(sample=='Sub-Saharan Africa' ~ 'Sub-Saharan Africans',
                                       sample=='The Americas' ~ 'Americans'),
              sample_label = paste0('Self−Reported Ancestral\n', sample_label, '\n(N = ',sum(n),')'))
  
  return(rd)
}

pd <- lapply(c('Sub-Saharan Africa','The Americas'), function(ga) prepare_irace(d, ga)) %>%
  bind_rows()
pd <- pd %>%
  mutate(skin_label = factor(skin_label, levels = rev(factor(pd$skin_label) %>% levels())[c(1:2,5:6,3:4,7:10)]))

# racial classification by skin tone category and genetic similarity proportion 
prepare_irace_bygsp <- function(data, ga_sample, focal_gsp) {
  
  if (ga_sample=='Sub-Saharan Africa') {
    y_label <- 'Fraction classified as Black'
  } else if (ga_sample=='The Americas') {
    y_label <- 'Fraction classified as non-White'
  } else {
    print('No such ancestral region category.')
  }
  
  # residualize the focal genetic similarity proportion on the others
  data <- data %>% filter(region==ga_sample)
  covX <- setdiff(paste0('X',c(1,3:4)), focal_gsp)
  unique_data <- data %>% select(aid, paste0('X',1:4)) %>% unique()
  mx <- lm(as.formula(paste(focal_gsp, '~', paste(covX, collapse = ' + '))), unique_data)
  unique_data$res_X0 <- mx$residuals + mean(unique_data[[focal_gsp]])
  data <- data %>% left_join(unique_data %>% select(aid, res_X0))
  my <- lm(as.formula(paste('y ~', paste(covX, collapse = ' + '))), data)
  data$res_y <- my$residuals + mean(data$y)
  
  # split each skin tone category into two equally size subgroups
  data <- data %>%
    arrange(skin, res_X0) %>%
    group_by(skin) %>%
    mutate(lmedian = ifelse(row_number()<=max(row_number())/2,1,0)) %>%
    ungroup()
  
  if (any(duplicated(data$aid))) {
    dat <- data %>%
      group_by(lmedian) %>%
      group_modify(~{
        m <- lm_robust(res_y ~ skin - 1, cluster = aid, .x)
        .x <- summary(m)$coefficients %>%
          as.data.frame() %>% 
          rownames_to_column('skin') %>%
          mutate(df.residual = m$df.residual)
      })
  } else {
    dat <- data %>%
      group_by(lmedian) %>%
      group_modify(~{
        m <- lm(res_y ~ skin - 1, .x)
        rd <- summary(m)$coefficients %>%
          cbind(confint(m))
        colnames(rd)[grepl('%', colnames(rd))] <- c('CI Lower', 'CI Upper')
        rd <- as.data.frame(rd) %>% 
          rownames_to_column('skin') %>%
          mutate(df.residual = m$df.residual)
        return(rd)
      })
  }
  
  dat <- dat %>%
    transmute(skin = str_remove(skin, 'skin'), y = Estimate, se = `Std. Error`,
              lower = `CI Lower`, upper = `CI Upper`, df.residual) %>%
    left_join(data %>% select(aid, lmedian, skin) %>% unique() %>%
                group_by(lmedian, skin) %>% summarise(n = max(row_number()))) %>%
    ungroup()
  
  # t test
  dat2 <- dat %>%
    pivot_wider(id_cols = skin,
                names_from = lmedian,
                values_from = c(y, se, df.residual),
                names_sep = '') %>%
    mutate(skin, diff_coef = y0 - y1, se_diff = sqrt(se0^2 + se1^2),
           t_stat = diff_coef/se_diff) %>%
    rowwise() %>%
    mutate(pvalue = 2*pt(abs(t_stat), df = min(df.residual0, df.residual1), lower.tail = F))
  
  rd <- dat %>%
    left_join(dat2 %>% select(skin, pvalue)) %>%
    group_by(lmedian) %>%
    transmute(skin, skin_prop = n/sum(n), skin_label = paste0(skin,' (',round(skin_prop*100),'%)'),
              y, se, upper, lower, y_label = y_label, sample = ga_sample, 
              sample_label = case_when(sample=='Sub-Saharan Africa' ~ 'Sub-Saharan Africans',
                                       sample=='The Americas' ~ 'Americans'),
              sample_label = paste0('Self−Reported Ancestral\n', sample_label, '\n(N = ',sum(n),')'), 
              pvalue = case_when(
                pvalue<0.001 ~ '***',
                pvalue>=0.001 & pvalue<0.01 ~ '**',
                pvalue>=0.01 & pvalue<0.05 ~ '*',
                pvalue>=0.05 & pvalue<0.1 ~ '.',
                T ~ ''
              )) %>%
    ungroup()
  
  return(rd)
}

pd2 <- prepare_irace_bygsp(d, 'Sub-Saharan Africa', 'X1') %>%
  mutate(panel = 1) %>%
  rbind(
    prepare_irace_bygsp(d, 'The Americas', 'X4') %>%
      mutate(panel = 2)
  ) %>%
  rbind(
    prepare_irace_bygsp(d, 'The Americas', 'X1') %>%
      mutate(panel = 3)
  ) %>%
  left_join(pd %>%
              transmute(skin_label2 = skin_label, skin, y_label, sample_label2 = sample_label)) %>%
  mutate(across(c(y,se,lower,upper), ~ ifelse(parse_number(as.character(skin_label2))<=3, NA, .)))

pd3 <- pd2 %>%
  group_by(panel, sample, skin_label2) %>%
  mutate(y = max(y) + 0.01) %>%
  transmute(panel, sample, skin_label2, y, pvalue, lmedian = 0) %>%
  unique() %>%
  filter(!pvalue %in% c('.',''))

my_irace_plot <- function(ga_sample, panel_num) {
  if (ga_sample=='Sub-Saharan Africa') {
    y_label <- 'Fraction classified as Black'
  } else if (ga_sample=='The Americas') {
    y_label <- 'Fraction classified as non-White'
  } else {
    print('No such ancestral region category.')
  }
  
  if (ga_sample=='Sub-Saharan Africa' | (ga_sample=='The Americas' & panel_num==3)) {
    color <- '#ffc300'
    anc <- 'AFR'
  } else if (ga_sample=='The Americas' & panel_num==2) {
    color <- '#e76f51'
    anc <- 'IAM'
  } else {
    print('No such specification available.')
  }
  
  p <- pd2 %>%
    filter(sample==ga_sample, panel==panel_num) %>%
    ggplot(aes(x = skin_label2, y = y, fill = as.factor(1-lmedian))) +
    geom_bar(stat = 'identity', position = 'dodge') +
    geom_text(aes(y = lower, label = sprintf('%.2f',y)), position = position_dodge(0.9), vjust = 1.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(0.9)) +
    geom_text(data = pd2 %>% filter(sample==ga_sample, panel==panel_num, grepl('White', skin)),
              aes(y = 1, label = sample_label2), size = 3.6) +
    geom_text(data = pd3 %>% filter(sample==ga_sample, panel==panel_num),
              aes(label = pvalue), color = '#bf0603', show.legend = FALSE) +
    scale_y_continuous(name = y_label, expand = c(0,0), breaks = seq(0,1,0.25)) +
    coord_cartesian(ylim = c(0, 1.2)) + 
    scale_x_discrete(name = '') +
    scale_fill_manual(name = bquote(italic('P')^.(anc)),
                      values = c(`0`= rgb(col2rgb(color)[1]/255, col2rgb(color)[2]/225, col2rgb(color)[3]/225, alpha = 0.4), `1`=color),
                      labels = c(`0`='< median', `1`='>=median')) +
    theme_custom() +
    theme(legend.position = 'top')
  
  return(p)
}

ggsave('figures/figure4.pdf', width = 8, height = 12,
       plot_grid(my_irace_plot('Sub-Saharan Africa',1), my_irace_plot('The Americas',2), 
                 my_irace_plot('The Americas',3), nrow = 3))

# Figure 19 (appendix): racial classification by skin tone
pd %>%
  ggplot(aes(x = skin_label, y = y)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_text(data = pd %>% filter(grepl('White', skin)), 
            aes(y = 1, label = sample_label), size = 3.6) +
  facet_wrap(~ sample, nrow = 2, strip.position = 'left', scales='free',
             labeller = as_labeller(c(`Sub-Saharan Africa` = 'Fraction classified as Black',
                                      `The Americas` = 'Fraction classified as non-White'))) +
  scale_y_continuous(name = '', expand = c(0,0), breaks = seq(0,1,0.25), limits = c(0,1.2)) +
  scale_x_discrete(name = '') +
  theme_custom() +
  theme(strip.background = element_blank(),
        strip.placement = 'outside')
ggsave('figures/appendix_nogsp_figure4.pdf', width = 8, height = 8)


# 3.4 Figure 2 (appendix): branching structure of unsupervised K = 2 ~ 6
us <- all[paste0('us',2:6)]
for (k in 2:6) {
  colnames(us[[paste0('us',k)]])[1:k] <- paste0('us', k, '_', colnames(us[[paste0('us',k)]])[1:k])
}

dat_us <- us$us2 %>%
  select(aid, region3, starts_with('us2_')) %>%
  left_join(us$us3 %>% select(aid, starts_with('us3_'))) %>%
  left_join(us$us4 %>% select(aid, starts_with('us4_'))) %>%
  left_join(us$us5 %>% select(aid, starts_with('us5_'))) %>%
  left_join(us$us6 %>% select(aid, starts_with('us6_'))) %>%
  mutate(sum1 = us3_X2 + us3_X3, sum2 = us4_X3 + us4_X4, sum3 = us5_X2 + us5_X5, sum4 = us6_X3 + us6_X6, 
         color = case_when(is.na(region3)|region3=='USA & Canada' ~ 'USA, Canada, & Unknown',
                           T ~ region3),
         color = factor(color, levels = c(region_label, 'USA, Canada, & Unknown'))) 

p1 <- my_scatter_rug_plot(dat_us %>% rename(x = us2_X2, y = sum1), expression(atop(italic(P)^2,(italic(K)==2))),
                          expression(atop(italic(P)^2 + italic(P)^3,(italic(K)==3))), F)
p2 <- my_scatter_rug_plot(dat_us %>% rename(x = us3_X3, y = sum2), expression(atop(italic(P)^3,(italic(K)==3))),
                          expression(atop(italic(P)^3 + italic(P)^4,(italic(K)==4))), F)
p3 <- my_scatter_rug_plot(dat_us %>% rename(x = us4_X2, y = sum3), expression(atop(italic(P)^2,(italic(K)==4))),
                          expression(atop(italic(P)^2 + italic(P)^5,(italic(K)==5))), F)
p4 <- my_scatter_rug_plot(dat_us %>% rename(x = us5_X3, y = sum4), expression(atop(italic(P)^3,(italic(K)==5))),
                          expression(atop(italic(P)^3 + italic(P)^6,(italic(K)==6))), T)
ggsave(paste0('figures/appendix_usadmx_branching',fig_format), plot_grid(p1,p2,p3,p4, ncol = 2),
       width = 12, height = 10)


# 3.5 Figure 3 (appendix): comparing supervised and unsupervised K = 4 estimates
k4 <- all[c('us4', 's4')]
colnames(k4$us4)[1:4] <- paste0('us4_', colnames(k4$us4)[1:4])
colnames(k4$s4)[1:4] <- paste0('s4_', colnames(k4$s4)[1:4])

dat_k4 <- k4$us4 %>%
  select(aid, region3, starts_with('us4_')) %>%
  left_join(k4$s4 %>% select(aid, starts_with('s4_'))) %>%
  mutate(color = case_when(is.na(region3)|region3=='USA & Canada' ~ 'USA, Canada, & Unknown',
                           T ~ region3),
         color = factor(color, levels = c(region_label, 'USA, Canada, & Unknown')))

p1 <- my_scatter_rug_plot(dat_k4 %>% select(x = s4_X1, y = us4_X1, color),
                          expression(italic(P)^AFR), expression(italic(P)^1), F)
p2 <- my_scatter_rug_plot(dat_k4 %>% select(x = s4_X2, y = us4_X2, color),
                          expression(italic(P)^EUR), expression(italic(P)^2), F)
p3 <- my_scatter_rug_plot(dat_k4 %>% select(x = s4_X3, y = us4_X3, color),
                          expression(italic(P)^EAS), expression(italic(P)^3), F)
p4 <- my_scatter_rug_plot(dat_k4 %>% select(x = s4_X4, y = us4_X4, color),
                          expression(italic(P)^IAM), expression(italic(P)^4), T)
ggsave(paste0('figures/appendix_compare_svsus',fig_format), plot_grid(p1,p2,p3,p4, nrow = 2),
       width = 12, height = 10)


# 3.6 Figure 5 (appendix): distribution of the 4 GSP variables in the US
s4 %>%
  gather(anc, X, paste0('X',1:4)) %>%
  mutate(anc = case_when(
    anc=='X1' ~ 'AFR', anc=='X2' ~ 'EUR', anc=='X3' ~ 'EAS', anc=='X4' ~ 'IAM',
    T ~ NA
  ) %>% factor(levels = c('AFR','EAS','IAM','EUR'))) %>%
  ggplot() +
  geom_histogram(aes(x = X, y = 0.05*..density.., fill = anc, weight = ipwgt),
                 binwidth = 0.05, color = 'gray60', center = 0.025) +
  facet_wrap(~ anc, nrow = 2, scale = 'free', strip.position = 'bottom', 
             labeller = label_bquote(cols = bold(italic(P)^.(as.character(anc))))) +
  scale_x_continuous('', expand = c(0,0), limit = c(0,1), breaks = seq(0,1,0.1)) +
  scale_y_continuous(, expand = c(0,0), breaks = seq(0,1,0.1)) +
  facetted_pos_scales(
    y = list(
      anc %in% c('AFR','EAS','IAM') ~ scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1,0.1)),
      T ~ scale_y_continuous(expand = c(0,0), limits = c(0,0.6), breaks = seq(0,0.6,0.1))
    )
  ) +
  scale_fill_manual(name = '', values = c(AFR = '#ffc300', EUR = '#219ebc', 
                                          EAS = '#a7c957', IAM = '#e76f51'),
                    labels = c(AFR = bquote(italic(P)^'AFR'), EUR = bquote(italic(P)^'EUR'),
                               EAS = bquote(italic(P)^'EAS'), IAM = bquote(italic(P)^'IAM'))) +
  labs(y = 'Fraction of Weighted Sample') +
  theme_custom() +
  theme(legend.position = 'top', legend.box = 'vertical', legend.spacing.y = unit(0.05,'cm'),
        strip.text.y = element_blank(), strip.placement = 'outside')
ggsave(paste0('figures/appendix_gsp_histo',fig_format), width = 6, height = 6)


# 3.7 Figure 7 (appendix): supervised K = 4 stacked bar plots by region
my_admx_plot(all$s4 %>% select(-subregion, -region) %>% rename(subregion = subregion2, region = region3),
             c('Europe', 'The Americas', 'Sub-Saharan Africa', 'Eastern Asia & Oceania'), 
             by_EUR = T, c('EUR', 'EAS', 'IAM', 'AFR'))
ggsave(paste0('figures/appendix_admixture_plot_multi',fig_format), width = 16, height = 25)


# 3.8 Figure 8 (appendix): separate figure for Middle East & Southern Asia 
my_admx_plot(all$s4 %>% select(-subregion, -region) %>% rename(subregion = subregion2, region = region3),
             'Middle East & Southern Asia', by_EUR = T, c('EUR', 'EAS', 'IAM', 'AFR'))
ggsave(paste0('figures/appendix_admixture_plot_mesa',fig_format), width = 16, height = 5)


# 3.9 Figure 9 (appendix): separate figure for all best-ancestry US/Canada before recategorization
my_admx_plot(all$s4, 'USA & Canada', by_EUR = T, c('EUR', 'EAS', 'IAM', 'AFR'))
ggsave(paste0('figures/appendix_admixture_plot_usacanada',fig_format), width = 16, height = 5)


# 3.10 Figure 10 (appendix): compare ADMIXTURE and Gnomix estimates
genanc <- s4 %>%
  gather(anc, P1, paste0('X',1:4)) %>%
  mutate(anc = case_when(
    anc=='X1' ~ 'AFR', anc=='X2' ~ 'EUR', anc=='X3' ~ 'EAS', anc=='X4' ~ 'IAM',
    T ~ NA
  )) %>%
  transmute(aid, FID, IID, anc, P1) %>%
  left_join(readRDS('dp_gsp.rds') %>% # from another project
              transmute(FID, IID, anc = X, P2 = P)) %>%
  mutate(anc = factor(anc, levels = ancestry_label2))

genanc_cor <- genanc %>%
  group_by(anc) %>%
  summarise(cor = cor(P1,P2)) %>%
  mutate(x = 0.02, y = 0.98, label = paste0('rho ~ "=" ~"', sprintf('%.3f', cor),'"'))

genanc %>%
  ggplot(aes(x = P1, y = P2, color = anc)) +
  geom_point(shape = 1, size = 2, stroke = 1, alpha = 0.4) +
  geom_smooth(method = 'lm', color = 'black') +
  geom_text(data = genanc_cor, aes(x = x, y = y, label = label), color = 'black', hjust = 0, parse = T) +
  geom_abline(linetype = 'dashed', color = 'gray75') +
  facet_wrap(~ anc, labeller = label_bquote(cols = bold(italic(P)^.(as.character(anc))))) +
  scale_color_manual(values = c(AFR = '#ffc300', EUR = '#219ebc', EAS = '#a7c957',
                                IAM = '#e76f51')) +
  labs(x = 'ADMIXTURE', y = 'Gnomix') +
  theme_custom() +
  theme(legend.position = 'none')
ggsave(paste0('figures/appendix_gsp_robust_gnomix',fig_format), width = 6, height = 6)


# 3.11 Figure 11 (appendix): comparing supervised K = 4 estimates derived using different references
dat_s4 <- readRDS('supervised4_limitedref.rds') %>%
  rename_with(~ paste0('old_',.), starts_with('X')) %>%
  left_join(all$s4 %>% select(aid, FID, IID, region3)) %>%
  left_join(readRDS('supervised_output.rds')$s4 %>% 
              mutate(across(c('FID','IID'), as.character))) %>%
  mutate(color = case_when(is.na(region3)|region3=='USA & Canada' ~ 'USA, Canada, & Unknown',
                           T ~ region3),
         color = factor(color, levels = c(region_label, 'USA, Canada, & Unknown')))

p1 <- my_scatter_rug_plot(dat_s4 %>% select(x = old_X1, y = X1, color),
                          expression(atop(italic(P)^AFR, '(Limited reference)')), 
                          expression(atop(italic(P)^AFR, '(Full reference)')), F)
p2 <- my_scatter_rug_plot(dat_s4 %>% select(x = old_X2, y = X2, color),
                          expression(atop(italic(P)^EUR, '(Limited reference)')), 
                          expression(atop(italic(P)^EUR, '(Full reference)')), F)
p3 <- my_scatter_rug_plot(dat_s4 %>% select(x = old_X3, y = X3, color),
                          expression(atop(italic(P)^EAS, '(Limited reference)')), 
                          expression(atop(italic(P)^EAS, '(Full reference)')), F)
p4 <- my_scatter_rug_plot(dat_s4 %>% select(x = old_X4, y = X4, color),
                          expression(atop(italic(P)^IAM, '(Limited reference)')), 
                          expression(atop(italic(P)^IAM, '(Full reference)')), T)
ggsave(paste0('figures/appendix_gsp_robust_ref',fig_format), plot_grid(p1,p2,p3,p4, nrow = 2),
       width = 12, height = 12)


# 3.12 Figure 13 (appendix): comparing with unsupervised K = 4 estimates derived using 2 equal-sized subsamples
dat_us5 <- all$us4 %>%
  select(aid, FID, IID, region3, paste0('X',1:4)) %>%
  rename_with(~ paste0('all_',.), starts_with('X')) %>%
  right_join(readRDS('unsupervised_subsample_output_k4.rds') %>% 
               bind_rows(.id = 'run') %>% 
               mutate(across(c('FID','IID'), as.character))) %>%
  mutate(color = case_when(is.na(region3)|region3=='USA & Canada' ~ 'USA, Canada, & Unknown',
                           T ~ region3),
         color = factor(color, levels = c(region_label, 'USA, Canada, & Unknown'))) %>%
  filter(type=='estimation')

p1 <- my_scatter_rug_plot(dat_us5 %>% select(x = all_X1, y = X1, color),
                          expression(atop(italic(P)^1, '(Full sample)')),
                          expression(atop(italic(P)^1, '(50% subsample)')), F)
p2 <- my_scatter_rug_plot(dat_us5 %>% select(x = all_X2, y = X2, color),
                          expression(atop(italic(P)^2, '(Full sample)')),
                          expression(atop(italic(P)^2, '(50% subsample)')), F)
p3 <- my_scatter_rug_plot(dat_us5 %>% select(x = all_X3, y = X3, color),
                          expression(atop(italic(P)^3, '(Full sample)')),
                          expression(atop(italic(P)^3, '(50% subsample)')), F)
p4 <- my_scatter_rug_plot(dat_us5 %>% select(x = all_X4, y = X4, color),
                          expression(atop(italic(P)^4, '(Full sample)')),
                          expression(atop(italic(P)^4, '(50% subsample)')), T)

plot_grid(p1,p2,p3,p4, nrow = 2)
ggsave('figures/appendix_gsp_robust_uns.pdf', width = 12, height = 12)


# 3.13 Figure 15 (appendix): ternary plots
s4 %>%
  filter(region3 %in% c('Sub-Saharan Africa','Europe','The Americas'),
         X3 > 0.05) %>%
  pull(region3) %>% 
  table()
# 21 Sub-Saharan Africa
# 38 Europe
# 51 The Americas

s4 %>%
  filter(region3 %in% c('Sub-Saharan Africa','Europe','The Americas'),
         X3 <= 0.05) %>%
  ggtern(aes(x = X1, y = X2, z = X4, color = region3)) + # X2+X3
  geom_point(alpha = 0.2, size = 3, show.legend = F) +
  geom_point(data = s4 %>%
               filter(region3 %in% c('Sub-Saharan Africa','Europe','The Americas'),
                      X3 <= 0.05) %>%
               group_by(region3) %>%
               summarise(across(paste0('X',1:4), mean)),
             shape = 1, color = 'black', size = 3, stroke = 1.5) +
  facet_wrap(~ region3,
             labeller = as_labeller(c(
               'Sub-Saharan Africa' = 'Self-reported\nancestral Sub-Saharan Africans',
               'Europe' = 'Self-reported\nancestral Europeans',
               'The Americas' = 'Self-reported\nancestral Americans'
             ))) +
  scale_color_manual(name = '', values = c(`Europe`='#007ea7', `The Americas`='#dd0426',
                                           `Sub-Saharan Africa`='#f9a620')) +
  labs(x = expression(italic(P)^AFR), y = expression(italic(P)^EUR), # italic(P)^EAS
       z = expression(italic(P)^IAM)) +
  theme_custom()
ggsave(paste0('figures/appendix_ternary_easfilter',fig_format), width = 15, height = 6)

