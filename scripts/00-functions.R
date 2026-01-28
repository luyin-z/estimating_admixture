# 0. Init ----------------------------------------------------------------------
library(tidyverse); library(RColorBrewer); library(modelsummary); library(gtsummary)
library(cowplot); library(gridExtra); library(Hmisc)

race_label <- c('NHW', 'NHB', 'AINA', 'AAPI', 'Hispanic')
region_label <- c('Europe', 'The Americas', 'Sub-Saharan Africa', 'Eastern Asia & Oceania', 'Middle East & Southern Asia')
subregion_label <- c('Sub-Saharan Africa', 'British Isles', 'East Europe', 'North Europe',
                  'South Europe', 'West Europe', 'East Asia', 'Oceania', 'Philippines',
                  'South Asia', 'Southeast Asia', 'Middle East', 'America', 'Caribbean',
                  'Central America', 'Canada', 'Mexico', 'Native America', 'South America',
                  'Other')

theme_custom <- function() {
  theme_classic() %+replace%
    theme(plot.margin = unit(c(5.5,10,5.5,5.5), 'pt'),
          plot.title.position = 'plot', 
          plot.caption.position = 'plot',
          strip.background = element_blank(), 
          strip.text = element_text(face='bold', size=rel(0.8), margin=margin(0,0,2,0)),
          plot.title = element_text(face='bold', size=rel(1.2), hjust=0,
                                    margin=margin(0,0,5.5,0)))
}



# 1. Admixture plot ------------------------------------------------------------
my_admx_plot <- function(data, regions, by_EUR, stacking_order) {
  # genetic similarity proportions in descending orders
  order <- data %>%
    filter(!is.na(region)) %>%
    group_by(region) %>%
    summarise(across(starts_with('X'), mean, na.rm = T)) %>%
    gather(anc, value, starts_with('X')) %>%
    arrange(region, -value)
  
  # prepare data for visualization
  plot_data <- list()
  for (reg in regions) {
    if (by_EUR) {
      plot_data[[reg]] <- data %>% 
        filter(region==reg) %>%
        arrange(region, -X2) %>% 
        mutate(n = row_number()) %>%
        gather('X', 'pc', starts_with('X'))
    } else {
      order_x <- order[order$region==reg,] %>% pull(anc)
      plot_data[[reg]] <- data %>% 
        filter(region==reg) %>%
        arrange(across(order_x, desc)) %>% 
        mutate(n = row_number()) %>%
        gather('X', 'pc', starts_with('X'))
    }
  }
  plot_data <- bind_rows(plot_data) %>%
    mutate(region = factor(region, levels = regions),
           X = case_when(X=='X1' ~ 'AFR',
                         X=='X2' ~ 'EUR',
                         X=='X3' ~ 'EAS',
                         X=='X4' ~ 'IAM', 
                         X=='X5' ~ 'MEA',
                         T ~ NA) %>% factor(levels = stacking_order))
  
  # add sample size
  text_data <- plot_data %>% 
    select(aid, region) %>%
    unique() %>%
    group_by(region) %>%
    summarise(n = max(row_number())) %>%
    mutate(text = paste0(region,'\n(N = ',n,')')) %>%
    arrange(-n)
  
  # plot
  admixture_plot <- plot_data %>%
    mutate(region = factor(region, levels = text_data$region, labels = text_data$text)) %>%
    ggplot() +
    geom_bar(aes(x = n, y = pc, fill = X), stat = 'identity', position = 'stack', width = 1) +
    facet_wrap(vars(region), nrow = 6, scales = 'free_x') +
    scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = '', y = 'Ancestry', fill = '') +
    scale_fill_manual(values = c(AFR = '#ffc300', EUR = '#219ebc', EAS = '#a7c957',
                                 IAM = '#e76f51', MEA = '#6a4c93'),
                      labels = c(AFR = bquote(italic(P)^'AFR'), EUR = bquote(italic(P)^'EUR'),
                                 EAS = bquote(italic(P)^'EAS'), IAM = bquote(italic(P)^'IAM'),
                                 MEA = bquote(italic(P)^'MEA'))) +
    geom_hline(yintercept = seq(0, 1, 0.125), color = 'gray', linetype = 'dotted', linewidth = 0.3) +
    theme_custom() +
    theme(legend.position = 'top', strip.text = element_text(size=rel(1.5)), legend.text = element_text(size=rel(1.5)))
  
  return(admixture_plot)
}



# 2. Scatter plots with marginal density ---------------------------------------------------
my_scatter_rug_plot <- function(data, x_label, y_label, legend_label = F) {
  data <- data %>% filter(!is.na(color), x>=0.01, y>=0.01)
  
  # create texts for correlation and beta coefficients
  n <- nrow(data)
  rho <- cor(data$x, data$y)
  beta <- coef(lm(y~x, data))[[2]]
  rho_exp <- paste0('rho ~ "=" ~"', sprintf('%.3f', rho),'"')
  beta_exp <- paste0('beta ~ "=" ~"', sprintf('%.3f', beta),'"')
  n_exp <- paste0('N ~ "=" ~"', n,'"')
  
  # scatter plot
  scatter <- data %>%
    ggplot(aes(x = x, y = y, color = color)) + 
    geom_point(shape = 1, size = 2, stroke = 1, alpha = 0.4) +
    scale_color_manual(values = c(`Europe`='#007ea7', `The Americas`='#dd0426', `Sub-Saharan Africa`='#f9a620',
                                  `Eastern Asia & Oceania`='#006400', `Middle East & Southern Asia`='#7678ed', 
                                  `USA, Canada, & Unknown`='black')) +
    geom_abline(slope = 1, color = 'black', linetype = 'dashed', linewidth = 1) +
    annotate('text', x = 0.1, y = 0.9, parse = T, label = rho_exp) +
    annotate('text', x = 0.1, y = 0.85, parse = T, label = beta_exp) +
    annotate('text', x = 0.1, y = 0.8, parse = T, label = n_exp) +
    scale_x_continuous(name = x_label, limits = c(0,1), expand = c(0,0)) +
    scale_y_continuous(name = y_label, limits = c(0,1), expand = c(0,0)) +
    theme_custom() +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
  
  # legend for the scatter plot
  if (legend_label) {
    scatter <- scatter + 
      theme(legend.position = c(0.8,0.3)) + 
      labs(color = '') +
      guides(color = guide_legend(override.aes = list(shape = 15,size = 5,alpha = 1)))
  } else {
    scatter <- scatter + 
      theme(legend.position = 'none')
  }
  
  # rug plot
  rug <- data %>%
    ggplot(aes(x = x, y = y, color = color)) +
    geom_segment(aes(color = color, xend = x,
                     y = 1+as.numeric(color)*0.03+0.01, 
                     yend = 1+as.numeric(color)*0.03-0.01), alpha = 0.1) + 
    scale_color_manual(values = c(`Europe`='#007ea7', `The Americas`='#dd0426', `Sub-Saharan Africa`='#f9a620',
                                  `Eastern Asia & Oceania`='#006400', `Middle East & Southern Asia`='#7678ed', `0`='black')) +
    scale_x_continuous(name = x_label, limits = c(0,1), expand = c(0,0)) +
    theme_custom() +
    theme(axis.title = element_blank(), axis.text = element_blank(), 
          axis.ticks = element_blank(), axis.line = element_line(color = 'white'),
          legend.position = 'none')
  
  # combine the scatter and rug plots
  if (length(unique(data$color))==1) {
    p <- scatter
  } else {
    p <- plot_grid(rug, scatter, nrow = 2, align = 'v', axis = 'l', rel_heights = c(1,5))
  }
  
  return(p)
}



# 3. Summary tables ------------------------------------------------------------
# 3.1 by single best geographic ancestry
# function for calculate the mean and sd for genetic similarity proportions by geographic ancestry
add_mean_sd_general <- function(data, variable, keep_var, ...) {
  data %>%
    group_by(.data[[variable]]) %>%
    arrange(.data[[variable]]) %>%
    summarise_all(list(`_mean` = mean, `_sd` = sd)) %>%
    gather(key = 'key', value = 'value', -by_var) %>%
    separate(key, c('key', 'stat'), '__') %>%
    spread(stat, value) %>%
    mutate(stat = sprintf('%.2f',mean)) %>%
    select(-c(mean, sd)) %>%
    spread(key, stat) %>%
    select(all_of(keep_var))
}

# function for creating summary statistic table based on geographic ancestry
summarytable <- function(data, results, drop_smallgroup, no_AmericaCanada,
                         single_anc, ancestry_subsample, by_var_label) {
  
  add_mean_sd <- partial(add_mean_sd_general, keep_var = results)
  
  if (drop_smallgroup) {
    data <- data %>% filter(!small_group)
  }
  if (no_AmericaCanada) {
    data <- data %>% filter(by_var!='America', by_var!='Canada')
  }
  if (single_anc) {
    data <- data %>% filter(single_anc==1)
  }
  if (ancestry_subsample) {
    data <- data %>% filter(ancestry_subsample==1)
  }
  # data for the summary statistic table
  tab_dat <- data %>% 
    # remove other and NA ancestries from the table
    filter(!(by_var %in% c('Other', NA))) %>%
    select(all_of(results), by_var) %>%
    group_by(by_var) %>%
    mutate(
      size = max(row_number()),
      size = ifelse(is.na(by_var), NA, size),
      by_var = ifelse(is.na(by_var), 'NA', as.character(by_var))
    ) 
  
  # rank geographic ancestry categories by frequency
  region_label_freq <- tab_dat[,c('size','by_var')] %>% 
    unique() %>%
    arrange(-size) %>% 
    pull(by_var) %>%
    as.character()
  tab_dat <- tab_dat %>%
    mutate(by_var = factor(by_var, levels = region_label_freq))
  
  # create table
  t <- tab_dat %>%
    # proportion of each geographic ancestry category
    tbl_summary(     
      include = by_var,
      statistic = all_categorical() ~ '{n} ({p}%)', # '{n}/n({p}%)
      digits = everything() ~ 0,
      label = by_var ~ by_var_label
    ) %>% 
    # add mean and sd for each genetic similarity proportion by geographic ancestry
    add_stat(
      all_categorical() ~ add_mean_sd,
      location = all_categorical() ~ 'level'
    ) %>%
    # remove footnote
    modify_footnote(c(all_stat_cols()) ~ NA)
  
  return(t)
}


# 3.2 by fractional geographic ancestry
# calculate average genetic similarity proportions by fractional geographic ancestry
create_mean <- function(data, results, by_var) {
  if (by_var=='Subregion') {
    prefix <- 'SUBREG '
  } else if (by_var=='Region') {
    prefix <- 'REGION '
  } else {
    print('no such region variables')
  }
  
  # average genetic similarity proportions
  s <- list()
  for (c in names(data %>% select(starts_with(prefix)))) {
    for (a in names(data %>% select(all_of(results)))) {
      dat <- data[data[[c]]>0,] %>% 
        select(c, a) %>%
        drop_na()
      wm <- weighted.mean(dat[[a]], dat[[c]])
      s[[c]][[a]] <- sprintf('%.2f', wm)
    }
    s[[c]] <- bind_rows(s[[c]])
  }
  s <- bind_rows(s, .id = by_var)
  
  # sample size
  td <- fam_origin %>%
    filter(aid %in% data$aid) %>%
    select(starts_with(prefix)) %>%
    colSums(na.rm = T) %>%
    as.data.frame() %>%
    rownames_to_column(var = by_var)
  names(td)[2] <- 'N'
  
  # merge together
  sum_tbl <- td %>% left_join(s)
  
  return(sum_tbl)
}

# function for creating summary statistic tables based on fractional geographic ancestry
summarytable2 <- function(data, results) {
  # subregion-level table
  sum_tbl <- create_mean(data, results, 'Subregion') %>%
    mutate(Subregion = str_remove(Subregion, 'SUBREG ')) %>%
    left_join(meta %>% select(-country) %>% unique() %>% rename_with(str_to_title)) %>%
    filter(!(Subregion %in% c('America', 'Canada', 'Other'))) %>%
    mutate(prop = N/sum(N))
  
  # region-level table
  data2 <- data %>% filter(`SUBREG Canada`==0 & `SUBREG America`==0 & `SUBREG Other`==0)
  sum_tbl2 <- create_mean(data2, results, 'Region') %>%
    mutate(Region = str_remove(Region, 'REGION ')) %>%
    mutate(Subregion = 'Total', prop = N/sum(N))
  
  # combine subregion & region tables
  total_tbl <- sum_tbl %>%
    rbind(sum_tbl2) %>%
    drop_na() %>%
    filter(Subregion!='Sub-Saharan Africa') %>%
    mutate(Subregion = factor(Subregion, levels = c('Total', subregion_label)),
           Region = factor(Region, levels = region_label)) %>%
    select(Region, Subregion, N, prop, all_of(results)) %>%
    arrange(Region, Subregion)
  names(total_tbl) <- gsub('N', paste0('N = ', sprintf('%.0f',sum(sum_tbl$N))), names(total_tbl))

  return(total_tbl)
}

# find the 10th and 90th percentiles of genetic similarity proportions by fractional geographic ancestry
create_percentile <- eval(parse(text = deparse(create_mean) %>%
                                  str_replace_all('weighted\\.mean\\(dat\\[\\[a\\]\\], dat\\[\\[c\\]\\]\\)', 
                                                  'wtd.quantile(dat[[a]], weights = dat[[c]], probs = c(0.1, 0.9))') %>%
                                  str_replace_all('sprintf\\(\\"%.2f\\", wm\\)', 
                                                  "paste0(\"(\", sprintf(\"%.2f\", wm[[1]]), \", \", sprintf(\"%.2f\", wm[[2]]), \")\")")))

summarytable3 <- eval(parse(text = deparse(summarytable2) %>%
                              str_replace_all('create_mean', 'create_percentile')))

