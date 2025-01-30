# 0. Init ----------------------------------------------------------------------
library(here); library(data.table); library(tidyverse); library(rnaturalearth); library(giscoR); library(sf)
library(cowplot); library(ggtext)

setwd(str_remove(here(),'/scripts'))



# 1. Import & prepare data -----------------------------------------------------
all <- readRDS('all_data.rds') %>%
  map(~{
    .x <- .x %>%
      mutate(`SUBREG US & Canada` = `SUBREG America` + `SUBREG Canada`) %>%
      select(-c(`SUBREG America`, `SUBREG Canada`))
  })
hapmap3 <- fread('hapmap3/relation_4a.txt', select = c(2,7), col.names = c('IID', 'pop')) %>%
  filter(pop != 'GIH')
hgdp <- fread('hgdp/hgdp_all.psam', select = c(1,3), col.names = c('IID', 'pop'))
ref <- fread('ah_hapmap3_hgdp_s5_nosex.fam', select = 2, col.names = 'IID')


# 1.1 create subregion-level ancestry proportions and sample size
freq_anc <- all$s5 %>%
  select(starts_with('SUBREG ')) %>%
  drop_na() %>%
  colSums(na.rm = T) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'subregion')
names(freq_anc) <- c('subregion', 'n')

s <- list()
for (c in names(all$s5 %>% select(starts_with('SUBREG ')))) {
  for (a in paste0('X',1:5)) {
    dat <- all$s5[all$s5[[c]]>0,] %>% 
      select(c, a) %>%
      drop_na()
    s[[c]][[a]] <- weighted.mean(dat[[a]], dat[[c]])
  }
  s[[c]] <- bind_rows(s[[c]])
}
freq_anc <- freq_anc %>% 
  left_join(bind_rows(s, .id = 'subregion')) %>%
  mutate(subregion = str_remove(subregion, 'SUBREG '), n = round(n)) %>%
  filter(!subregion %in% c('Other','US & Canada'))


# 1.2 create country-level sample size
cntr_freq <- all$s5 %>%
  select(paste('SUBREG', c('Sub-Saharan Africa', 'British Isles', 'South America', 'Central America')),
         starts_with('CNTR ')) %>%
  mutate(`CNTR United States` = `CNTR America`, `CNTR Native America` = `CNTR American Indian`,
         `CNTR China` = `CNTR China` + `CNTR Taiwan (China)`) %>%
  select(-c(`CNTR America`, `CNTR American Indian`, `CNTR Taiwan (China)`)) %>%
  drop_na() %>%
  colSums(na.rm = T) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'country') %>%
  mutate(country = country %>% str_remove('CNTR ') %>% str_remove('SUBREG '))
names(cntr_freq) <- c('country', 'n')


# 1.3 prepare subregion-level geometry data
# assign geounit into 'subregion'-level categories
subregion_grp <- ne_countries(scale = 'medium', returnclass = 'sf',
                           # use map_units so French Guiana is marked as South America
                           type = 'map_units') %>%
  mutate(
    subregion = case_when(
      # Middle East
      geounit=='Iran' ~ 'Middle East',
      # British Isles
      geounit %in% c('England', 'Guernsey', 'Isle of Man', 'Jersey', 'Northern Ireland',
                     'Scotland', 'Wales', 'Ireland') ~ 'British Isles',
      # Oceania
      subregion %in% c('Australia and New Zealand', 'Melanesia', 'Micronesia', 'Polynesia') ~ 'Oceania',
      # Sub-Saharan Africa
      region_wb=='Sub-Saharan Africa' | geounit=='Djibouti' ~ 'Sub-Saharan Africa',
      # Mexico & Philippines
      geounit %in% c('Mexico', 'Philippines') ~ geounit,
      subregion=='Seven seas (open ocean)' ~ region_wb,
      # US & Canada
      geounit %in% c('United States of America', 'Canada') ~ 'US & Canada',
      T ~ subregion),
    subregion = case_when(
      subregion %in% c('Northern Africa', 'Western Asia') ~ 'Middle East',
      subregion=='South-Eastern Asia' ~ 'Southeast Asia',
      subregion %in% c('Central Asia', 'Southern Asia') ~ 'South Asia',
      T ~ subregion),
    subregion = str_remove_all(subregion, 'ern')
  )

# general geometry data
world_data <- gisco_get_countries() %>%
  # remove Antarctica
  filter(ISO3_CODE!='ATA') %>%
  left_join(subregion_grp %>%
              as.data.frame() %>%
              select(geounit, ISO3_CODE = iso_a3, subregion)) %>%
  mutate(
    subregion = case_when(
      ISO3_CODE %in% c('ATG') ~ 'Caribbean',
      ISO3_CODE %in% c('GEO') ~ 'East Europe',
      ISO3_CODE %in% c('BEL', 'FRA') ~ 'West Europe',
      ISO3_CODE %in% c('BIH','SRB', 'PRT') ~ 'South Europe',
      ISO3_CODE %in% c('PSE') ~ 'Middle East',
      ISO3_CODE %in% c('PNG') ~ 'Oceania',
      ISO3_CODE %in% c('GBR') ~ 'British Isles',
      ISO3_CODE %in% c('XC', 'XE') ~ 'East Asia',
      ISO3_CODE %in% c('XU', 'XG') ~ 'Sub-Saharan Africa',
      ISO3_CODE %in% c('XH', 'XD') ~ 'South Asia',
      T ~ subregion
    )) %>%
  filter(subregion!='Antarctica')

# geometry data for French Guiana
french_guiana <- subregion_grp %>%
  filter(name=='French Guiana') %>%
  mutate(geometry = st_cast(geometry, 'GEOMETRY'))

# geometry data for New Mexico (NM) & Arizona (AZ)
native_america <- ne_states(country = 'united states of america', returnclass = 'sf') %>%
  filter(name %in% c('New Mexico', 'Arizona', 'Utah', 'Colorado')) %>% # Nevada
  summarise() %>%
  mutate(subregion = 'Native America')

# geometry data for US & Canada combined
USC_data <- world_data %>%
  filter(subregion=='US & Canada') %>%
  summarise()

# create 'subregion'-level geometry data
subregion_data <- world_data %>%
  group_by(subregion) %>%
  summarise()
subregion_data[subregion_data$subregion=='West Europe',]$geometry <- 
  # remove French Guiana from the 'West Europe' area
  st_difference(subregion_data[subregion_data$subregion=='West Europe',]$geometry, french_guiana$geometry)

# get the centroid coordinates from the 'Native America' geometry data
NA_coord <- st_coordinates(st_centroid(native_america)) %>%
  cbind(native_america)


# 1.4 prepare country-level geometry data
country_grp <- world_data %>%
  # assign geounit into 'country'-level categories
  mutate(country = case_when(
    NAME_ENGL %in% c('Bahrain', 'Iran', 'Iraq', 'Israel', 'Jordan', 'Kuwait', 'Palestine', 
                     'Oman', 'Qatar', 'Saudi Arabia', 'Syria', 'Turkey', 'United Arab Emirates',
                     'Yemen') ~ 'Other Middle East',
    subregion %in% c('Sub-Saharan Africa', 'British Isles', 'Philippines', 'Mexico', 'Canada', 
                  'South America', 'Central America') ~ subregion,
    NAME_ENGL %in% c('Czechia', 'Slovakia') ~ 'Czech & Slovakia',
    NAME_ENGL %in% c('Albania', 'Bulgaria', 'Estonia', 'Latvia', 'Lithuania', 'Romania') ~ 'Other East Europe',
    NAME_ENGL=='Russian Federation' ~ 'Russia',
    NAME_ENGL=='Iceland' ~ 'Iceland and Atlantic Islands',
    NAME_ENGL %in% c('Belgium', 'Luxembourg', 'Netherlands') ~ 'Belgium, Luxembourg, and Netherlands',
    NAME_ENGL %in% c('Bangladesh', 'Nepal', 'Pakistan', 'Sri Lanka') ~ 'Bangladesh, Nepal, Pakistan, and Sri Lanka',
    NAME_ENGL %in% c('North Korea', 'South Korea') ~ 'North and South Korea',
    NAME_ENGL=='Hong Kong' ~ 'Hong Kong (China)',
    NAME_ENGL %in% c('Fiji', 'Guam', 'New Zealand', 'Tonga', 'American Samoa', 'Australia', 
                     'Azores', 'Bermuda', 'Falkland Islands') ~ 'Pacific Islands, Australia, and Atlantic Islands',
    NAME_ENGL %in% c('Brunei', 'Cambodia', 'Indonesia', 'Laos', 'Malaysia', 'Myanmar/Burma', 
                     'Singapore', 'Thailand') ~ 'Other Southeast Asia',
    NAME_ENGL %in% c('Antigua and Barbuda', 'Aruba', 'Bahamas', 'Barbados', 'Cayman Islands', 'Dominica', 
                     'Dominican Republic', 'Grenada', 'Montserrat', 'Saint Kitts and Nevis', 'Saint Lucia',
                     'Saint Vincent and The Grenadines', 'Trinidad and Tobago', 'Us Virgin Islands') ~ 'Other West Indies',
    T ~ NAME_ENGL
  ))

# create 'country'-level geometry data
country_data <- country_grp %>%
  group_by(country) %>%
  summarise()
country_data[country_data$country=='France',]$geometry <- 
  # remove French Guiana from the France area
  st_difference(country_data[country_data$country=='France',]$geometry, french_guiana$geometry)

# country without sample coverage
country_NA <- country_data %>%
  left_join(cntr_freq) %>%
  filter(is.na(n)) %>%
  summarise()


# 1.5 prepare coordinates for HapMap3 & HGDP samples
ref <- ref %>% 
  inner_join(hgdp %>% rbind(hapmap3)) %>%
  group_by(pop) %>%
  summarise(n = n()) %>%
  mutate(ancestry = case_when(
    pop=='YRI' ~ 'Sub-Saharan African',
    pop=='CEU' ~ 'European',
    pop=='CHB' ~ 'East Asian',
    pop %in% c('Maya', 'Pima') ~ 'Indigenous American',
    T ~ 'Middle Eastern'
  ),
  country = case_when(
    pop=='YRI' ~ 'Nigeria',
    pop=='CEU' ~ 'United States',
    pop=='CHB' ~ 'China',
    pop %in% c('Maya', 'Pima') ~ 'Mexico',
    T ~ 'Israel'
  ),
  state = case_when(
    pop=='YRI' ~ 'Oyo', # Ibadan
    pop=='CEU' ~ 'Utah',
    pop=='CHB' ~ 'Beijing',
    pop=='Bedouin' ~ 'HaDarom', # Negev (Southern Israel)
    pop=='Druze' ~ 'HaZafon', # Carmel (Northern Israel)
    pop=='Palestinian' ~ 'HaMerkaz', # Central Israel
    pop=='Maya' ~ 'Campeche',
    pop=='Pima' ~ 'Sonora'
  ))

oyo <- ne_states(country = 'nigeria', returnclass = 'sf') %>%
  filter(name=='Oyo')
utah <- ne_states(country = 'united states of america', returnclass = 'sf') %>%
  filter(name=='Utah')
beijing <- ne_states(country = 'china', returnclass = 'sf') %>%
  filter(name=='Beijing')
israel <- ne_states(country = 'israel', returnclass = 'sf') %>%
  filter(name %in% c('HaDarom', 'HaZafon', 'HaMerkaz'))
mexico <- ne_states(country = 'mexico', returnclass = 'sf') %>%
  filter(name %in% c('Campeche', 'Sonora'))

ref_state <- rbind(oyo, utah, beijing, israel, mexico) %>%
  select(state = name, geometry) %>%
  right_join(ref)
ref_state <- st_coordinates(st_centroid(ref_state)) %>%
  cbind(ref_state)



# 2. World map -----------------------------------------------------------------
# 2.1 merge geometry data with attributes
world_data <- world_data %>% 
  left_join(freq_anc) %>%
  select(names(native_america)) %>%
  rbind(native_america)
french_guiana <- french_guiana %>% left_join(freq_anc)
native_america <- native_america %>% left_join(freq_anc)
NA_coord <- NA_coord %>% left_join(freq_anc)
subregion_data <- subregion_data %>% left_join(freq_anc)


# 2.2 visualize average GSPs
my_gsp_map <- function(gsp, anc_label, color) {
  gsp_map <- ggplot() +
    geom_sf(data = subregion_data, aes(fill = !!sym(gsp))) +
    geom_sf(data = french_guiana, aes(fill = !!sym(gsp))) +
    geom_point(data = NA_coord, aes(x = X, y = Y, fill = !!sym(gsp)), size = 4, shape = 21) +
    geom_sf(data = country_NA, fill = 'gray75') +
    scale_fill_gradient(name = anc_label,
                        low = 'white', high = color, na.value = 'gray75', limits = c(0,1)) +
    theme_void() +
    theme(legend.position = 'bottom')
  gsp_map <- plot_grid(gsp_map + theme(legend.position = 'none'), 
                       ggdraw(get_plot_component(gsp_map,'guide-box-bottom',T)), ncol = 1, rel_heights = c(5,1))
  return(gsp_map)
}

maps <- list()
maps$afr_map <- my_gsp_map('X1', expression(italic('P')^AFR~'(Sub-Saharan African)'), '#ffc300')
maps$eur_map <- my_gsp_map('X2', expression(italic('P')^EUR~'(European)'), '#219ebc')
maps$eas_map <- my_gsp_map('X3', expression(italic('P')^EAS~'(East Asian)'), '#a7c957')
maps$iam_map <- my_gsp_map('X4', expression(italic('P')^IAM~'(Indigenous American)'), '#e76f51')
maps$mea_map <- my_gsp_map('X5', expression(italic('P')^MEA~'(Middle Eastern)'), '#6a4c93')

  
# 2.3 location of HapMap3 & HGDP references
maps$ref_map <- ggplot() +
  geom_sf(data = subregion_data, fill = 'white') +
  geom_point(data = ref_state, aes(x = X, y = Y, color = ancestry, size = n)) +
  scale_color_manual(name = 'Reference\npopulation', drop = F, 
                     values = c(`Sub-Saharan African`='#ffc300', `European`='#219ebc',
                                `East Asian`='#a7c957', `Indigenous American`='#e76f51',
                                `Middle Eastern`='#6a4c93'),
                     labels = c(`Sub-Saharan African`='Yoruba in Nigeria',
                                `European`='Northern & Western European in Utah',
                                `East Asian`='Han Chinese in Beijing', 
                                `Indigenous American`='Pima & Mayan in Mexico',
                                `Middle Eastern`='Palestinian, Bedouin, & Druze, in Israel')) +
  scale_size_continuous(name = 'Sample size', breaks = c(25, 50, 100)) +
  theme_void() +
  theme(legend.position = 'bottom', legend.box = 'vertical') +
  guides(color = guide_legend(nrow=2,order=2), size = guide_legend(order=1))
maps$ref_map <- plot_grid(maps$ref_map + theme(legend.position = 'none'), 
                          ggdraw(get_plot_component(maps$ref_map,'guide-box-bottom',T)),
                          ncol = 1, rel_heights = c(5,1))


# 2.4 sample coverage
subregion_levels <- c('fake - Sub-Saharan Africa', 'Sub-Saharan Africa', paste('PH1',1:4), 
                   'fake - Europe', 'British Isles', 'West Europe', 'South Europe', 'East Europe', 'North Europe',
                   'fake - Middle East & South Asia', 'Middle East', 'South Asia', paste('PH2',1:3),
                   'fake - Eastern Asia & Oceania', 'East Asia', 'Southeast Asia', 'Philippines', 'Oceania', 'PH3',
                   'fake - The Americas', 'Mexico', 'Caribbean', 'Central America', 'South America', 'Native America')

N <- freq_anc[match(subregion_levels, freq_anc$subregion), ]$n
N <- ifelse(is.na(N), '', paste0('(N = ',N,')'))

lbl <- c('Sub-Saharan Africa', 'Sub-Saharan Africa', rep('',4),
         'Europe', 'British Isles', 'West Europe', 'South Europe', 'East Europe', 'North Europe',
         'Middle East & South Asia', 'Middle East', 'South Asia', rep('',3),
         'Eastern Asia & Oceania', 'East Asia', 'Southeast Asia', 'Philippines', 'Oceania', '',
         'The Americas', 'Mexico', 'Caribbean', 'Central America', 'South America', 'Native America')
lbl <- paste(lbl, N)

lbl_all <- c(expression(bold('Sub-Saharan Africa')), lbl[2:6],
             expression(bold('Europe')), lbl[8:12],
             expression(bold('Middle East & South Asia')), lbl[14:18],
             expression(bold('Eastern Asia & Oceania')), lbl[20:24],
             expression(bold('The Americas')), lbl[26:30])

col <- c(NA, '#ec9f05', rep(NA,4), NA, '#911eb4', '#ea638c', '#fabed4', '#b8b8ff',
         '#dcbeff', NA, '#023e8a', '#0096c7', rep(NA,3), NA, '#00b4d8', '#90e0ef', 
         '#03045e', '#faf0ca', NA, NA, '#43aa8b', '#80ed99', '#a7c957', '#ecf39e', '#83c5be')

maps$subregion_map <- ggplot() +
  geom_sf(data = world_data, aes(fill = factor(subregion, levels = subregion_levels),
                                 color = factor(subregion, levels = subregion_levels))) +
  geom_sf(data = world_data, fill = NA) +
  geom_sf(data = french_guiana, aes(fill = factor(subregion, levels = subregion_levels),
                                    color = factor(subregion, levels = subregion_levels))) +
  geom_sf(data = french_guiana, fill = NA) +
  geom_sf(data = world_data %>% filter(subregion=='North America'), fill = 'gray75') +
  geom_sf(data = country_NA, fill = 'gray75') +
  geom_sf(data = USC_data, fill = 'gray75') +
  geom_point(data = NA_coord, aes(x = X, y = Y), fill = '#83c5be', size = 4, shape = 21) +
  scale_fill_manual(name = '', values = col, label = lbl_all, na.value = 'white', drop = F) +
  scale_color_manual(name = '', values = col, label = lbl_all, na.value = 'white', drop = F) +
  theme_void() +
  theme(legend.key = element_rect(color = NA),
        legend.margin = margin(0, 10, 0, 0),
        legend.key.size = unit(0.15, 'inches'),
        legend.position = 'bottom',
        legend.spacing.x = unit(0.15, 'inches')) +
  guides(fill = guide_legend(nrow = 6))



# 3. Export --------------------------------------------------------------------
ggsave('figures/figure1.pdf', width = 16, height = 15,
       plot_grid(maps$afr_map, maps$eur_map, maps$eas_map, maps$iam_map, maps$mea_map, maps$ref_map, ncol = 2, nrow = 3))
ggsave('figures/appendix_figure1.pdf', maps$subregion_map, width = 16, height = 12)

