library(here); library(tidyverse); library(bigsnpr); library(data.table)
library(rhdf5)

setwd(str_remove(here(), '/scripts'))

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

# read in 5 PCs
ah_PCs <- readRDS('ah_projected_PC5.rds')

# read in supervised ADMIXTURE 5
ah_s5 <- readRDS('all_data.rds')$s5

# merge into one dataset
ah <- ah_PCs %>%
  mutate(across(c(FID,IID), function(x) x = as.character(x))) %>% 
  left_join(ah_s5)

# standardize PCs based on weighted means and SDs
mean <- sd <- list()
for (i in 1:5) {
  l <- paste0('PC',i)
  td <- ah %>% filter(!is.na(ipwgt)) # 7884 left -- 2090 missing, majority due to missing racec (N = 1614)
  mean[[l]] <- weighted.mean(td[[l]], td$ipwgt)
  sd[[l]] <- sqrt(sum((td[[l]] - mean[[l]])^2*td$ipwgt)/(sum(td$ipwgt)-1))
  ah[[l]] <- (ah[[l]] - mean[[l]])/sd[[l]]
}

# d <- ah %>%
#   filter(include_ssa==1) %>%
#   mutate(region = 'Sub-Saharan Africa') %>%
#   rbind(ah %>%
#           filter(include_america==1) %>%
#           mutate(region = 'The Americas')) %>%
#   rbind(ah %>%
#           filter(include_europe==1) %>%
#           mutate(region = 'Europe')) %>%
#   rename(color = region)
d <- ah %>% rename(color = region3)

# define a function for comparing one GSP and one PC
my_compare_PCX <- function(data, x_label, y_label, color_label, base_color) {
  data <- data %>% filter(!is.na(color), y>=0.01)
  n_rho <- data %>%
    group_by(color, group) %>%
    summarise(n = n(),
              rho = cor(.data$x, .data$y)) %>%
    mutate(text = paste0("atop(rho ~ \"=\" ~ ", sprintf("%.3f", rho), "*\",\"~ N ~ \"=\" ~", n, ")"),
           y = c(0.95, 0.9))
  
  scatter <- data %>%
    ggplot(aes(x = x, y = y, color = group)) + 
    geom_point(shape = 1, size = 2, stroke = 1) +
    geom_text(data = n_rho, aes(x = -Inf, y = y, label = text), parse = T, hjust = -0.2, show.legend = F) +
    scale_x_continuous(name = x_label, expand = c(0,0)) +
    scale_y_continuous(name = y_label, limits = c(0,1), expand = c(0,0)) +
    scale_color_manual(color_label, values = c(base_color, 'black')) +
    theme_custom() +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
          legend.position=c(0.9,0.4))
  
  rug <- data %>%
    ggplot(aes(x = x, y = y)) +
    geom_segment(aes(xend = x, y = 1, yend = 1.2), color = base_color, alpha = 0.1) + 
    scale_x_continuous(name = x_label, expand = c(0,0)) +
    theme_custom() +
    theme(axis.title = element_blank(), axis.text = element_blank(), 
          axis.ticks = element_blank(), axis.line = element_line(color = 'white'),
          legend.position = 'none')
  
  p <- plot_grid(rug, scatter, nrow = 2, align = 'v', axis = 'l', rel_heights = c(1,8))

  return(p)
}

p1 <- my_compare_PCX(d %>% rename(x = PC1, y = X1) %>% 
                       filter(color=='Sub-Saharan Africa') %>%
                       mutate(group = X4>=0.05),
                     'PC1', expression(italic(P)^AFR), expression(italic(P)^IAM>=0.05), '#f9a620')
p2 <- my_compare_PCX(d %>% rename(x = PC2, y = X4) %>% 
                       filter(color=='The Americas') %>%
                       mutate(group = X1>=0.05),
                     'PC2', expression(italic(P)^IAM), expression(italic(P)^AFR>=0.05), '#dd0426')
p3 <- my_compare_PCX(d %>% rename(x = PC1, y = X1) %>%
                       filter(color=='The Americas') %>%
                       mutate(group = X4>=0.05),
                     'PC1', expression(italic(P)^AFR), expression(italic(P)^IAM>=0.05), '#dd0426')
ggsave('figures/appendix_figure6.pdf', plot_grid(p1,p2,p3, nrow = 2), width = 10, height = 10) 

