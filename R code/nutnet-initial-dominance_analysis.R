## Load Libraries

## data wrangling
library(tidyverse)
library(broom)

## analysis
library(codyn)
library(nlme)
library(lme4)
library(emmeans)
library(partR2)
library(car)
library(MuMIn)
#library(RRPP)
library(merTools)
library(trajr)

## graphics
library(ggplot2)
library(cowplot)
library(GGally)
library(RColorBrewer)
library(gridExtra)
library(grid)

theme_figs <- theme_classic() +
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 24),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = 'right', legend.justification = c(1, 0.4),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size=18,face='bold'),legend.text = element_text(size = 18),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        strip.text =element_text(size=18,face='bold'),
        axis.line = element_line(colour = 'grey60', size = 0.35), axis.ticks = element_line(colour = 'grey60', size = 0.35)) 
theme_set(theme_figs)

anova.tablefunc <- function(mod) {  
  dat <- data.frame(anova(mod)) %>%
    tibble::rownames_to_column(var = 'Effect') %>%
    dplyr::mutate(F.value = round(F.value, digits = 2), 
                  p.value = round(p.value, digits = 3)) %>%
    dplyr::mutate(p.value = ifelse(p.value == 0.000, '< 0.001', p.value)) %>%
    tidyr::unite(DF, numDF, denDF, sep = ',') %>%
    dplyr::rename(`F value` = F.value, `P value` = p.value)
  return(dat)
}

summary.tablefunc <- function(mod) {  
  dat <- data.frame(summary(mod)$tTable) %>%
    tibble::rownames_to_column(var = 'Effect') %>%
    rename_with(stringr::str_replace, 
                pattern = "-", replacement = ".", 
                matches("Length")) %>% 
    dplyr::mutate(Estimate = signif(Value, digits = 3),
                  Std.Error = signif(Std.Error, digits = 2),
                  t.value = signif(t.value, digits = 2),
                  p.value = signif(p.value, digits = 2)) %>%
    dplyr::mutate(p.value = ifelse(p.value <= 0.001, '< 0.001', p.value)) %>% 
    dplyr::select(-Value) %>% 
    relocate(Estimate,.before = Std.Error)
  return(dat)
}


dominant_pop <- read_csv('./Data/Dominants-through-time_2022-10-25.csv')

#dom_comm <- read_csv('./Data/RAC-subordinate_2022-07-20.csv')


dominant_pop

dominant_pop$local_lifespan <- gsub('ANNUAL','Annual',dominant_pop$local_lifespan)
dominant_pop$local_lifespan <- gsub('PERENNIAL','Perennial',dominant_pop$local_lifespan)


######## Decay of cover predicting cover over time (relative) ####
# R2out <- NULL
# 
# for(i in 1:10){
#   dom_temp <- dominant_pop %>%
#     filter(year_trt == i) # %>%
#     #filter(max_year > 10)
#   for(j in unique(dominant_pop$trt)) {
#     trt_temp <- dom_temp[dom_temp$trt == j,]
#     mod_temp_rel <- lmer(cover ~ initial_rel_cover + (1|site_code), data = trt_temp)
#     coef <- summary(mod_temp_rel)$coefficients[2,1]
#     R2_temp <- partR2(mod_temp_rel, data = dom_temp, R2_type = "marginal", nboot = 10)
#     R2_df <- data.frame(R2_temp$R2,trt = j, year_trt = i,coef = coef)
#     R2out <- rbind(R2out,R2_df)
#   }
#   }
# 
# R2out
# ## direction of effect is always positive

display.brewer.pal(11, "Blues")
brewer.pal(11, "Blues")

pos <- position_dodge(0.5)

ctrl_pal <- c('black',"#D53E4F","#6BAED6","#5E4FA2")
# 
# gg_cover_year_R2 <- ggplot(R2out, aes(x=year_trt,y=estimate,col = trt)) +
#   #geom_ribbon(aes(ymin = CI_lower,ymax= CI_upper,fill=trt),alpha=0.15,color='white') +
#   #geom_errorbar(aes(group = trt,ymin = CI_lower,ymax= CI_upper),alpha =0.5,color='black',position = pos, width = 0.3,size=1.05) +
#   geom_errorbar(aes(ymin = CI_lower,ymax= CI_upper),alpha = 0.3,position = pos, width = 0.3) +
#   geom_line(position = pos,size=1.2,alpha=0.9) +
#   geom_point(aes(group = trt),position = pos,size=1.7,color='black') +
#   geom_point(position = pos,size = 1.4) +
#   scale_color_manual(values = ctrl_pal,guide='none') +
#   #guides(color = guide_legend(title = "Treatment")) +
#   xlab('Year post-treatment') + ylab('Initial cover (pseudo-R2)')
# 
# gg_cover_year_R2
# #
# #
# #
# gg_cover_year_coef <- ggplot(R2out, aes(x=year_trt,y=coef,col = trt)) +
#   #geom_ribbon(aes(ymin = CI_lower,ymax= CI_upper,fill=trt),alpha=0.15,color='white') +
#   #geom_errorbar(aes(group = trt,ymin = CI_lower,ymax= CI_upper),alpha =0.5,color='black',position = pos, width = 0.3,size=1.05) +
#   #geom_errorbar(aes(ymin = CI_lower,ymax= CI_upper),alpha = 0.3,position = pos, width = 0.3) +
#   geom_line(position = pos,size=1.2,alpha=0.9) +
#   geom_point(aes(group = trt),position = pos,size=1.7,color='black') +
#   geom_point(position = pos,size = 1.4) +
#   scale_color_manual(values = ctrl_pal) +
#   guides(color = guide_legend(title = "Treatment")) +
#   xlab('Year post-treatment') + ylab('Initial cover (coefficient)')
# 
# gg_cover_year_coef
# 
# gg_cover_cover <- plot_grid(gg_cover_year_R2,gg_cover_year_coef, ncol = 2,rel_widths = c(1, 1.25))
# gg_cover_cover


# ggsave(paste0('Graphs/cover-v-cover-over-time_',Sys.Date(),'.png'),
#        gg_cover_cover, width = 12, height = 8, units = "in", dpi = 600)

year_hist <- dominant_pop %>%
  group_by(site_code) %>%
  mutate(max_year = max(year_trt)) %>% 
  ungroup %>% 
  distinct(site_code, .keep_all = TRUE)

tapply(year_hist$max_year,year_hist$max_year,length)

scale_col <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# scale everything for graphing purposes

dominant_pop$ppt_var_scale <- scale_col(dominant_pop$MAP_VAR_v2)
dominant_pop$ppt_scale <- scale_col(dominant_pop$MAP_v2)
dominant_pop$rich_scale <- scale_col(dominant_pop$site_richness)
dominant_pop$cov_scale <- scale_col(dominant_pop$initial_rel_cover)

#remove year 0

#dominant_pop <- dominant_pop %>% filter(!year_trt %in% c(0)) #remove pre-treatment year

unique(dominant_pop$site_code[is.na(dominant_pop$MAP_v2)])

cor(dominant_pop[,c('ppt_scale','ppt_var_scale','site_richness')])
cor(dominant_pop[,c('ppt_scale','ppt_var_scale','rich_scale','cov_scale','plotfreq')])

# ggsave(paste0('Graphs/FigS1-predictor-corr_',Sys.Date(),'.png'),
#        ggpairs(dominant_pop[,c('initial_rel_cover','site_richness','MAP_v2','MAP_VAR_v2')],
#         columnLabels = c('Initial\nRelative Cover','Site Richness','Mean Annual\nPrecipitation','CV of MAP'),
#         axisLabels="none"),
#        width = 10, height = 10)

pairs(dominant_pop[,c('initial_rel_cover','site_richness','MAP_v2','MAP_VAR_v2')], pch = 21,  cex = 1,
      bg = 'cadetblue3', col = 'black',
      labels = c('Initial\nRelative Cover','Site Richness','Mean Annual\nPrecipitation','CV of MAP'),
      lower.panel=NULL)

#### Model of rank decay ####

# asinTransform <- function(p) { asin(sqrt(p)) }

logit

# dominant_pop$rank_as <- asinTransform(dominant_pop$perc_rank)
dominant_pop$rank_adj <- ifelse(dominant_pop$perc_rank == 1, 0.99,
                                ifelse(dominant_pop$perc_rank == 0, 0.01, dominant_pop$perc_rank))
dominant_pop$rank_logit <- logit(dominant_pop$rank_adj)

#### reduce functional groups to two roughly equal (taxon-wise) groups for simplicity
dominant_pop$func_simple <- ifelse(dominant_pop$functional_group == 'GRAMINOID','GRAMINOID','NON-GRAMINOID')

dom_complete <- dominant_pop %>% filter(local_provenance != 'UNK')

mod_rank <- lme(rank_logit ~ year_trt * initial_rel_cover * NPK * Fence +
                                   local_lifespan * year_trt * NPK * Fence  +
                                   local_provenance * year_trt * NPK * Fence  +
                                   func_simple * year_trt * NPK * Fence +
                                   MAP_v2 * year_trt * NPK * Fence +
                                   MAP_VAR_v2 * year_trt * NPK * Fence +
                                   site_richness * year_trt * NPK * Fence,
                                   random = ~ 1|site_code/plot,
                                 data = dom_complete)

# mod_rank_lmer <- lmer(rank_logit ~ year_trt * initial_rel_cover * NPK * Fence +
#                   local_lifespan * year_trt * NPK * Fence  +
#                   local_provenance * year_trt * NPK * Fence  +
#                   func_simple * year_trt * NPK * Fence +
#                   MAP_v2 * year_trt * NPK * Fence +
#                   MAP_VAR_v2 * year_trt * NPK * Fence +
#                   site_richness * year_trt * NPK * Fence + 
#                     (1|site_code/plot),
#                 data = dom_complete)

# mod_rank_yr <- lme(perc_rank ~  initial_rel_cover * year_trt * NPK * Fence +
#                   local_lifespan * year_trt * NPK * Fence  +
#                   ppt_scale * year_trt * NPK * Fence +
#                   rich_scale * year_trt * NPK * Fence,
#                 random = ~ 1|site_code/plot,
#                 data = dom_pop_yr[dom_pop_yr$year_trt != 0,])
# 
summary(mod_rank)$tTable
anova(mod_rank)

rank_table <- summary.tablefunc(mod_rank)

rank_table <- rank_table %>% 
  mutate(order = if_else(Effect %in% grep('initial',rank_table$Effect, value = TRUE),2,
                   if_else(Effect %in% grep('lifespan',rank_table$Effect, value = TRUE),3,
                     if_else(Effect %in% grep('provenance',rank_table$Effect, value = TRUE),4,
                       if_else(Effect %in% grep('richness',rank_table$Effect, value = TRUE),5,
                         if_else(Effect %in% grep('func',rank_table$Effect, value = TRUE),6,
                          if_else(Effect %in% grep('MAP_v2',rank_table$Effect, value = TRUE),7,
                           if_else(Effect %in% grep('MAP_VAR_v2',rank_table$Effect, value = TRUE),8,1)))))))) %>% 
  arrange(order)

rank_table

cor(dominant_pop[,c('MAP_v2','site_richness','MAP_VAR_v2','MAT_v2','initial_rel_cover')])

r.squaredGLMM(mod_rank)
# 
#         R2m       R2c
# [1,] 0.2403656 0.5550991

# write_csv(rank_table,
#           paste0('./Tables/rank-decay-logit-table_',
#                  Sys.Date(),'.csv'))


#### graphing rank decay ####

## define function for determining graphical display of continuous gradients

# upper <- function(x) max(x)
# lower <- function(x) min(x)
upper <- function(x) quantile(x, 0.95)
lower <- function(x) quantile(x, 0.05)
# upper <- function(x) mean(x) + (sd(x)*1.5)
# lower <- function(x) mean(x) - (sd(x)*1.5)


upper(dominant_pop$initial_rel_cover)


#create discrete categories for graphing purporses
dominant_pop <- dominant_pop %>% 
  mutate(
    initial_cover  = if_else(initial_rel_cover < median(initial_rel_cover), 'Low Initial Cover','High Initial Cover'),
    richness  = if_else(site_richness < median(site_richness), 'Low site richness','High site richness'),
    ppt  = if_else(MAP_v2 < median(MAP_v2), 'Low MAP','High MAP'),
    pptvar  = if_else(MAP_VAR_v2 < median(MAP_VAR_v2), 'Low PPT Variation','High PPT Variation')
    )


#* lifespan effect --------


rank_table %>% filter(Effect %in% grep('lifespan',rank_table$Effect, value = TRUE))

### all sig

rank_life_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | local_lifespan * year_trt,
                                            at = list(year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_life_emm <- rank_life_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_table %>% filter(Effect %in% grep('lifespan',rank_table$Effect, value = TRUE))

#rank_life_emm$sig <- ifelse(rank_life_emm$trt == 'NPK+Fence','No','Yes')
rank_life_emm$sig <- 'Yes'

# rank_life_emm$em_btrans <- sin(rank_life_emm$emmean)^2
# rank_life_emm$low_btrans <- sin(rank_life_emm$lower.CL)^2
# rank_life_emm$up_btrans <- sin(rank_life_emm$upper.CL)^2

rank_life_emm$em_btrans <- invlogit(rank_life_emm$emmean)
rank_life_emm$low_btrans <- invlogit(rank_life_emm$lower.CL)
rank_life_emm$up_btrans <- invlogit(rank_life_emm$upper.CL)


gg_rank_year_life_trt <- ggplot(rank_life_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.2, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.15,color=NA) +
  geom_line(aes(group = trt), size = 2) +
  #geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  #scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  facet_grid(~local_lifespan)
gg_rank_year_life_trt

# ggsave(paste0('Graphs/Presentation/rank-decay-lifespan_',Sys.Date(),'.png'),
#        gg_rank_year_life_trt, width = 12, height = 8, units = "in", dpi = 600)


#* provenance effect --------

rank_table %>% filter(Effect %in% grep('provenance',rank_table$Effect, value = TRUE))
### NPK+fence

rank_prov_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | local_provenance * year_trt,
                                            at = list(year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_prov_emm <- rank_prov_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_prov_emm$sig <- ifelse(rank_life_emm$trt %in% c('NPK+Fence'),'No','Yes')

# rank_life_emm$em_btrans <- sin(rank_life_emm$emmean)^2
# rank_life_emm$low_btrans <- sin(rank_life_emm$lower.CL)^2
# rank_life_emm$up_btrans <- sin(rank_life_emm$upper.CL)^2

rank_prov_emm$em_btrans <- invlogit(rank_prov_emm$emmean)
rank_prov_emm$low_btrans <- invlogit(rank_prov_emm$lower.CL)
rank_prov_emm$up_btrans <- invlogit(rank_prov_emm$upper.CL)


gg_rank_year_prov_trt <- ggplot(rank_prov_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dom_complete, aes(x=year_trt, y = perc_rank),  width=0.2, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.15,color=NA) +
  #geom_line(aes(group = trt), size = 2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  facet_grid(~local_provenance)
gg_rank_year_prov_trt


#* functional group effect --------


rank_table %>% filter(Effect %in% grep('func',rank_table$Effect, value = TRUE))
### only fence is sig (npk+fence marginal)

rank_func_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | func_simple * year_trt,
                                            at = list(year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_func_emm <- rank_func_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_func_emm$sig <- ifelse(rank_func_emm$trt %in% c('NPK+Fence','NPK'),'No','Yes')

# rank_life_emm$em_btrans <- sin(rank_life_emm$emmean)^2
# rank_life_emm$low_btrans <- sin(rank_life_emm$lower.CL)^2
# rank_life_emm$up_btrans <- sin(rank_life_emm$upper.CL)^2

rank_func_emm$em_btrans <- invlogit(rank_func_emm$emmean)
rank_func_emm$low_btrans <- invlogit(rank_func_emm$lower.CL)
rank_func_emm$up_btrans <- invlogit(rank_func_emm$upper.CL)


gg_rank_year_func_trt <- ggplot(rank_func_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dom_complete, aes(x=year_trt, y = perc_rank),  width=0.2, alpha = 0.8, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.15,color=NA) +
  #geom_line(aes(group = trt), size = 2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  facet_grid(~func_simple)
gg_rank_year_func_trt


#* initial cover effect ----

rank_table %>% filter(Effect %in% grep('initial',rank_table$Effect, value = TRUE))
### all sig

# cov.low <- mean(dominant_pop$initial_rel_cover) - 2 * sd(dominant_pop$initial_rel_cover)
# cov.high <- mean(dominant_pop$initial_rel_cover) + 2 * sd(dominant_pop$initial_rel_cover)

cov.low <- lower(dominant_pop$initial_rel_cover)
cov.high <- upper(dominant_pop$initial_rel_cover)


rank_init_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | initial_rel_cover * year_trt,
                                            at = list(initial_rel_cover = seq(from = cov.low, to = cov.high, by =cov.high-cov.low),
                                                      year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_init_emm <- rank_init_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_init_emm$initial_cover <- factor(ifelse(rank_init_emm$initial_rel_cover == cov.low, 'Low Initial Cover','High Initial Cover'), levels = c('Low Initial Cover','High Initial Cover'))
#rank_init_emm$upper.CL <- ifelse(rank_init_emm$upper.CL > 1.05, 1.05, rank_init_emm$upper.CL)
#rank_init_emm$lower.CL <- ifelse(rank_init_emm$lower.CL < 0, 0, rank_init_emm$lower.CL)

# rank_init_emm$em_btrans <- sin(rank_init_emm$emmean)^2
# rank_init_emm$low_btrans <- sin(rank_init_emm$lower.CL)^2
# rank_init_emm$up_btrans <- sin(rank_init_emm$upper.CL)^2
#### should negative arcsin values be set to 0 on backtransformation (otherwise artificial uptick due to sin)?

rank_init_emm$em_btrans <- invlogit(rank_init_emm$emmean)
rank_init_emm$low_btrans <- invlogit(rank_init_emm$lower.CL)
rank_init_emm$up_btrans <- invlogit(rank_init_emm$upper.CL)



gg_rank_year_init_trt <- ggplot(rank_init_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.11, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.25,color=NA) +
  geom_line(aes(group = trt), size = 2, color= 'black') +
  geom_line(aes(group = trt), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  #ylim(0,1.05) +
  facet_wrap(~initial_cover)

gg_rank_year_init_trt

# ggsave(paste0('Graphs/Presentation/rank-decay-initcover_',Sys.Date(),'.png'),
#        gg_rank_year_init_trt, width = 12, height = 8, units = "in", dpi = 600)

#* richness effect -----

rank_table %>% filter(Effect %in% grep('richness',rank_table$Effect, value = TRUE))

### only npk sig

rich.low <- lower(dominant_pop$site_richness)
rich.high <- upper(dominant_pop$site_richness)


rank_rich_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | site_richness * year_trt,
                                            at = list(site_richness = seq(from = rich.low, to = rich.high, by = (rich.high-rich.low)),
                                                      year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_rich_emm <- rank_rich_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_rich_emm$sig <- ifelse(rank_rich_emm$trt %in% c('Fence','NPK+Fence'),'No','Yes')


rank_rich_emm$richness <- factor(ifelse(rank_rich_emm$site_richness == rich.low, 'Low site richness','High site richness'), levels = c('Low site richness','High site richness'))
# rank_rich_emm$upper.CL <- ifelse(rank_rich_emm$upper.CL > 1.05, 1.05, rank_rich_emm$upper.CL)
# rank_rich_emm$lower.CL <- ifelse(rank_rich_emm$lower.CL < 0, 0, rank_rich_emm$lower.CL)

# rank_rich_emm$em_btrans <- sin(rank_rich_emm$emmean)^2
# rank_rich_emm$low_btrans <- sin(rank_rich_emm$lower.CL)^2
# rank_rich_emm$up_btrans <- sin(rank_rich_emm$upper.CL)^2

rank_rich_emm <-  mutate(rank_rich_emm,
                          em_btrans = invlogit(emmean), 
                          low_btrans = invlogit(lower.CL),
                          up_btrans = invlogit(upper.CL))


gg_rank_year_rich_trt <- ggplot(rank_rich_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.11, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.25,color=NA) +
  #geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  ylim(0,1.05) +
  facet_wrap(~richness)

gg_rank_year_rich_trt

# ggsave(paste0('Graphs/Presentation/rank-decay-initcover_',Sys.Date(),'.png'),
#        gg_rank_year_init_trt, width = 12, height = 8, units = "in", dpi = 600)


#* MAP effect ----

rank_table %>% filter(Effect %in% grep('MAP_v2',rank_table$Effect, value = TRUE))

#fence sig

map.low <- lower(dominant_pop$MAP_v2)
map.high <- upper(dominant_pop$MAP_v2)


rank_map_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | MAP_v2 * year_trt,
                                            at = list(MAP_v2 = seq(from = map.low, to = map.high, by = map.high-map.low),
                                                      year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_map_emm <- rank_map_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))


rank_map_emm$sig <- ifelse(rank_map_emm$trt %in% c('NPK','NPK+Fence'),'No','Yes')


rank_map_emm$ppt <- factor(ifelse(rank_map_emm$MAP_v2 == map.low, 'Low MAP','High MAP'), levels = c('Low MAP','High MAP'))
# rank_map_emm$upper.CL <- ifelse(rank_map_emm$upper.CL > 1.05, 1.05, rank_map_emm$upper.CL)
#rank_map_emm$lower.CL <- ifelse(rank_map_emm$lower.CL < 0, 0, rank_map_emm$lower.CL)

# rank_map_emm$em_btrans <- sin(rank_map_emm$emmean)^2
# rank_map_emm$low_btrans <- sin(rank_map_emm$lower.CL)^2
# rank_map_emm$up_btrans <- sin(rank_map_emm$upper.CL)^2

rank_map_emm <-  mutate(rank_map_emm,
                         em_btrans = invlogit(emmean), 
                         low_btrans = invlogit(lower.CL),
                         up_btrans = invlogit(upper.CL))

gg_rank_year_ppt_trt <- ggplot(rank_map_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.11, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.25,color=NA) +
  #geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  ylim(0,1.05) +
  facet_wrap(~ppt)

gg_rank_year_ppt_trt


#* MAP VAR effect ----

rank_table %>% filter(Effect %in% grep('MAP_VAR_v2',rank_table$Effect, value = TRUE))
# fence is sig, NPK+fence marginal

mvar.low <- lower(dominant_pop$MAP_VAR_v2)
mvar.high <- upper(dominant_pop$MAP_VAR_v2)


rank_pptvar_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | MAP_VAR_v2 * year_trt,
                                           at = list(MAP_VAR_v2 = seq(from = mvar.low, to = mvar.high, by = mvar.high-mvar.low),
                                                     year_trt = seq(from = 0, to = 15, by = 1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))

rank_pptvar_emm <- rank_pptvar_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_pptvar_emm$sig <- ifelse(rank_pptvar_emm$trt %in% c('NPK','NPK+Fence'),'No','Yes')


rank_pptvar_emm$pptvar <- factor(ifelse(rank_pptvar_emm$MAP_VAR_v2 == mvar.low, 'Low PPT Variation','High PPT Variation'), levels = c('Low PPT Variation','High PPT Variation'))
# rank_pptvar_emm$upper.CL <- ifelse(rank_pptvar_emm$upper.CL > 1.05, 1.05, rank_pptvar_emm$upper.CL)
# rank_pptvar_emm$lower.CL <- ifelse(rank_pptvar_emm$lower.CL < 0, 0, rank_pptvar_emm$lower.CL)

# rank_pptvar_emm$em_btrans <- sin(rank_pptvar_emm$emmean)^2
# rank_pptvar_emm$low_btrans <- sin(rank_pptvar_emm$lower.CL)^2
# rank_pptvar_emm$up_btrans <- sin(rank_pptvar_emm$upper.CL)^2

rank_pptvar_emm <-  mutate(rank_pptvar_emm,
                        em_btrans = invlogit(emmean), 
                        low_btrans = invlogit(lower.CL),
                        up_btrans = invlogit(upper.CL))


gg_rank_year_pptvar_trt <- ggplot(rank_pptvar_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.11, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.35,color=NA) +
  #geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  ylim(0,1.05) +
  facet_wrap(~pptvar)

gg_rank_year_pptvar_trt

#### Fig 2 ####
#### stitch them together

rank_legend <- get_legend(gg_rank_year_func_trt + 
                            theme(legend.direction = "horizontal", legend.box = "horizontal"))

gg_rank_init <- gg_rank_year_init_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
                                             
gg_rank_life <- gg_rank_year_life_trt + theme(legend.position="none", axis.title.y=element_blank())
                                              #axis.title.x=element_blank(),axis.text.x=element_blank())

# Fig 2 with initial rank and lifespan
# Fig2 <- ggdraw() +
#   draw_plot(rank_legend, x = -0.3, y = 0.95, width = 1, height = .05) +
#   draw_plot(gg_rank_init, x = 0.03, y = 0.52, width = .97, height = .43) +
#   draw_plot(gg_rank_life, x = 0.03, y = 0, width = .97, height = .5) +
#   draw_plot_label(label = c("A", "B"), size = 20,
#                     x = c(0.03, 0.03), y = c(0.97, 0.51)) +
#   draw_text('Rank Percentile',angle = 90, x = 0.015, size = 25)

#Fig 2 just initial cover

Fig2 <- gg_rank_year_init_trt

Fig2


# Fig 3 - species chars

gg_rank_life <- gg_rank_year_life_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
gg_rank_func <- gg_rank_year_func_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
gg_rank_prov <- gg_rank_year_prov_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())

#axis.title.x=element_blank(),axis.text.x=element_blank())

Fig3 <- ggdraw() +
  draw_plot(rank_legend, x = -0.3, y = 0.95, width = 1, height = .05) +
  draw_plot(gg_rank_life, x = 0.03, y = 0.65, width = .97, height = .3) +
  draw_plot(gg_rank_func, x = 0.03, y = 0.35, width = .97, height = .3) +
  draw_plot(gg_rank_prov, x = 0.03, y = 0.05, width = .97, height = .3) +
  draw_plot_label(label = c("A", "B", "C"), size = 20,
                    x = c(0.03, 0.03, 0.03), y = c(0.95, 0.65, 0.35)) +
  draw_text('Rank Percentile',angle = 90, x = 0.015, size = 25) +
  draw_text('Year after treatment', x = 0.35, y= 0.02, size = 25)

Fig3

# Fig 4 - species chars

gg_rank_rich <- gg_rank_year_rich_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
gg_rank_ppt <- gg_rank_year_ppt_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
gg_rank_pptvar <- gg_rank_year_pptvar_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())

#axis.title.x=element_blank(),axis.text.x=element_blank())

Fig4 <- ggdraw() +
  draw_plot(rank_legend, x = -0.3, y = 0.95, width = 1, height = .05) +
  draw_plot(gg_rank_rich, x = 0.03, y = 0.65, width = .97, height = .3) +
  draw_plot(gg_rank_ppt, x = 0.03, y = 0.35, width = .97, height = .3) +
  draw_plot(gg_rank_pptvar, x = 0.03, y = 0.05, width = .97, height = .3) +
  draw_plot_label(label = c("A", "B", "C"), size = 20,
                  x = c(0.03, 0.03, 0.03), y = c(0.95, 0.65, 0.35)) +
  draw_text('Rank Percentile',angle = 90, x = 0.015, size = 25) +
  draw_text('Year after treatment', x = 0.35, y= 0.02, size = 25)

Fig4

# gg_rank_freq <- gg_rank_year_life_trt + theme(legend.position="none")
gg_rank_rich <- gg_rank_year_rich_trt + 
  theme(legend.position="none",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank()) 
# + 
  #labs(tag = 'C')
gg_rank_ppt <- gg_rank_year_ppt_trt + theme(legend.position="none",axis.title.y=element_blank(),axis.title.x=element_blank())
                                            # plot.tag.position = c(0.08,1)) 
# + 
#   
#   labs(tag = 'D')
gg_rank_pptvar <- gg_rank_year_pptvar_trt + theme(legend.position="none",axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank()) 
  #labs(tag = 'E')

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

Fig2.rank <- grid.arrange(gg_rank_init,blankPlot,
                          rank_legend,
                          gg_rank_life,blankPlot,
                          gg_rank_rich,
                          gg_rank_ppt, blankPlot,
                          gg_rank_pptvar,
             ncol = 3,nrow=3,
             widths = c(3,0.1,2.5), heights = c(3,3,3.2),
             left = textGrob("Rank percentage", rot = 90, gp=gpar(fontsize=24)),
             bottom = textGrob("Year after treatment", gp=gpar(fontsize=24)))


# ggdraw() +
#   draw_plot(gg_rank_init, x = 0, y = 0, width = .5, height = .5) +
#   draw_plot(gg_rank_life, x = .5, y = 0.25, width = .25, height = .25) +
#   draw_plot(gg_rank_rich, x = .75, y = 0.25, width = .25, height = 0.25) +
#   draw_plot(gg_rank_ppt, x = .5, y = 0, width = .25, height = 0.25) +
#   draw_plot(gg_rank_pptvar, x = .75, y = 0, width = .25, height = 0.25) 
# 
# 
# Fig2 <- ggdraw() +
#   draw_plot(rank_legend, x = -0.3, y = 0.95, width = 1, height = .05) +
#   draw_plot(gg_rank_init, x = 0.03, y = 0.52, width = .97, height = .43) +
#   draw_plot(gg_rank_life, x = 0.03, y = 0.27, width = .5, height = .25) +
#   draw_plot(gg_rank_rich, x = .53, y = 0.27, width = .47, height = 0.25) +
#   draw_plot(gg_rank_ppt, x = 0.03, y = 0.02, width = .5, height = 0.25) +
#   draw_plot(gg_rank_pptvar, x = .53, y = 0.02, width = .47, height = 0.25)+
#   draw_plot_label(label = c("A", "B", "C", "D",'E'), size = 15,
#                   x = c(0.1, 0.1, 0.54,0.1,0.54), y = c(0.89, 0.46, 0.46,0.21,0.21)) + 
#   draw_label("Years after treatment", colour = "black", size = 25, x = 0.5, y = 0.02) +
#   draw_label("Rank Percentile", colour = "black", size = 25, angle = 90, x = 0.01, y = 0.5)




# save_plot(paste0('Graphs/Fig2-rank-time_',Sys.Date(),'.png'),
#           Fig2,
#           base_height = 8, base_width = 12)

# AS PDF
# save_plot(paste0('Graphs/Fig2-rank-time_',Sys.Date(),'.pdf'),
#           Fig2,
#           base_height = 8, base_width = 12)


# Defining inertia #####

### first can we create an 'effect size' table of when model upper CI drop below 0.95 as a proxy for inertia?
### on closer examination - looking at the extreme values of each continuous variable means wider confidence intervals, 

emm_options(rg.limit = 200000)

cf = mean(dom_complete$perc_rank[dom_complete$year_trt != 0])

inertia_npk <-  data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence |  year_trt,
                                                          at = list(year_trt = seq(from = 0, to = 15, by = 0.1)),
                                                          cov.reduce = TRUE,
                                                          type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL),Predictor = 'Average') %>% 
  group_by(NPK,Fence) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = Predictor) 

inertia_npk

inertia_init <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | initial_rel_cover * year_trt,
#                                           at = list(initial_rel_cover = seq(from = cov.low, to = cov.high, by = (cov.high-cov.low)/2),
                                            at = list(initial_rel_cover = seq(from = cov.low, to = cov.high, by = (cov.high-cov.low)),
                                                                         year_trt = seq(from = 0, to = 15, by = 0.1)),
                                                               cov.reduce = TRUE,
                                                               type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, initial_rel_cover) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  mutate(initial_cover = factor(if_else(initial_rel_cover == cov.low, 'Low Initial Cover',
                                        if_else(initial_rel_cover == cov.high,'High Initial Cover', 'Average Initial Cover')))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = initial_cover) %>% 
  rename(Predictor = 1) 

inertia_init

inertia_map <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | MAP_v2 * year_trt,
                                           # at = list(MAP_v2 = seq(from = map.low, to = map.high, by = (map.high-map.low)/2),
                                                     at = list(MAP_v2 = seq(from = map.low, to = map.high, by = (map.high-map.low)),
                                                     year_trt = seq(from = 0, to = 15, by = 0.1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, MAP_v2) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  mutate(MAP = factor(if_else(MAP_v2 == map.low, 'Low MAP',
                                        if_else(MAP_v2 == map.high,'High MAP', 'Average MAP')))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = MAP) %>% 
  rename(Predictor = 1)

inertia_map

inertia_mapvar <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | MAP_VAR_v2 * year_trt,
                                          # at = list(MAP_VAR_v2 = seq(from = mvar.low, to = mvar.high, by = (mvar.high-mvar.low)/2),
                                                    at = list(MAP_VAR_v2 = seq(from = mvar.low, to = mvar.high, by = (mvar.high-mvar.low)),
                                                    year_trt = seq(from = 0, to = 15, by = 0.1)),
                                          cov.reduce = TRUE,
                                          type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, MAP_VAR_v2) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  mutate(ppt_variation = factor(if_else(MAP_VAR_v2 == mvar.low, 'Low PPT Variation',
                              if_else(MAP_VAR_v2 == mvar.high,'High PPT Variation', 'Average PPT Variation')))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = ppt_variation)%>% 
  rename(Predictor = 1)

inertia_mapvar

inertia_rich <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | site_richness * year_trt,
                                             # at = list(site_richness = seq(from = rich.low, to = rich.high, by = (rich.high-rich.low)/2),
                                                       at = list(site_richness = seq(from = rich.low, to = rich.high, by = (rich.high-rich.low)),
                                                       year_trt = seq(from = 0, to = 15, by = 0.1)),
                                             cov.reduce = TRUE,
                                             type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, site_richness) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  mutate(Site_Richness = factor(if_else(site_richness == rich.low, 'Low Site Richness',
                                        if_else(site_richness == rich.high,'High Site Richness', 'Average Site Richness')))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = Site_Richness)%>% 
  rename(Predictor = 1)

inertia_rich

inertia_life <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | local_lifespan * year_trt,
                                           at = list(year_trt = seq(from = 0, to = 15, by = 0.1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, local_lifespan) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = local_lifespan)%>% 
  rename(Predictor = 1)

inertia_life  

inertia_func <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | func_simple * year_trt,
                                           at = list(year_trt = seq(from = 0, to = 15, by = 0.1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, func_simple) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = func_simple)%>% 
  rename(Predictor = 1)

inertia_func
  
inertia_prov <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | local_provenance * year_trt,
                                           at = list(year_trt = seq(from = 0, to = 15, by = 0.1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, local_provenance) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = local_provenance) %>% 
  rename(Predictor = 1)

inertia_prov

inertia_tab <- bind_rows(inertia_npk, inertia_init,inertia_life,inertia_func,inertia_prov,inertia_map,inertia_mapvar,inertia_rich)

inertia_tab <-inertia_tab %>% mutate(fence_per = -100*(1-Fence/Control),
                       npk_per = -100*(1-NPK/Control),
                       npkf_per = -100*(1-`NPK+Fence`/Control)) %>% 
  mutate(Fence = paste0(Fence, ' (',signif(fence_per,3),'%)'),
         NPK = paste0(NPK, ' (',signif(npk_per,3),'%)'),
         `NPK+Fence` = paste0(`NPK+Fence`, ' (',signif(npkf_per,3),'%)')) %>% 
  dplyr::select(-c(fence_per,npk_per,npkf_per))
  
inertia_tab

#write_csv(inertia_tab,'./Tables/Inertia-effect-size.csv')

####

dom_pop_yr_pred  <- dom_complete %>%
  bind_cols(.,predictInterval(mod_rank_lmer,dom_complete,
                              level = 0.90, n.sims = 1000,
                              stat = "median", type="linear.prediction",
                              include.resid.var = TRUE))


# only works for lmer models

dom_complete$pred_logit <- predict(mod_rank,dom_complete)
dom_complete$pred_invlogit <- invlogit(dom_complete$pred_logit)

#### Figure 2 Alt - individual linear models ####

## define year span for rank linear rank loss
yr_cut <- 5

dom_yr <- dominant_pop %>%
  group_by(site_code) %>%
  mutate(max_year = max(year_trt)) %>%
  ungroup() %>%
  filter(max_year >= yr_cut)

length(unique(dom_yr$site_code))

rank_loss_yr <- split(dom_yr, ~ site_code + plot) %>% 
  discard(~nrow(.)==0) %>%
  map(~lm(perc_rank-1 ~ year_trt-1, data = .,na.action=na.omit))  %>%  #one way to fix intercept at '1'
  map_df(tidy,.id='plot') %>%
  rename(est_rank = estimate)
  

rank_loss_yr

#rank_loss_yr %>% filter(site_code == 'burrawan.au' & plot == 1)

rank_loss_yr_graph <- rank_loss_yr %>% 
  left_join(.,dominant_pop %>% dplyr::select(site_code, plot, Taxon, initial_rel_cover, perc_rank, year_trt, trt) %>% 
              mutate(id = paste0(site_code,'.',plot)),
            by = c('plot'='id')) %>% 
  mutate(pred = 1 + (year_trt*est_rank))
              
rank_loss_yr_graph

rank_loss_yr_stats <- rank_loss_yr %>% 
  left_join(.,dominant_pop %>% 
              dplyr::select(site_code, plot, Taxon, initial_rel_cover, MAP_v2,MAP_VAR_v2,site_richness,local_lifespan,plotfreq, trt, NPK, Fence) %>% 
              mutate(id = paste0(site_code,'.',plot)) %>% 
              distinct(id,.keep_all = TRUE),
            by = c('plot'='id'))

lme_yr_rank_slope <- lme(est_rank ~ NPK*Fence*initial_rel_cover +
                           NPK*Fence*local_lifespan +
                           plotfreq * NPK * Fence +
                           MAP_v2 * NPK * Fence +
                           MAP_VAR_v2 * NPK * Fence +
                           site_richness * NPK * Fence,
                         random = ~1|site_code/Taxon,
                         data = rank_loss_yr_stats)

summary(lme_yr_rank_slope)
anova(lme_yr_rank_slope)

## relative rank effect

cov.low <- mean(dominant_pop$initial_rel_cover) - 1.5 * sd(dominant_pop$initial_rel_cover)
cov.high <- mean(dominant_pop$initial_rel_cover) + 1.5 * sd(dominant_pop$initial_rel_cover)

yr_init_emm <- data.frame(summary(emmeans(lme_yr_rank_slope, pairwise ~ NPK * Fence | initial_rel_cover,
                                            at = list(initial_rel_cover = seq(from = 0.2, to = 1, by = 0.8)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

yr_init_emm <- yr_init_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

yr_init_emm

yr_init_emm$initial_cover <- factor(ifelse(yr_init_emm$initial_rel_cover == 0.2, 'Low Initial Cover','High Initial Cover'), levels = c('Low Initial Cover','High Initial Cover'))
#yr_init_emm$upper.CL <- ifelse(yr_init_emm$upper.CL > 1.05, 1.05, yr_init_emm$upper.CL)
#yr_init_emm$lower.CL <- ifelse(yr_init_emm$lower.CL < 0, 0, yr_init_emm$lower.CL)

range(rank_loss_yr_graph$initial_rel_cover[rank_loss_yr_graph$year_trt == 0])
mean(rank_loss_yr_graph$initial_rel_cover)
rank_loss_yr_graph$initial_cover <- ifelse(rank_loss_yr_graph$initial_rel_cover < mean(rank_loss_yr_graph$initial_rel_cover),
                                           'Low Initial Cover','High Initial Cover')


rank_loss_yr_graph$initial_cover <- factor(rank_loss_yr_graph$initial_cover, levels = c('Low Initial Cover','High Initial Cover'))

gg_rank_year_init_trt_alt <- ggplot(rank_loss_yr_graph, aes(x = year_trt, y = perc_rank, col = trt)) +
  geom_jitter(width=0.1, alpha = 0.8, col = 'grey90') +
  #geom_ribbon(aes(ymin = lower.CL+1,ymax= upper.CL+1,fill=trt),alpha=0.25,color=NA) +
  geom_abline(data = yr_init_emm, aes(intercept = 1, slope = emmean, col = trt), size = 2) +
  #geom_line(aes(group = trt), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  ylim(0,1.05) +
  facet_wrap(~initial_cover)

gg_rank_year_init_trt_alt

#### illustration of linear fit problem

ggplot(dominant_pop %>% filter(site_code == 'arch.us' & plot == 28 & Taxon == 'AXONOPUS FURCATUS'),
       aes(x=year_trt,y = perc_rank - 1)) + geom_point() + geom_line() + 
      geom_smooth(method = 'lm', se = FALSE, formula=y~x+0) + ylim(-1,0)
ggplot(dominant_pop %>% filter(site_code == 'hall.us' & plot == 3 & Taxon == 'ANDROPOGON GERARDII'),
       aes(x=year_trt,y = perc_rank - 1)) + geom_point() + geom_line() + 
  geom_smooth(method = 'lm', se = FALSE, formula=y~x+0) + ylim(-1,0)
ggplot(dominant_pop %>% filter(site_code == 'azi.cn' & plot == 18 & Taxon == 'AGROSTIS STOLONIFERA'),
       aes(x=year_trt,y = perc_rank - 1)) + geom_point() + geom_line() + 
  geom_smooth(method = 'lm', se = FALSE, formula=y~x+0) + ylim(-1,0)



## defining inertia (logit model) ####

## define year span for inertia component
yr_cut <- 5

dom_pop_yr <- dominant_pop %>%
  group_by(site_code) %>%
  mutate(max_year = max(year_trt)) %>%
  ungroup() %>%
  filter(max_year >= yr_cut) 

# %>% 
#   filter(year_trt <= yr_cut)

length(unique(dom_pop_yr$site_code))
#25 for 10 year cutoff
#57 for 5 year cutoff

### fit logit model to all species across time and extract first year where confidence interval doesn't include 1
# 
inertia_mod <- glmer(perc_rank ~  year_trt +
                                      (1|site_code/plot),
                                    data = dom_pop_yr[dom_pop_yr$year_trt != 0,],
                                    family = binomial, control = glmerControl(optimizer = "bobyqa",
                                                                              optCtrl=list(maxfun=2e6)))


dom_pop_yr_pred  <- dom_pop_yr[dom_pop_yr$year_trt != 0,] %>%
  bind_cols(.,predictInterval(inertia_mod,dom_pop_yr[dom_pop_yr$year_trt != 0,],
                level = 0.90, n.sims = 1000,
                stat = "median", type="probability",
                include.resid.var = TRUE))


hist(dom_pop_yr_pred$upr, nclass=100)


ggplot(dom_pop_yr_pred, aes(x=perc_rank,y=upr,col=year_trt)) + geom_point() + geom_hline(yintercept = 0.99) +
  xlab('Rank percentile') + ylab('Upper CI (logit model estimate)')

# first try - upper bound is arbitrary still :(
inertia_resistance <- dom_pop_yr_pred %>%
  group_by(site_code,plot) %>%
  arrange(site_code,plot,year_trt) %>%
  slice(if(any(upr < 0.99)) which.max(upr < 0.99) else which.max(year_trt == max(year_trt))) %>%
  mutate(inertia_mod = if_else(perc_rank == 1, year_trt, year_trt-1)) %>%
  left_join(dom_pop_yr_pred %>%
              group_by(site_code,plot) %>%
              summarize(resistance = min(perc_rank)),
            by = c('site_code','plot'))

hist(inertia_resistance$inertia_mod, xlab = 'Inertia: Year to rank loss')
hist(inertia_resistance$resistance, xlab = 'Resistance')

cor(inertia_resistance[,c('inertia_mod','resistance')])

plot(inertia_mod~resistance,data=inertia_resistance)


### defining inertia (empirical) ####

# first year of rank loss and last year of being 1
inertia_rank_loss <- dom_pop_yr %>% 
  group_by(site_code,plot) %>% 
  arrange(site_code,plot,year_trt) %>% 
  mutate(inertia_last = if_else(any(perc_rank == 1),year_trt[max(which(perc_rank==1))],0)) %>% 
  slice(if(any(perc_rank != 1)) which.max(perc_rank != 1) else which.max(year_trt == max(year_trt))) %>% 
  mutate(inertia_first = if_else(perc_rank == 1, year_trt, year_trt-1)) %>% 
  mutate(inertia = (inertia_last + inertia_first)/2)

ggplot(inertia_rank_loss, aes(x= inertia_first, y=inertia_last)) + geom_jitter()
### many species that lose rank in year 1 recover to have rank = 1 in later years.


hist(inertia_rank_loss$inertia_last)
hist(inertia_rank_loss$inertia_first)
hist(inertia_rank_loss$inertia)


### defining resistance ####

#### lowest rank reached
# rank_loss_components <- dom_pop_yr %>% 
#   group_by(site_code,plot,Taxon) %>% 
#   filter(perc_rank == min(perc_rank)) %>% 
#   distinct(site_code,plot,Taxon,.keep_all = TRUE) %>% 
#   mutate(resistance = perc_rank) %>% 
#   select(site_code,plot,Taxon, resistance) %>% 
#   right_join(inertia_rank_loss,by=c('site_code','plot','Taxon'))

#### average of two lowest years (avoid outliers)

rank_loss_components <- dom_pop_yr %>% 
  group_by(site_code,plot,Taxon) %>%
  slice_min(perc_rank,n=2) %>% 
  summarize(resistance = mean(perc_rank)) %>% 
  right_join(inertia_rank_loss,by=c('site_code','plot','Taxon'))

## how does logit model defined inertia compare to empirically defined inertia

iner_res <- inertia_resistance %>% 
  dplyr::select(site_code, plot, Taxon, inertia_mod, inertia_mod, resistance) %>% 
  left_join(rank_loss_components, by = c('site_code','plot','Taxon'))

cor(iner_res[!is.na(iner_res$inertia_first),c('resistance.x','resistance.y','inertia_mod','inertia_first','inertia_last')])

#                 resistance.x resistance.y inertia_mod inertia_first inertia_last
# resistance.x     1.0000000    0.9625569   0.9099374     0.6524125    0.6628405
# resistance.y     0.9625569    1.0000000   0.9028225     0.6179325    0.7087079
# inertia_mod      0.9099374    0.9028225   1.0000000     0.7279518    0.7352506
# inertia_first    0.6524125    0.6179325   0.7279518     1.0000000    0.5054031
# inertia_last     0.6628405    0.7087079   0.7352506     0.5054031    1.0000000

# So the logit model 'peeling away' from 1 is more linked to the minimum rank than time to rank loss. While this appears to best capture both aspects, the overall gain is minimal and the approach is convoluted, so best to drop.

#### correlation of inertia and resistance

summary(lm(rank_loss_components$inertia_last ~ rank_loss_components$resistance)) 
##R2 w 10 yrs = 0.48
summary(lm(rank_loss_components$inertia_first ~ rank_loss_components$resistance)) 
##R2 = 0.38
summary(lm(rank_loss_components$inertia ~ rank_loss_components$resistance)) 
##R2 = 0.58

plot(rank_loss_components$inertia_last ~ rank_loss_components$resistance, pch = 16)
cor(rank_loss_components$inertia_last, rank_loss_components$resistance,method = 'pearson') #0.71


gg_res_iner_corr <- ggplot(rank_loss_components,aes(x = resistance, y = inertia_last)) + 
  geom_jitter(alpha = 0.4) +
  geom_smooth(method = 'lm', se= FALSE, size = 1.5) +
  scale_color_manual(values = ctrl_pal) + scale_fill_manual(values = ctrl_pal) +
  guides(color = guide_legend(title = "Treatment"),fill = guide_legend(title = "Treatment")) + 
  ylab('Biological inertia\n(% of years until rank loss)') + xlab('Resistance')
gg_res_iner_corr

# definitely correlated, but some distinct information here

# ggsave(paste0('Graphs/resistance-v-inertia_',Sys.Date(),'.png'),
#        gg_res_iner_corr, width = 12, height = 8, units = "in", dpi = 600)


#### defining oscillations ####

osc <- split(dom_yr, ~ site_code + plot + Taxon) %>% 
  discard(~nrow(.)==0) %>%
  discard(~mean(.$perc_rank) == 1) %>% 
  map(~arima(.$perc_rank, c(1,1,0),include.mean=FALSE,method="ML"))  %>%  #### (AR1 - (1,0,0) causes error, not sure if adding MA1 is valid though)
  map_df(tidy,.id='plot') %>%
  filter(term == 'ar1')

# dom_yr <- dom_yr %>% mutate(tax_id = paste(site_code, plot, Taxon, sep = '.'))
# 
# for(i in unique(dom_yr$tax_id)) {
#   df <- dom_yr[dom_yr$tax_id == i,]
#   #if(i == "kilp.fi.17.FESTUCA OVINA") next
#   if(mean(df$perc_rank) == 1) next
#   arima(as.numeric(df$perc_rank), c(1,1,0), include.mean=FALSE, method = 'ML')
#   }

range(osc$estimate)

## high phi
osc %>% arrange(desc(estimate))

par(mfrow = c (2,3))
plot(dom_yr$perc_rank[dom_yr$site_code == 'arch.us' & dom_yr$plot == 28 & dom_yr$Taxon == 'AXONOPUS FURCATUS'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt') 
plot(dom_yr$perc_rank[dom_yr$site_code == 'frue.ch' & dom_yr$plot == 9], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt') 
plot(dom_yr$perc_rank[dom_yr$site_code == 'kbs.us' & dom_yr$plot == 20], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt') 
plot(dom_yr$perc_rank[dom_yr$site_code == 'rook.uk' & dom_yr$plot == 2 & dom_yr$Taxon == 'GALIUM SAXATILE'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt') 
plot(dom_yr$perc_rank[dom_yr$site_code == 'hall.us' & dom_yr$plot == 3 & dom_yr$Taxon == 'ANDROPOGON GERARDII'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt') 
plot(dom_yr$perc_rank[dom_yr$site_code == 'ping.au' & dom_yr$plot == 19 & dom_yr$Taxon == 'AVENA BARBATA'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt') 



#low phi
osc %>% arrange(estimate)
plot(dom_yr$perc_rank[dom_yr$site_code == 'azi.cn' & dom_yr$plot == 18 & dom_yr$Taxon == 'AGROSTIS STOLONIFERA'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'shps.us' & dom_yr$plot == 13 & dom_yr$Taxon == 'ELYMUS SPICATUS'], ylim = c(0.9,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'rook.uk' & dom_yr$plot == 24 & dom_yr$Taxon == 'GALIUM SAXATILE'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'kibber.in' & dom_yr$plot == 21 & dom_yr$Taxon == 'POLYGONUM SP.'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'azi.cn' & dom_yr$plot == 8 & dom_yr$Taxon == 'POTENTILLA ANSERINA'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'cdpt.us' & dom_yr$plot == 41 & dom_yr$Taxon == 'CAREX FILIFOLIA'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')


# median phi
osc[abs(osc$estimate) < .05,] %>% arrange(abs(estimate))
par(mfrow = c(2,3))
plot(dom_yr$perc_rank[dom_yr$site_code == 'doane.us' & dom_yr$plot == 9], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'hero.uk' & dom_yr$plot == 3], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'konz.us' & dom_yr$plot == 29], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'marc.ar' & dom_yr$plot == 11], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'kbs.us' & dom_yr$plot == 40], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')
plot(dom_yr$perc_rank[dom_yr$site_code == 'pinj.au' & dom_yr$plot == 19], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt')


# middling ranges
osc[abs(osc$estimate) < .4,] %>% arrange(desc(abs(estimate)))

par(mfrow = c(2,3))
plot(dom_yr$perc_rank[dom_yr$site_code == 'cdcr.us' & dom_yr$plot == 18], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt', main = 'phi = -.04') #-0.04
plot(dom_yr$perc_rank[dom_yr$site_code == 'chilcas.ar' & dom_yr$plot == 45 & dom_yr$Taxon == 'PANICUM BERGII'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt', main = 'phi= -.20') # -0.2
plot(dom_yr$perc_rank[dom_yr$site_code == 'yarra.au' & dom_yr$plot == 2 & dom_yr$Taxon == 'CYNODON DACTYLON'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt', main = 'phi= -.40') # -.4
plot(dom_yr$perc_rank[dom_yr$site_code == 'kiny.au' & dom_yr$plot == 30], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt', main = 'phi = .05') # 0.05 - this looks more oscillation-y
plot(dom_yr$perc_rank[dom_yr$site_code == 'spin.us' & dom_yr$plot == 9], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt', main = 'phi = .20') # 0.2
plot(dom_yr$perc_rank[dom_yr$site_code == 'mtca.au' & dom_yr$plot == 38 & dom_yr$Taxon == 'STIPA NITIDA_TRICHOPHYLLA'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt', main = 'phi = .40') # 0.4


#.4/-.4
plot(dom_yr$perc_rank[dom_yr$site_code == 'yarra.au' & dom_yr$plot == 2 & dom_yr$Taxon == 'CYNODON DACTYLON'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt') # -.4
plot(dom_yr$perc_rank[dom_yr$site_code == 'mtca.au' & dom_yr$plot == 38 & dom_yr$Taxon == 'STIPA NITIDA_TRICHOPHYLLA'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt') # 0.4
plot(dom_yr$perc_rank[dom_yr$site_code == 'shps.us' & dom_yr$plot == 14 & dom_yr$Taxon == 'STIPA COMATA'], ylim = c(0,1), ylab= 'Rank', xlab = 'Year_trt') # -.4



rank_loss_components <- rank_loss_components %>% 
  mutate(id = paste(site_code,plot,Taxon,sep='.')) %>% 
  left_join(.,osc %>% mutate(phi = estimate) %>% dplyr::select(plot, phi), by = c('id'='plot'))



plot(phi ~ inertia_last, data = rank_loss_components)

cor(rank_loss_components[,c('resistance','phi','inertia_first','inertia_last','inertia')], use = 'complete.obs')

rank_loss_components %>% filter(max_year == 15 & abs(phi) < 0.05) %>% print(n=60)
plot(dom_yr$perc_rank[dom_yr$site_code == 'bnch.us' & dom_yr$plot == 26])


#### assign plots that never lose rank as 0 phi (revisit this choice later)

rank_loss_components[is.na(rank_loss_components$phi),]$phi <- 0


#### univariate exploration of rank components ####
#### simple is just main effects
#### decay is based on sig effects from linear model

#### Inertia mod ####

rank_loss_components$rich_scale <- scale_col(rank_loss_components$site_richness)
rank_loss_components$cov_scale <- scale_col(rank_loss_components$initial_rel_cover)
#rank_loss_components$freq_scale <- scale_col(rank_loss_components$plotfreq) # too skewed to get sensible range


pois_inertia_decay <- glmer(
  inertia_first ~ NPK * Fence * cov_scale + local_lifespan * NPK + local_lifespan * Fence + rich_scale * NPK + ppt_scale * Fence  + ppt_var_scale +
    (1|site_code), data = rank_loss_components, family = poisson, control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e6)))

pois_inertia_simple <- glmer(
  inertia_first ~ NPK * Fence + cov_scale + local_lifespan + rich_scale + ppt_scale  + ppt_var_scale +
    (1|site_code), data = rank_loss_components, family = poisson)

R2_iner <- partR2(pois_inertia_simple, data = rank_loss_components,
                  partvars = c("NPK", "Fence","cov_scale", "local_lifespan","plotfreq","rich_scale","ppt_scale","ppt_var_scale"), 
                  R2_type = "marginal", nboot = 10, max_level = 1,
                  olre= FALSE)                               

R2_iner

summary(pois_inertia_simple)
summary(pois_inertia_decay)

r.squaredGLMM(pois_inertia_simple)
# R2m       R2c
# delta     0.2773455 0.7552135
# lognormal 0.2864548 0.7800181
# trigamma  0.2661154 0.7246338

r.squaredGLMM(pois_inertia_decay)


iner_table <- data.frame(summary(pois_inertia_simple)$coef) %>% 
  tibble::rownames_to_column(var = 'Effect') %>%
  dplyr::mutate(Estimate = round(Estimate, digits = 2),
                Std.Error = round(Std..Error, digits = 2),
                z.value = round(z.value, digits = 2), 
                p.value = round(Pr...z.., digits = 3)
  ) %>%
  dplyr::mutate(p.value = ifelse(p.value == 0.000, '< 0.001', p.value)) %>%
  dplyr::select(Effect,Estimate,Std.Error,z.value,p.value) %>% 
  mutate(Effect = gsub('Perennial','',Effect)) %>% 
  left_join(.,R2_iner$R2 %>% 
              mutate(Est_R2 = round(estimate,digits = 2)) %>% 
              mutate(Est_R2 = ifelse(Est_R2 == 0, '< 0.01',Est_R2)) %>% 
              dplyr::select(term,Est_R2),
            by = c('Effect' = 'term')
  )

iner_table

# write_csv(iner_table,
#           paste0('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NutNet manuscripts/Initial dominance/Tables/inertia-response-table_',
#                  Sys.Date(),'.csv'))


iner_trt_emm <- data.frame(summary(emmeans(pois_inertia_simple, pairwise ~ NPK * Fence, type = 'response')$emmeans))
iner_trt_emm <- iner_trt_emm %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

iner_init_emm <- data.frame(summary(emmeans(pois_inertia_simple, ~cov_scale,
                        at = list(cov_scale = seq(from = -1.5, to = 1.5, by = 3)),
                        cov.reduce = TRUE,
                        type = "response")))

iner_rich_emm <- data.frame(summary(emmeans(pois_inertia_simple, ~rich_scale,
                                            at = list(rich_scale = seq(from = -1.5, to = 1.5, by = 3)),
                                            cov.reduce = TRUE,
                                            type = "response")))

iner_freq_emm <- data.frame(summary(emmeans(pois_inertia_simple, ~plotfreq,
                                            at = list(plotfreq = seq(from = 0.5, to = 1, by = 0.5)),
                                            cov.reduce = TRUE,
                                            type = "response")))

iner_ppt_emm <- data.frame(summary(emmeans(pois_inertia_simple, ~ppt_scale,
                                            at = list(ppt_scale = seq(from = -1.5, to = 1.5, by = 3)),
                                            cov.reduce = TRUE,
                                            type = "response")))

iner_pptvar_emm <- data.frame(summary(emmeans(pois_inertia_simple, ~ppt_var_scale,
                                           at = list(ppt_var_scale = seq(from = -1, to = 1.5, by = 2.5)),
                                           cov.reduce = TRUE,
                                           type = "response")))


iner_life_emm <- data.frame(summary(emmeans(pois_inertia_simple, ~ local_lifespan,
                        type = "response")))

iner_trt_emm
iner_life_emm
iner_init_emm$initial_cover <- factor(ifelse(iner_init_emm$cov_scale == -1.5, 'Low Initial Cover','High Initial Cover'), levels = c('Low Initial Cover','High Initial Cover'))
iner_rich_emm$site_rich <- factor(ifelse(iner_rich_emm$rich_scale == -1.5, 'Low Site Richness','High Site Richness'), levels = c('Low Site Richness','High Site Richness'))
iner_freq_emm$frequency <- factor(ifelse(iner_freq_emm$plotfreq == 0.5, 'Low Initial Frequency','High Initial Frequency'), levels = c('Low Initial Frequency','High Initial Frequency'))
iner_ppt_emm$ppt <- factor(ifelse(iner_ppt_emm$ppt_scale == -1.5, 'Low MAP','High MAP'), levels = c('Low MAP','High MAP'))
iner_pptvar_emm$pptvar <- factor(ifelse(iner_pptvar_emm$ppt_var_scale == -1, 'Low Precipitation Variation','High Precipitation Variation'), levels = c('Low Precipitation Variation','High Precipitation Variation'))


iner_effects <- bind_rows(iner_trt_emm %>% 
                            rename(emmean = rate, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,trt) %>% 
                            rename(predictor = trt),
                          iner_init_emm %>% 
                            rename(emmean = rate, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,initial_cover) %>% 
                            rename(predictor = initial_cover),
                          iner_freq_emm %>% 
                            rename(emmean = rate, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,frequency) %>% 
                            rename(predictor = frequency),
                          iner_life_emm %>% 
                            rename(emmean = rate, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,local_lifespan) %>% 
                            rename(predictor = local_lifespan),
                          iner_rich_emm %>% 
                            rename(emmean = rate, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,site_rich) %>% 
                            rename(predictor = site_rich),
                          iner_ppt_emm %>% 
                            rename(emmean = rate, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,ppt) %>% 
                            rename(predictor = ppt),
                          iner_pptvar_emm %>% 
                            rename(emmean = rate, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,pptvar) %>% 
                            rename(predictor = pptvar)
                          )




iner_effects$predictor <-factor(iner_effects$predictor, 
                                levels = c('High Precipitation Variation','Low Precipitation Variation',
                                           'High MAP','Low MAP',
                                           'High Site Richness','Low Site Richness',
                                           'High Initial Frequency','Low Initial Frequency',
                                           'High Initial Cover','Low Initial Cover',
                                           'Annual','Perennial',
                                           'NPK+Fence','NPK','Fence','Control'))
  
iner_sig <- data.frame(
  x = c(14.5,11.5,9.5,7.5),
  y = c(3.5,3.5),
  label = c("*", "*"))

gg_iner <- ggplot(iner_effects,aes(y = emmean,x = predictor)) +
  geom_point(col = 'black',size = 4) +
  geom_errorbar(aes(ymax= upper.CL, ymin = lower.CL),width=0) +
  coord_flip() + geom_vline(xintercept = c(12.5,10.5,8.5,6.5,4.5,2.5),size = 2) +
  xlab('') + ylab('Inertia\nEstimated years to loss of rank') +
  scale_y_continuous(breaks = c(log(2),log(4),log(8),log(16)),labels=c(2,4,8,16))+
  #scale_y_continuous(breaks = log(pretty(exp(iner_effects$emmean))), labels = pretty(exp(iner_effects$emmean))) +
  geom_text(data=iner_sig, aes( x=x, y=y, label=label),
            color="darkred", 
            size=12 , angle=45, fontface="bold" )

gg_iner

### resistance main effects ####

logit_resistance_simple <- glmer(
  resistance ~ NPK * Fence + cov_scale + local_lifespan + rich_scale + plotfreq + ppt_scale + ppt_var_scale +
    (1|site_code/Taxon), data = rank_loss_components,
  family = binomial, control = glmerControl(optimizer = "bobyqa",
                                            optCtrl=list(maxfun=2e6)))


logit_resistance_decay <- glmer(resistance ~ 
                            NPK * Fence * cov_scale + local_lifespan * NPK + local_lifespan * Fence + plotfreq + rich_scale * NPK + ppt_scale * Fence  + ppt_var_scale +
                                 (1 | site_code), data = rank_loss_components,
                               family = binomial, control = glmerControl(optimizer = "bobyqa",
                                                                         optCtrl=list(maxfun=2e6)))

# mod_resistance <- lmer(resistance ~ (NPK + Fence + initial_rel_cover)^3 + (NPK + Fence + local_lifespan)^3 +
#                       (1 | site_code), data = rank_loss_components)

anova(logit_resistance_simple)
summary(logit_resistance_simple)

anova(logit_resistance_decay)
summary(logit_resistance_decay)

r.squaredGLMM(logit_resistance_simple)

R2_res <- partR2(logit_resistance_simple, data = rank_loss_components,
                  partvars = c("NPK", "Fence","cov_scale", "local_lifespan","plotfreq","rich_scale","ppt_scale","ppt_var_scale"), 
                  R2_type = "marginal", nboot = 10, max_level = 1,
                  olre= FALSE)                               

R2_res$R2

res_table <- data.frame(summary(logit_resistance_simple)$coef) %>% 
  tibble::rownames_to_column(var = 'Effect') %>%
  dplyr::mutate(Estimate = round(Estimate, digits = 2),
                Std.Error = round(Std..Error, digits = 2),
                z.value = round(z.value, digits = 2), 
                p.value = round(Pr...z.., digits = 3)
                ) %>%
  dplyr::mutate(p.value = ifelse(p.value == 0.000, '< 0.001', p.value)) %>%
  dplyr::select(Effect,Estimate,Std.Error,z.value,p.value) %>% 
  mutate(Effect = gsub('Perennial','',Effect)) %>% 
  left_join(.,R2_res$R2 %>% 
              mutate(Est_R2 = round(estimate,digits = 2)) %>% 
              mutate(Est_R2 = ifelse(Est_R2 == 0, '< 0.01',Est_R2)) %>% 
              dplyr::select(term,Est_R2),
            by = c('Effect' = 'term')
            )

res_table

# write_csv(res_table,
#           paste0('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NutNet manuscripts/Initial dominance/Tables/resistance-response-table_',
#                  Sys.Date(),'.csv'))


res_trt_emm <- data.frame(summary(emmeans(logit_resistance_simple, pairwise ~ NPK * Fence, type = 'response')$emmeans))
res_trt_emm <- res_trt_emm %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

res_init_emm <- data.frame(summary(emmeans(logit_resistance_simple, ~cov_scale,
                                            at = list(cov_scale = seq(from = -1.5, to = 1.5, by = 3)),
                                            cov.reduce = TRUE,
                                            type = "response")))

res_rich_emm <- data.frame(summary(emmeans(logit_resistance_simple, ~rich_scale,
                                            at = list(rich_scale = seq(from = -1.5, to = 1.5, by = 3)),
                                            cov.reduce = TRUE,
                                            type = "response")))

res_freq_emm <- data.frame(summary(emmeans(logit_resistance_simple, ~plotfreq,
                                            at = list(plotfreq = seq(from = 0.5, to = 1, by = 0.5)),
                                            cov.reduce = TRUE,
                                            type = "response")))

res_ppt_emm <- data.frame(summary(emmeans(logit_resistance_simple, ~ppt_scale,
                                           at = list(ppt_scale = seq(from = -1.5, to = 1.5, by = 3)),
                                           cov.reduce = TRUE,
                                           type = "response")))

res_pptvar_emm <- data.frame(summary(emmeans(logit_resistance_simple, ~ppt_var_scale,
                                              at = list(ppt_var_scale = seq(from = -1, to = 1.5, by = 2.5)),
                                              cov.reduce = TRUE,
                                              type = "response")))


res_life_emm <- data.frame(summary(emmeans(logit_resistance_simple, ~ local_lifespan,
                                            type = "response")))

res_trt_emm
res_life_emm
res_init_emm$initial_cover <- factor(ifelse(res_init_emm$cov_scale == -1.5, 'Low Initial Cover','High Initial Cover'), levels = c('Low Initial Cover','High Initial Cover'))
res_rich_emm$site_rich <- factor(ifelse(res_rich_emm$rich_scale == -1.5, 'Low Site Richness','High Site Richness'), levels = c('Low Site Richness','High Site Richness'))
res_freq_emm$frequency <- factor(ifelse(res_freq_emm$plotfreq == 0.5, 'Low Initial Frequency','High Initial Frequency'), levels = c('Low Initial Frequency','High Initial Frequency'))
res_ppt_emm$ppt <- factor(ifelse(res_ppt_emm$ppt_scale == -1.5, 'Low MAP','High MAP'), levels = c('Low MAP','High MAP'))
res_pptvar_emm$pptvar <- factor(ifelse(res_pptvar_emm$ppt_var_scale == -1, 'Low Precipitation Variation','High Precipitation Variation'), levels = c('Low Precipitation Variation','High Precipitation Variation'))


res_effects <- bind_rows(res_trt_emm %>% 
                            rename(emmean = prob, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,trt) %>% 
                            rename(predictor = trt),
                          res_init_emm %>% 
                            rename(emmean = prob, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,initial_cover) %>% 
                            rename(predictor = initial_cover),
                          res_freq_emm %>% 
                            rename(emmean = prob, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,frequency) %>% 
                            rename(predictor = frequency),
                          res_life_emm %>% 
                            rename(emmean = prob, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,local_lifespan) %>% 
                            rename(predictor = local_lifespan),
                          res_rich_emm %>% 
                            rename(emmean = prob, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,site_rich) %>% 
                            rename(predictor = site_rich),
                          res_ppt_emm %>% 
                            rename(emmean = prob, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,ppt) %>% 
                            rename(predictor = ppt),
                          res_pptvar_emm %>% 
                            rename(emmean = prob, lower.CL = asymp.LCL, upper.CL = asymp.UCL) %>% 
                            dplyr::select(emmean,lower.CL,upper.CL,pptvar) %>% 
                            rename(predictor = pptvar)
)




res_effects$predictor <-factor(res_effects$predictor, 
                                levels = c('High Precipitation Variation','Low Precipitation Variation',
                                           'High MAP','Low MAP',
                                           'High Site Richness','Low Site Richness',
                                           'High Initial Frequency','Low Initial Frequency',
                                           'High Initial Cover','Low Initial Cover',
                                           'Annual','Perennial',
                                           'NPK+Fence','NPK','Fence','Control'))

res_sig <- data.frame(
  x = c(14.5,11.5, 9.5),
  y = c(0.8,0.8,0.8),
  label = c("*", "*", "*"))

summary(logit_resistance_simple)

gg_res <- ggplot(res_effects,aes(y = emmean,x = predictor)) +
  geom_point(col = 'black',size = 4) +
  geom_errorbar(aes(ymax= upper.CL, ymin = lower.CL),width=0) +
  coord_flip() + geom_vline(xintercept = c(12.5,10.5,8.5,6.5,4.5,2.5),size = 2) + ylim(0,1) +
  xlab('') + ylab('Resistance\nEstimated minimum rank') + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_text(data=res_sig, aes( x=x, y=y, label=label),                 , 
            color="darkred", 
            size=12 , angle=45, fontface="bold" )

gg_res






#### Fig 3 ####

gg_iner_res <- plot_grid(gg_iner,gg_res,ncol = 2,rel_widths = c(0.65,0.35))

gg_iner_res

# ggsave(paste0('Graphs/Fig-3-resistance-inertia-response_',Sys.Date(),'.png'),
#        gg_iner_res, width = 15, height = 8, units = "in", dpi = 600)


### oscillation main effects ####

linear_osc_simple <- lmer(
  phi ~ NPK * Fence + cov_scale + local_lifespan + rich_scale + ppt_scale + ppt_var_scale +
    (1|site_code/Taxon), data = rank_loss_components)

linear_osc_decay <- lmer(
  phi ~ NPK * Fence * cov_scale + NPK * Fence * local_lifespan + NPK * Fence *rich_scale + ppt_scale + ppt_var_scale +
    (1|site_code/Taxon), data = rank_loss_components)

summary(linear_osc_simple)
anova(linear_osc_simple)

r.squaredGLMM(linear_osc_simple)
r.squaredGLMM(linear_osc_decay)


#### categorize species and look at changes ####

rank_loss_cats <- rank_loss_components %>% 
  mutate(traj_cat = if_else(
    
  ))

#### predicting changes in community responses (biomass, rich, RAC) ####

dom_comm

#following Seabloom 2020 for extracting rates of change in biomass and diversity
#log tran

dom_comm_10 <- dom_comm %>% 
  mutate(mass = if_else(is.na(vascular_live_mass),unsorted_live_mass,vascular_live_mass)) %>% 
  mutate(mass.lg = log10(mass), rich.lg = log10(subord_rich), rac.lg = log10(rank_change), yr.lg = log10(year_trt+1)) %>% 
  filter(site_code %in% dom_pop_yr$site_code, year_trt <= 10)
  #mutate(LRR_mass = if_else(is.na(live_mass),log(unsorted_mass/ctrl_mass),log(live_mass/ctrl_mass)))
       
  dom_comm_10 <- do.call(data.frame, lapply(dom_comm_10, function(x) replace(x, is.infinite(x), NA)))
         
# npp.simp <- lme(mass.lg ~ trt*yr.lg, random = ~1|site_code, data=dom_comm_10, na.action=na.omit)
# summary(npp.simp)
# #lm.npp2t.big <- lme(mass.lg ~ trt*yr.lg, random = ~1 + yr.lg+trt|site_code/plot, data=dom_comm_10, na.action=na.omit)
# # does not converge
# #summary(lm.npp2t.big)
# 
# #slightly simpler
# lm.npp2t.big <- lme(mass.lg ~ trt*yr.lg, random = ~1 + yr.lg+trt|site_code, data=dom_comm_10, na.action=na.exclude)
# summary(lm.npp2t.big)
# # converges!
# anova(npp.simp,lm.npp2t.big)
# 
# rich.simp <- lme(rich.lg ~ trt*yr.lg, random = ~1|site_code, data=dom_comm_10, na.action=na.omit)
# summary(rich.simp)
# lm.rich.big <- lme(rich.lg ~ trt*yr.lg, random = ~1 + yr.lg+trt|site_code, data=dom_comm_10, na.action=na.exclude)
# summary(lm.rich.big)
# anova(rich.simp, lm.rich.big)
# 
# rac.simp <- lme(rac.lg ~ trt*yr.lg, random = ~1|site_code, data=dom_comm_10, na.action=na.omit)
# summary(rac.simp)
# lm.rac.big <- lme(rac.lg ~ trt*yr.lg, random = ~1 + yr.lg+trt|site_code, data=dom_comm_10, na.action=na.exclude)
# summary(lm.rac.big)
# anova(rac.simp, lm.rac.big)
# 
# # Get plot level predictions and take site level means
# dom_comm_10$mass.lg.pred <- predict(lm.npp2t.big)  # live mass
# dom_comm_10$rich.lg.pred <- predict(lm.rich.big) # subordinate richness
# dom_comm_10$rac.lg.pred <- predict(lm.rac.big)   # subordinate rank change
# 
# # Backtransform to raw data 
# dom_comm_10$mass.pred <- 10^dom_comm_10$mass.lg.pred
# dom_comm_10$rich.pred <- 10^dom_comm_10$rich.lg.pred
# dom_comm_10$rac.pred <- 10^dom_comm_10$rac.lg.pred
# 
# 
# ggplot(dom_comm_10[dom_comm_10$site_code == 'cdcr.us',], aes(x=year_trt,y=mass.pred, col=factor(trt))) + geom_point()
# ggplot(dom_comm_10[dom_comm_10$site_code == 'cdcr.us',], aes(x=yr.lg,y=mass.lg.pred, col=factor(trt))) + geom_line()
# 
# coef(lm.npp2t.big)


#### look at whether rates of dominance loss correlate to rates of community changes in plots
  
  plot_rank_loss <- split(dom_pop_yr, ~ site_code + plot) %>% 
    discard(~nrow(.)==0) %>%
    map(~lm(perc_rank-1 ~ year_trt-1, data = .,na.action=na.omit))  %>%  #one way to fix intercept at '1'
    map_df(tidy,.id='plot') %>%
    rename(est_rank = estimate) 
  
  range(plot_rank_loss$est_rank)
  hist(plot_rank_loss$est_rank,nclass= 100)
  
filter(plot_rank_loss, est_rank == 0)

plot_mass_loss <- filter(dom_comm_10,!is.na(mass)) %>%
  split(dom_comm_10, ~ site_code + plot) %>%
  discard(~nrow(.)==0) %>%
  map(~lm(mass ~ year_trt-1, data = .,na.action=na.omit))  %>%  #one way to fix intercept at '1'
  map_df(tidy,.id='plot') %>%
  rename(est_mass = estimate)

plot_rich_loss <- split(dom_comm_10, ~ site_code + plot) %>% 
  discard(~nrow(.)==0) %>%
  map(~lm(subord_rich ~ year_trt, data = .,na.action=na.omit))  %>% 
  map_df(tidy,.id='plot') %>%
  filter(term != '(Intercept)') %>% 
  rename(est_rich = estimate) 

rich_rank <- comm_rank_change <- rank_loss_components %>% 
  mutate(id = paste0(site_code,'.',plot)) %>% 
  left_join(.,plot_rank_loss %>% 
              dplyr::select(plot,est_rank),by=c('id'='plot')) %>% 
  left_join(.,plot_rich_loss %>% 
              dplyr::select(plot,est_rich),by=c('id'='plot'))

rich_rank_mod <- lme(est_rich ~ trt * est_rank,
                    random = (~1|site_code),
                    data = rich_rank,
                    na.action = na.omit)

summary(rich_rank_mod)


gg_rich_rank <- ggplot(rich_rank, aes(x= est_rank, y = est_rich, col = trt)) +
  geom_point(alpha = 0.4) + geom_smooth(method = 'lm', size = 2) +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) + 
  xlab('Estimated rank decay slope') + ylab('Estimated change in\n subordinate richness with time')
gg_rich_rank

rich_rank %>% filter(est_rank < -.1 & est_rich < -1 & trt == 'NPK') %>% dplyr::select(est_rank) # marc.ar plot 18

marc18.rank <- filter(dom_pop_yr, site_code == 'marc.ar' & plot == 18)
gg_marc_rank <- ggplot(marc18.rank, aes(x = year_trt, y = perc_rank)) + geom_point(col= 'black',fill='black', size = 3) + 
  geom_smooth(method = 'lm', se = FALSE) +
  xlab('Year after treatment') + ylab('Rank Percentile') + annotate("text",x = 8, y = 0.8,label = "Mar Chiquita - plot 18\n Bromus catharticus")
gg_marc_rank

marc18.rich <- filter(dom_comm_10, site_code == 'marc.ar' & plot == 18)
gg_marc_rich <- ggplot(marc18.rich, aes(x = year_trt, y = subord_rich)) + geom_point(col= 'black',fill='black', size = 3) + geom_smooth(method = 'lm', se = FALSE) +
  xlab('Year after treatment') + ylab('Subordinate richness')
gg_marc_rich

  
dom_comm_LRR <- dom_comm_10 %>% 
  group_by(site_code,year_trt) %>% 
    mutate(ctrl_mass = if_else(is.na(vascular_live_mass),mean(unsorted_live_mass[trt == 'Control'],na.rm=T),mean(vascular_live_mass[trt == 'Control'],na.rm=T)),
           ctrl_RAC = mean(rank_change[trt == 'Control'],na.rm=T),
           ctrl_rich = mean(subord_rich[trt == 'Control'],na.rm=T)) %>% 
    mutate(LRR_mass = if_else(is.na(vascular_live_mass),log10(unsorted_live_mass/ctrl_mass),log10(vascular_live_mass/ctrl_mass)),
           LRR_RAC = log10(rank_change/ctrl_RAC),
           LRR_rich = log10(subord_rich/ctrl_rich)) 

# Replace Inf in data by NA
dom_comm_LRR <- do.call(data.frame,                      
                   lapply(dom_comm_LRR,
                          function(x) replace(x, is.infinite(x), NA)))

LRR_mass <- filter(dom_comm_LRR,site_code %in% dom_pop_yr$site_code,year_trt <= yr_cut,!is.na(LRR_mass),!is.infinite(LRR_mass)) %>% 
    split(., ~ site_code + plot) %>% 
    discard(~nrow(.)==0) %>%
    map(~lm(LRR_mass ~ yr.lg - 1, data = .,na.action=na.omit))  %>% 
    map_df(tidy,.id='plot') %>%
    rename(est_bio = estimate)

LRR_rich <- filter(dom_comm_LRR,site_code %in% dom_pop_yr$site_code,year_trt <= yr_cut,!is.na(LRR_rich),!is.infinite(LRR_rich)) %>% 
  split(., ~ site_code + plot) %>% 
  discard(~nrow(.)==0) %>%
  map(~lm(LRR_rich ~ yr.lg - 1, data = .,na.action=na.omit))  %>% 
  map_df(tidy,.id='plot') %>%
  rename(est_rich = estimate)

LRR_rac <- filter(dom_comm_LRR,site_code %in% dom_pop_yr$site_code,year_trt <= yr_cut,!is.na(LRR_RAC),!is.infinite(LRR_RAC)) %>% 
  split(., ~ site_code + plot) %>% 
  discard(~nrow(.)==0) %>%
  map(~lm(LRR_RAC ~ yr.lg - 1, data = .,na.action=na.omit))  %>% 
  map_df(tidy,.id='plot') %>%
  rename(est_rac = estimate)

comm_rank_change <- rank_loss_components %>% 
  mutate(id = paste0(site_code,'.',plot)) %>% 
  left_join(.,LRR_mass %>% 
              dplyr::select(plot,est_bio),by=c('id'='plot')) %>% 
  left_join(.,LRR_rich %>% 
              dplyr::select(plot,est_rich),by=c('id'='plot')) %>% 
  left_join(.,LRR_rac %>% 
              dplyr::select(plot,est_rac),by=c('id'='plot')) %>% 
  left_join(.,dom_comm_LRR %>% 
              mutate(id = paste0(site_code,'.',plot)) %>% 
              group_by(id) %>% 
              summarize(max_LRR_mass = max(LRR_mass[is.finite(LRR_mass)],na.rm=TRUE),
                        min_LRR_rich = min(LRR_rich[is.finite(LRR_rich)],na.rm=TRUE),
                        max_LRR_rac = max(LRR_RAC[is.finite(LRR_RAC)],na.rm=TRUE)),
            by = 'id') %>% 
  mutate(yr.lg = log10(year_trt + 1))



LRR_components <- rank_loss_components %>% 
  dplyr::select(-year_trt) %>% 
  left_join(.,dom_comm_LRR %>% 
              dplyr::select(year_trt,yr.lg,site_code,plot,Taxon,richness_change,rank_change,losses,LRR_mass,LRR_RAC,LRR_rich),
            by = c('site_code','plot','Taxon')) %>% 
  filter(trt != 'Control')

LRR_bio_components_mod <- lme(LRR_mass ~ trt * yr.lg* inertia_last +
                              trt * yr.lg * resistance +
                              trt * yr.lg * initial_rel_cover - 1,
                        random = (~1|site_code/plot),
                        data = LRR_components,
                        na.action = na.omit)

summary(LRR_bio_components_mod)
anova(LRR_bio_components_mod)
r.squaredGLMM(LRR_bio_components_mod)


LRR_rich_components_mod <- lme(LRR_rich ~ trt * year_trt * inertia_last +
                                trt * year_trt * resistance+
                                 trt * year_trt * initial_rel_cover - 1,
                              random = (~1|site_code/plot),
                              data = LRR_components,
                              na.action = na.omit)

summary(LRR_rich_components_mod)
r.squaredGLMM(LRR_rich_components_mod)


rac_components_mod <- lme(rank_change ~ trt * year_trt * inertia_last +
                                 trt * year_trt * resistance +
                                trt * year_trt * initial_rel_cover - 1,
                               random = (~1|site_code/plot),
                               data = LRR_components,
                               na.action = na.omit)

summary(rac_components_mod)
r.squaredGLMM(rac_components_mod)


summary.tablefunc(LRR_bio_components_mod)
summary.tablefunc(LRR_rich_components_mod)
summary.tablefunc(rac_components_mod)

# write_csv(summary.tablefunc(LRR_bio_components_mod),
#           paste0('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NutNet manuscripts/Initial dominance/Tables/LRR-bio-response_',
#                  Sys.Date(),'.csv'))
# 
# write_csv(summary.tablefunc(LRR_rich_components_mod),
#           paste0('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NutNet manuscripts/Initial dominance/Tables/LRR-rich-response_',
#                  Sys.Date(),'.csv'))
# 
# write_csv(summary.tablefunc(rac_components_mod),
#           paste0('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NutNet manuscripts/Initial dominance/Tables/LRR-RAC-response_',
#                  Sys.Date(),'.csv'))


iner_rich_emm <- data.frame(summary(emmeans(LRR_rich_components_mod, pairwise ~ trt | inertia_last * year_trt,
                                             at = list(year_trt = seq(from = 1, to = 13, by = 1),
                                                       inertia_last = seq(from = 0, to = 10, by = 10)),
                                             cov.reduce = TRUE,
                                             type = "response")$emmeans))

iner_rich_emm$inertia <- factor(ifelse(iner_rich_emm$inertia_last == 0, 'Low Inertia','High Inertia'), levels = c('Low Inertia','High Inertia'))
# rank_init_emm$upper.CL <- ifelse(rank_init_emm$upper.CL > 1.05, 1.05, rank_init_emm$upper.CL)
# rank_init_emm$lower.CL <- ifelse(rank_init_emm$lower.CL < 0, 0, rank_init_emm$lower.CL)


gg_iner_richLRR <- ggplot(iner_rich_emm, aes(x = year_trt, y = emmean, col = trt)) +
  #geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.1, alpha = 0.8, col = 'grey90') +
  geom_ribbon(aes(ymin = lower.CL,ymax= upper.CL,fill=trt),alpha=0.35,color=NA) +
  geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt), size = 3) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('LRR Richness\n(Initial subordinates)') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal[2:4]) +
  scale_fill_manual(values = ctrl_pal[2:4]) +
  ylim(-1,0.5) +
  geom_hline(yintercept=0, lty = 'dashed') +
  facet_wrap(~inertia)

gg_iner_richLRR

# ggsave(paste0('Graphs/Fig-4a-iner-rich_',Sys.Date(),'.png'),
#        gg_cov_rac, width = 12, height = 8, units = "in", dpi = 600)



##### rac response



cov_rac_emm <- data.frame(summary(emmeans(rac_components_mod, pairwise ~ trt | initial_rel_cover * year_trt,
                                           at = list(year_trt = seq(from = 1, to = 13, by = 1),
                                                     initial_rel_cover = seq(from = 0.2, to = 1, by = 0.8)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))

cov_rac_emm$cover <- factor(ifelse(cov_rac_emm$initial_rel_cover == 0.2, 'Low Initial Cover','High Initial Cover'), levels = c('Low Initial Cover','High Initial Cover'))
# rank_init_emm$upper.CL <- ifelse(rank_init_emm$upper.CL > 1.05, 1.05, rank_init_emm$upper.CL)
# rank_init_emm$lower.CL <- ifelse(rank_init_emm$lower.CL < 0, 0, rank_init_emm$lower.CL)


gg_cov_rac <- ggplot(cov_rac_emm, aes(x = year_trt, y = emmean, col = trt)) +
  #geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.1, alpha = 0.8, col = 'grey90') +
  geom_ribbon(aes(ymin = lower.CL,ymax= upper.CL,fill=trt),alpha=0.35,color=NA) +
  geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt), size = 3) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank Abundance Change\n(Initial subordinates)') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal[2:4]) +
  scale_fill_manual(values = ctrl_pal[2:4]) +
  ylim(0,0.5) +
  #geom_hline(yintercept=0, lty = 'dashed') +
  facet_wrap(~cover)

gg_cov_rac


# ggsave(paste0('Graphs/Fig-4b-rac-cover_',Sys.Date(),'.png'),
#        gg_cov_rac, width = 12, height = 8, units = "in", dpi = 600)


### Fig 4 ####



LRR_legend <- get_legend(gg_cov_rac + theme(legend.direction = "horizontal"))

gg_iner_fig <- gg_iner_richLRR + theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank(),
                               plot.tag.position = c(0.08,1)) + 
  labs(tag = 'A')

gg_cov_fig <- gg_cov_rac + theme(legend.position="none",
                             plot.tag.position = c(0.08,1)) + 
  labs(tag = 'B')

lay <- rbind(c(1,1,1),
             c(2,2,2),
             c(3,4,3))

Fig4.LRR <- grid.arrange(gg_iner_fig,
                         gg_cov_fig,
                         blankPlot,
                         LRR_legend,
                         layout_matrix = lay,
                          heights = c(2.5,3,0.5))
                          #left = textGrob("Rank percentage", rot = 90, gp=gpar(fontsize=24)),
                          #bottom = textGrob("Year after treatment", gp=gpar(fontsize=24)))

# save_plot(paste0('Graphs/Fig4-System-response_',Sys.Date(),'.png'),
#           Fig4.LRR,
#           base_height = 10, base_width = 10)



######## deprecated 

res_rich_emm <- data.frame(summary(emmeans(LRR_rich_components_mod, pairwise ~ trt | resistance * year_trt,
                                           at = list(year_trt = seq(from = 1, to = 13, by = 1),
                                                     resistance = seq(from = 0, to = 1, by = 1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))

res_rich_emm$resist <- factor(ifelse(res_rich_emm$resistance == 0, 'Low Resistance','High Resistance'), levels = c('Low Resistance','High Resistance'))
# rank_init_emm$upper.CL <- ifelse(rank_init_emm$upper.CL > 1.05, 1.05, rank_init_emm$upper.CL)
# rank_init_emm$lower.CL <- ifelse(rank_init_emm$lower.CL < 0, 0, rank_init_emm$lower.CL)


gg_res_richLRR <- ggplot(res_rich_emm, aes(x = year_trt, y = emmean, col = trt)) +
  #geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.1, alpha = 0.8, col = 'grey90') +
  geom_ribbon(aes(ymin = lower.CL,ymax= upper.CL,fill=trt),alpha=0.35,color=NA) +
  geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt), size = 3) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('LRR Richness') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal[2:4]) +
  scale_fill_manual(values = ctrl_pal[2:4]) +
  ylim(-1,0.5) +
  geom_hline(yintercept=0, lty = 'dashed') +
  facet_wrap(~resist)

gg_res_richLRR

cov_rich_emm <- data.frame(summary(emmeans(LRR_rich_components_mod, pairwise ~ trt | initial_rel_cover * year_trt,
                                           at = list(year_trt = seq(from = 1, to = 13, by = 1),
                                                     initial_rel_cover = seq(from = 0.2, to = 1, by = 0.8)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))

cov_rich_emm$cover <- factor(ifelse(cov_rich_emm$initial_rel_cover == 0.2, 'Low Initial Cover','High Initial Cover'), levels = c('Low Initial Cover','High Initial Cover'))
# rank_init_emm$upper.CL <- ifelse(rank_init_emm$upper.CL > 1.05, 1.05, rank_init_emm$upper.CL)
# rank_init_emm$lower.CL <- ifelse(rank_init_emm$lower.CL < 0, 0, rank_init_emm$lower.CL)


gg_cov_richLRR <- ggplot(cov_rich_emm, aes(x = year_trt, y = emmean, col = trt)) +
  #geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.1, alpha = 0.8, col = 'grey90') +
  geom_ribbon(aes(ymin = lower.CL,ymax= upper.CL,fill=trt),alpha=0.35,color=NA) +
  geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt), size = 3) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('LRR Richness') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal[2:4]) +
  scale_fill_manual(values = ctrl_pal[2:4]) +
  ylim(-1,0.5) +
  geom_hline(yintercept=0, lty = 'dashed') +
  facet_wrap(~cover)

gg_cov_richLRR



iner_rac_emm <- data.frame(summary(emmeans(rac_components_mod, pairwise ~ trt | inertia_last * year_trt,
                                           at = list(year_trt = seq(from = 1, to = 13, by = 1),
                                                     inertia_last = seq(from = 0, to = 10, by = 10)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))

iner_rac_emm$inertia <- factor(ifelse(iner_rac_emm$inertia_last == 0, 'Low Inertia','High Inertia'), levels = c('Low Inertia','High Inertia'))
# rank_init_emm$upper.CL <- ifelse(rank_init_emm$upper.CL > 1.05, 1.05, rank_init_emm$upper.CL)
# rank_init_emm$lower.CL <- ifelse(rank_init_emm$lower.CL < 0, 0, rank_init_emm$lower.CL)


gg_iner_racLRR <- ggplot(iner_rac_emm, aes(x = year_trt, y = emmean, col = trt)) +
  #geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.1, alpha = 0.8, col = 'grey90') +
  geom_ribbon(aes(ymin = lower.CL,ymax= upper.CL,fill=trt),alpha=0.35,color=NA) +
  geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt), size = 3) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank Abundance Change') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal[2:4]) +
  scale_fill_manual(values = ctrl_pal[2:4]) +
  ylim(0,0.5) +
  #geom_hline(yintercept=0, lty = 'dashed') +
  facet_wrap(~inertia)

gg_iner_racLRR



res_rac_emm <- data.frame(summary(emmeans(rac_components_mod, pairwise ~ trt | resistance * year_trt,
                                          at = list(year_trt = seq(from = 1, to = 13, by = 1),
                                                    resistance = seq(from = 0, to = 1, by = 1)),
                                          cov.reduce = TRUE,
                                          type = "response")$emmeans))

res_rac_emm$resist <- factor(ifelse(res_rac_emm$resistance == 0, 'Low Resistance','High Resistance'), levels = c('Low Resistance','High Resistance'))
# rank_init_emm$upper.CL <- ifelse(rank_init_emm$upper.CL > 1.05, 1.05, rank_init_emm$upper.CL)
# rank_init_emm$lower.CL <- ifelse(rank_init_emm$lower.CL < 0, 0, rank_init_emm$lower.CL)


gg_res_rac <- ggplot(res_rac_emm, aes(x = year_trt, y = emmean, col = trt)) +
  #geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.1, alpha = 0.8, col = 'grey90') +
  geom_ribbon(aes(ymin = lower.CL,ymax= upper.CL,fill=trt),alpha=0.35,color=NA) +
  geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt), size = 3) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('LRR racness') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal[2:4]) +
  scale_fill_manual(values = ctrl_pal[2:4]) +
  ylim(0,0.5) +
  #geom_hline(yintercept=0, lty = 'dashed') +
  facet_wrap(~resist)

gg_res_rac


