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

## graphics
library(ggplot2)
library(cowplot)
library(GGally)


theme_figs <- theme_classic() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
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

## Read in data
setwd("~/Dropbox/NutNet data")
files<-list.files()
# coverfile<-files[grep("full-cover-", files)]
# coverfile<-last(coverfile)
# original_coverdat <- read_csv(coverfile)

original_coverdat <- read.csv('/Users/wilf0020/Dropbox/Peter/Cover and biomass by date/full-cover-by-date-2022-10-10.csv')
coverdat <- original_coverdat # backup

dim(coverdat)
# 305412   x  21 - sorted by date

comb_file<-files[grep("comb-by-plot-clim-soil-diversity", files)]
comb_file<-last(comb_file)
original_comb <- read_csv(comb_file)
comb <- original_comb # backup

cover_fix <- read_csv('https://github.com/jon-bakker/NutNet_TaxonomicAdjustments/raw/main/taxonomic-adjustments-2022-05-02.csv')

# choose all experimental sites, with at least an N treatment, 1 year post-treatment:
sites <-unique(comb$site_code[comb$trt %in% c("NPK", "Fence") # at least an NPK treatment or fence
                               & (!comb$experiment_type %in% c("Observational")) # not observational site
                               &  comb$year_trt >=1 # has at least one year post-treatment data
                               ]) 


length(sites)
sort(sites)


# 95
# just the sites we want
working_coverdat <- coverdat[coverdat$site_code %in% sites,]

### lagoas has off-season pre-treatment year, but data looks sufficient to use after checking

working_coverdat$year_trt <- ifelse(working_coverdat$site_code == 'lagoas.br' & working_coverdat$year_trt == -1,
                                    0, working_coverdat$year_trt)


## need to make sure they have pre-treatment data as well
site.miss <- sites[!sites %in%  unique(working_coverdat$site_code[working_coverdat$year_trt == 0])]

sites <- sort(unique(working_coverdat$site_code[working_coverdat$year_trt == 0]))

working_coverdat <- coverdat[coverdat$site_code %in% sites,]


# write.csv(sites,
#           '/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NutNet manuscripts/Initial dominance/Data/site list.csv')



working_coverdat <- droplevels(working_coverdat) %>% 
  filter(live == 1, !functional_group %in% c('BRYOPHYTE','LICHEN','LIVERWORT'), max_cover != 0) %>% 
  group_by(site_code) %>% 
  mutate(max_year = max(year_trt))  ##  add in max years

working_coverdat[working_coverdat$Taxon == 'RUMEX SP.',]$local_lifespan <- 'PERENNIAL'

working_coverdat$functional_group <- ifelse(working_coverdat$functional_group %in% c('GRASS','GRAMINOID'), 'GRAMINOID',working_coverdat$functional_group)

### doane.us only sampled block 1 after pre-treatment year

working_coverdat <- working_coverdat %>% 
  filter(!(site_code == 'doane.us' & block != 1))

#### Jena was only sampled in september in pre-treatment, so only using sept data thereafter

unique(working_coverdat[working_coverdat$site_code == 'jena.de',]$date)
dim(working_coverdat[working_coverdat$site_code == 'jena.de',])

working_coverdat <- working_coverdat %>% 
  filter(!(site_code == 'jena.de' & substr(date,6,7) %in% c('05','06')))


### bnch.us: per Seabloom combine all Carex into single taxa
working_coverdat %>%
  filter(site_code == 'bnch.us', substr(Taxon,1,5) == 'CAREX') %>% 
    count(year,Taxon) 

working_coverdat[working_coverdat$site_code == 'bnch.us' & substr(working_coverdat$Taxon,1,5) == 'CAREX',]$Taxon <- 'CAREX SP.'

working_coverdat %>%
  filter(site_code == 'hopl.us', substr(Taxon,1,5) == 'AVENA') %>% 
  count(year,Taxon) 

# hopl.us	all	Avena barbata	I would check to see whether the zeroes are due to misidentifying as A. fatua or being labeled as Avena species
# plot 2

working_coverdat %>%
  filter(site_code == 'hopl.us', substr(Taxon,1,5) == 'AVENA',plot == 2) %>% 
  select(plot, year_trt, Taxon, max_cover) %>% 
  print(n=60)

### sier.us	all	Bromus diandrus	We had problems early on separating B. diandrus and B. sterilis. You might need to combine these species, due uncertainty early in experiment

working_coverdat %>%
  filter(site_code == 'sier.us', Taxon %in%  c('BROMUS DIANDRUS','BROMUS STERILIS')) %>% 
  count(year,Taxon) %>% 
  print(n=60)

working_coverdat %>%
  filter(site_code == 'sier.us', Taxon %in%  c('BROMUS DIANDRUS','BROMUS STERILIS'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)

### for now we will combine both into BROMUS SP. and later take the maximum value of duplicates (needed to get around Jena fix - see below)
working_coverdat[working_coverdat$site_code == 'sier.us' & working_coverdat$Taxon %in%  c('BROMUS DIANDRUS','BROMUS STERILIS'),]$Taxon <- 'BROMUS SP.'


# doane.us	8,10, 19	Sporabolus compositus	This was from the pre-treatment year from summer of 2012 when Nebraska was facing a drought, so I wonder if I misidentified this species as Sporabolus, but it was Andropogon gerardii, which is more common there.	

working_coverdat %>%
  filter(site_code == 'doane.us', Taxon %in%  c('SPOROBOLUS COMPOSITUS','ANDROPOGON GERARDII')) %>% 
  count(year,Taxon) %>% 
  print(n=60)
# sporobolus completely disappears in Y2

### for now I'll turn SC into AG and later take the maximum value of duplicates (needed to get around Jena fix - see below)
working_coverdat[working_coverdat$site_code == 'doane.us' & working_coverdat$Taxon %in%  c('SPOROBOLUS COMPOSITUS','ANDROPOGON GERARDII'),]$Taxon <- 'ANDROPOGON GERARDII'


### cbgb.us	61	Schizachyrium scoparium	It should be there, but I checked the original field dataasheet for 2012 (year_trt = 3) and it is not listed, probably overlooked at that time

working_coverdat %>%
  filter(site_code == 'cbgb.us', substr(Taxon,1,5) == 'SCHIZ', plot == 61, year_trt %in% c(2,3,4)) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover)

working_coverdat <- rbind(working_coverdat,
                          working_coverdat %>%
  filter(site_code == 'cbgb.us', substr(Taxon,1,5) == 'SCHIZ', plot == 61, year_trt ==2) %>% 
  mutate(year_trt = 3, max_cover = 60))

working_coverdat %>%
  filter(site_code == 'cbgb.us', substr(Taxon,1,5) == 'SCHIZ', plot == 61, year_trt %in% c(2,3,4)) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover, local_lifespan)

#### take only maximum cover value when multiple are present in a year (for some taxa fixes and sites with multiple sampling dates)

cover_max <- working_coverdat %>% 
  group_by(site_code,year_trt,plot,Taxon) %>% 
  slice_max(order_by = max_cover,with_ties = FALSE)

working_coverdat %>%
  filter(site_code == 'sier.us', Taxon %in%  c('BROMUS SP.'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)

cover_max %>%
  filter(site_code == 'sier.us', Taxon %in%  c('BROMUS SP.'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)
### looks good in this example


### known site with two sampling dates per year
working_coverdat %>%
  filter(site_code == 'bayr.de', Taxon %in%  c('FESTUCA RUBRA'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)


cover_max %>%
  filter(site_code == 'bayr.de', Taxon %in%  c('FESTUCA RUBRA'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)

## now apply Bakker fixes - previous script ID'ed that 8 initial dominants were in this fix (very low number)

cover_max_fix <- left_join(cover_max,cover_fix %>% dplyr::select(site_code,Taxon,NewTaxon),
                              by = c('site_code','Taxon')) %>% 
  mutate(NewTaxon = ifelse(is.na(NewTaxon), Taxon, NewTaxon)) %>%
  ungroup() %>% 
  dplyr::select(-Taxon) %>% 
  rename(Taxon = NewTaxon) %>% 
  group_by(site_code, block, plot, subplot, trt, year, year_trt, max_year, Taxon,live) %>%
  summarize(max_cover = sum(max_cover), .groups = "keep")

cover_max <- left_join(cover_max_fix,cover_max %>% 
                         ungroup() %>% 
                         dplyr::select(Taxon,local_provenance,functional_group,local_lifespan) %>% 
                         distinct(Taxon,.keep_all = TRUE),
                       by='Taxon')

### uses relative cover here

working_coverdat <- cover_max %>% 
  group_by(site_code, year_trt, plot) %>% 
  mutate(rel_cover = max_cover / sum(max_cover)) %>%  
  mutate(cover = rel_cover) %>%  #use relative cover
  #mutate(cover = max_cover) %>%  # use absolute cover
  mutate(rank = min_rank(desc(rel_cover)), rich = n()) %>% 
  mutate(perc_rank = (rich + 1 - rank)/rich) %>% 
  filter(year_trt >= 0) %>% 
  group_by(site_code) %>% 
  mutate(nplots = n_distinct(plot)) %>% 
  ungroup() %>%   
  group_by(site_code, Taxon,year_trt) %>% 
  mutate(nplot = n_distinct(plot)) %>% 
  mutate(plotfreq_yr = nplot/nplots)

#224758

### calculate the frequency of a species in all pre-treatment plots
cover_pre <- working_coverdat[working_coverdat$year_trt == 0 & working_coverdat$live == 1,] %>% 
  filter(year_trt == 0) %>% 
  group_by(site_code, Taxon) %>% 
  mutate(avgCover = mean(cover)) %>% 
  mutate(avgRankPerc = mean(perc_rank)) %>% 
  mutate(plotfreq = plotfreq_yr) %>% 
  ungroup

range(cover_pre$avgCover)
hist(cover_pre$avgCover)

range(cover_pre$avgRankPerc)
hist(cover_pre$avgRankPerc)


### create initial dominance data.frame

### this code selects initial dominants based on single most dominant species
initial_dominants <- working_coverdat %>%
  filter(year_trt == 0) %>%
  group_by(site_code,plot)  %>%
  slice_max(order_by = max_cover,with_ties = TRUE) %>%   #2947
  mutate(initial_rel_cover = cover, initial_abs_cover = max_cover) %>% 
  dplyr::select(site_code,plot,Taxon,initial_rel_cover,initial_abs_cover) 
  # summarize(Taxon = Taxon[which.max(cover)], initial_abs_cover = max_cover[which.max(cover)], 
  #           initial_rel_cover = cover[which.max(cover)])
  # 2715 with_ties = FALSE
  # ~200 ties. Hmm

initial_dominants %>% mutate(id = paste0(site_code,'.',plot)) %>% 
  group_by(id) %>% filter(n() > 1) %>% 
  dplyr::select(site_code,plot,Taxon,initial_rel_cover) %>% print(n=400) 
## certainly seems that some sites are more likely to have 'ties' than others.
## Trying to cast forward to select 'dominant' introduces bias in subsequent analysis, but is that problematic?

initial_dominants %>% mutate(id = paste0(site_code,'.',plot)) %>% 
  group_by(id) %>% filter(n() > 1) %>% 
  left_join(.,working_coverdat %>% filter(year_trt == 1) %>% 
              #mutate(id = paste0(site_code,'.',plot)) %>% 
              ungroup() %>% 
              dplyr::select(site_code,plot,Taxon,max_cover),
            by = c('site_code','plot','Taxon')
            ) %>% 
  slice_max(order_by = max_cover,with_ties = TRUE) %>% 
  group_by(id) %>% filter(n() > 1) %>% 
  dplyr::select(site_code,plot,Taxon,initial_rel_cover) %>% print(n=400) 
### 19 repeats remain in year 2
### Going to keep repeats for now. Will need to repair code down the line

### are any of these species in Jon Bakker's taxonomy fixing code?
cover_fix %>% inner_join(.,initial_dominants,by=c('site_code','Taxon')) %>% 
                           distinct(Taxon,.keep_all = TRUE)
# now fixed earlier so zero
#8 species (35 plots due to repeats)
# much more manageable than expected!



# initial_dominants <- left_join(initial_dominants,cover_fix %>% dplyr::select(site_code,Taxon,NewTaxon),
#                                by = c('site_code','Taxon')) %>% 
#   mutate(NewTaxon = ifelse(is.na(NewTaxon), Taxon, NewTaxon)) %>%
#   select(-Taxon) %>% 
#   rename(Taxon = NewTaxon)


## merge to multi-year data and expand grid to include zero for years where spp is missing
dominant_year <- droplevels(initial_dominants[,c('site_code','plot','Taxon')]) %>% 
  left_join(working_coverdat[,c('site_code','plot','Taxon','cover','max_cover','rank','perc_rank','year_trt','max_year')],by = c('site_code','plot','Taxon')) %>% 
  group_by(site_code, plot) %>% 
  complete(year_trt = 0:max_year,nesting(Taxon,max_year), fill = list(cover =0, max_cover = 0, perc_rank=0)) %>%  # this adds about 3000 entries (from 16,000), which seems reasonable but needs a little testing
  filter(!is.na(site_code))


## some sites skipped a sampling year, need to weed those years out (will be all zeros now)
zero_year <- dominant_year %>% 
  group_by(site_code,year_trt) %>% 
  summarize(total_cov = sum(cover)) %>% 
  filter(total_cov == 0) %>% 
  mutate(id = paste(site_code,year_trt,sep='_'))

dominant_year <- dominant_year %>% 
  mutate(id = paste(site_code,year_trt,sep='_')) %>% 
  filter(!id %in% zero_year$id) %>% 
  left_join(
    (working_coverdat[,c('site_code','Taxon','year_trt','plotfreq_yr')] %>% 
       mutate(temp_id = paste(site_code,Taxon, year_trt,sep='_')) %>% 
       filter(!duplicated(temp_id))), by = c('site_code','year_trt','Taxon')) %>% 
  mutate(plotfreq_yr = ifelse(is.na(plotfreq_yr),0,plotfreq_yr))


# 19,822 (outdated) with single initial dominants (these were randomly chosen so not that valid)
# 22,597 with multiple initial dominants and Bakker fixs

scale_col <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

comb_soil <- comb %>% 
  dplyr::select(c('site_code','year_trt',"pct_C",'pct_N','ppm_P','ppm_K','ppm_Ca')) %>% 
  filter(year_trt == 0) %>% 
  mutate(c = scale_col(pct_C),n = scale_col(pct_N),p = scale_col(ppm_P),
         k = scale_col(ppm_K), ca = scale_col(ppm_Ca)) %>%
  group_by(site_code) %>% 
  summarize(c = var(c, na.rm = TRUE),
            n = var(n, na.rm = TRUE),
            p = var(p, na.rm = TRUE),
            k = var(k, na.rm = TRUE),
            ca = var(ca, na.rm = TRUE)) %>% 
  mutate(soil_var = rowMeans(dplyr::select(.,c(c,n,p,k,ca)),na.rm = TRUE))

### put it all together
dominant_pop <- initial_dominants %>% 
  #mutate(id = paste(site_code,plot,Taxon,sep='.')) %>% 
  left_join(.,dominant_year[,c('site_code','plot','Taxon','max_cover','cover','perc_rank','year_trt','plotfreq_yr')],
            by=c('site_code'='site_code','plot'='plot','Taxon'='Taxon')) %>% 
  left_join(.,cover_pre[, c('site_code','plot','trt','Taxon',"local_provenance", "local_lifespan", "functional_group",'plotfreq','avgCover')],
            by=c('site_code'='site_code','plot'='plot','Taxon'='Taxon')) %>% 
  left_join(.,comb[!duplicated(comb$site_code),
                   c('site_code','site_richness',"MAP_VAR_v2","MAP_v2")],
            by = 'site_code') %>% 
  left_join(.,comb[comb$year_trt == 0,c('site_code','plot','vascular_live_mass','litter_mass')], by = c('site_code','plot')) %>% 
    rename(initial_live_mass = vascular_live_mass, initial_litter_mass = litter_mass) %>% 
  left_join(.,comb_soil[,c('site_code','soil_var')],by='site_code') %>% 
  filter(trt %in% c('Control','NPK','Fence','NPK+Fence')) %>% 
  mutate(NPK = if_else(trt %in% c('NPK','NPK+Fence'),1,0), Fence = if_else(trt %in% c('Fence','NPK+Fence'),1,0))


#domMissByYear <- dominant_pop[dominant_pop$year1_change == 0,]

dominant_pop[dominant_pop$local_lifespan == 'NULL',] # 293

#what's going on with Lolium multiflorum in frue.ch?
working_coverdat %>% filter(site_code == 'frue.ch' & plot == 20) %>%
  dplyr::select(year, Taxon, max_cover) %>% arrange(Taxon) %>%    print(n=50)
### looks like it really might disappear in year 1

unique(dominant_pop$Taxon[dominant_pop$local_lifespan == 'NULL']) # 5


dominant_pop[dominant_pop$Taxon %in% c('NARDUS STRICTA','LOLIUM MULTIFLORUM',"ELYMUS SPICATUS"),]$local_lifespan <- 'PERENNIAL'
dominant_pop[dominant_pop$Taxon %in% c('DAUCUS CAROTA'),]$local_lifespan <- 'BIENNIAL'


#for now. poa sp. at koffler is likely perennial (poa pratensis or trivialis)
dominant_pop$local_lifespan <- ifelse(dominant_pop$local_lifespan %in% c('PERENNIAL','INDETERMINATE','NULL'), 'PERENNIAL','ANNUAL')
dominant_pop$functional_group <- ifelse(dominant_pop$functional_group %in% c('GRASS','GRAMINOID'), 'GRAMINOID',dominant_pop$functional_group)

### 7778 as of 2022-02-02
### 8908 with multiple initial dominants as of 2022-07-30

write.csv(dominant_pop,
          paste0('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NutNet manuscripts/Initial dominance/Project/initial-dominance/Data/Dominants-through-time_',
          Sys.Date(),'.csv'),
          row.names = F)
# 
# #### Community changes ####
# ### calculate RAC curves of communities INDEPENDENT of the initial dominant
# 
# working_coverdat$id <- paste(working_coverdat$site_code,working_coverdat$plot,sep='_')
# 
# rich_subord <- working_coverdat %>%
#   filter(trt %in% c('Control','NPK','Fence','NPK+Fence')) %>%
#   left_join(.,dominant_pop %>% mutate(dom = 'YES') %>%
#               dplyr::select(site_code,plot,year_trt,Taxon,dom),
#             by = c('site_code','plot','year_trt','Taxon')) %>%
#   mutate(dom = if_else(is.na(dom),'NO',dom)) %>%
#   filter(!dom == 'YES') %>%
#   group_by(site_code,plot,year_trt) %>%
#   summarize(subord_rich = n())
# 
# cover_subord <- working_coverdat %>% 
#   filter(trt %in% c('Control','NPK','Fence','NPK+Fence')) %>% 
#   left_join(.,dominant_pop %>% mutate(dom = 'YES') %>% dplyr::select(site_code,plot,Taxon,dom),
#             by = c('site_code','plot','Taxon')) %>% 
#   mutate(dom = if_else(is.na(dom),'NO',dom)) %>%
#   filter(!dom == 'YES')
# 
# rac_out <- RAC_change(df = cover_subord,
#                       species.var = 'Taxon',
#                       abundance.var = 'max_cover',
#                       replicate.var = 'id',
#                       time.var = 'year_trt',
#                       reference.time = 0)
# 
# rich_subord %>% filter(site_code == 'cdcr.us' & plot == 1) %>% dplyr::select(year_trt,subord_rich) %>% print(n=1000)
# working_coverdat %>% filter(site_code == 'cdcr.us' & plot == 1) %>% dplyr::select(year_trt,trt,Taxon,max_cover) %>% arrange(year_trt,max_cover) %>% print(n=1000)
# 
# cor(rac_out[!is.na(rac_out$evenness_change),4:8])
# 
# ggpairs(rac_out[!is.na(rac_out$evenness_change),4:8],upper = list(continuous = wrap("cor", size = 10)))
# 
# dom_comm <- dominant_pop %>% 
#   filter(year_trt == 0) %>% 
#   dplyr::select(-year_trt) %>% 
#   mutate(id = paste(site_code,plot, sep = '_')) %>% 
#   left_join(rac_out, by='id') %>% 
#   filter(trt %in% c('Control','NPK','Fence','NPK+Fence')) %>% 
#   mutate(NPK = if_else(trt %in% c('NPK','NPK+Fence'),1,0), Fence = if_else(trt %in% c('Fence','NPK+Fence'),1,0)) %>% 
#   mutate(ppt_var = as.numeric(MAP_VAR_v2)) %>% 
#   dplyr::select(-year_trt) %>% rename(year_trt = year_trt2) %>% 
#   left_join(.,comb %>% 
#               dplyr::select(site_code,plot,year_trt,vascular_live_mass,unsorted_live_mass,rich,evenness),
#             by=c('site_code','plot','year_trt')) %>% 
#   left_join(.,rich_subord %>% dplyr::select(site_code,plot,year_trt,subord_rich),
#             by=c('site_code','plot','year_trt'))
# 
# dominant_pop %>% 
#   filter(year_trt == 0) %>% 
#   mutate(id = paste(site_code,plot, sep = '_')) %>% 
#   filter(!id %in% rac_out$id)
# 
# working_coverdat %>% filter(id == 'ethamc.au_1') %>% dplyr::select(year_trt,Taxon) %>% print(n=100)
# # so plots with only one species in pre-treatment are getting filtered out here. 29 total, going to go forward without them
# 
# # write.csv(dom_comm,
# #           paste0('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NutNet manuscripts/Initial dominance/Project/initial-dominance/Data/RAC-subordinate_',Sys.Date(),'.csv'),
# #           row.names = F)


