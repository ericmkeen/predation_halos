setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(swfscMisc)
library(tidyr)
library(ggtext)

# Trip summary key
(trips <- readr::read_delim('tripSummary.txt', delim=';'))
trips %>% as.data.frame %>% head
trips_all <- trips
trips %>% nrow
trips$year %>% range

trips$duration %>% hist
trips$maxdist %>% range
sort(trips$maxdist)

tsumm <-
  trips %>%
  #mutate(colony_size = round(colony_size)) %>%
  group_by(colony_size) %>%
  summarize(site= site[1],
            n_sites = length(unique(site)),
            range = median(maxdist, na.rm=TRUE),
            range80 = quantile(maxdist, .9, na.rm=TRUE),
            range_mean = mean(maxdist, na.rm=TRUE),
            range_sd = sd(maxdist, na.rm=TRUE),
            hours = median(duration, na.rm=TRUE),
            hours_mean = mean(duration, na.rm=TRUE),
            hours_sd = sd(duration, na.rm=TRUE)) %>%
  mutate(range_cv = range_sd / range_mean) %>%
  mutate(hours_cv = hours_sd / hours_mean) %>%
  mutate(daily_trips = (18/2) / hours_mean)

tsumm %>% head(15)
tsumm %>% tail(15)
tsumm %>% arrange(range_mean) %>% head
tsumm %>% arrange(range_mean) %>% tail
tsumm %>% arrange(range_mean) %>% pull(range_mean)

ggplot(tsumm, aes(x=colony_size,
                  y=range80)) +
  geom_point() +
  scale_x_continuous(trans='log') +
  scale_y_continuous(trans='log')

ggplot(tsumm, aes(x=colony_size,
                  y=hours)) +
  geom_point() +
  scale_x_continuous(trans='log') +
  scale_y_continuous(trans='log')

ggplot(tsumm, aes(x=colony_size,
                  y=daily_trips)) +
  geom_point() +
  scale_x_continuous(trans='log') +
  scale_y_continuous(trans='log')

(trip2join <- tsumm)


################################################################################
################################################################################
################################################################################
# TRACKS

load('tracks.rds')

tracks <-
  tracks %>%
  mutate(id = paste(site, species, bird_id, deployment, tripID, sep='_')) %>%
  left_join(trip2join, by='site') %>%
  mutate(colony_size = round(colony_size/100))

tracks %>% head

ggplot(tracks,
       aes(x=lon,
           y=lat,
           group=id)) +
  geom_path(alpha=.2) +
  geom_point(mapping=aes(x=dep_lon,
                         y=dep_lat),
             color='dodgerblue', alpha=.2) +
  facet_wrap(~colony_size, scales='free', ncol=8) +
  xlab(NULL) + ylab(NULL) +
  labs(title='GPS tracks for each colony size sampled (hundreds of pairs)') +
  theme_light()

ggsave('fig_tracks.png', width=12, height=6, bg='white')

################################################################################
################################################################################
################################################################################
# Poop track examples

poop_map <- function(i){
  (lf <- list.files('poop_tracks2/'))
  #(i = sample(1:length(lf), size=1))
  message(i)
  (lfi <- paste0('poop_tracks2/',lf[i]))
  load(lfi)
  trackexp %>% head
  ptit <-
    paste0('Site: ', trackexp$site[1],' (', trackexp$year[1], ')')
  #ptit <-
  #  paste0('Site: ', trackexp$site[1],' (',
  #       lubridate::date(trackexp$time[1]), ')')

  psubtit <-
    paste0(trackexp$species[1], ' bird ID = ',
           trackexp$bird_id[1], ', trip ID = ',
           trackexp$tripID[1], ')')

  p <-
    ggplot(trackexp) +
    geom_path(mapping=aes(x=x, y)) +
    geom_point(mapping=aes(x=x[1], y=y[1]), pch=16, size=1.5) +
    geom_point(mapping=aes(x=x[nrow(trackexp)], y=y[nrow(trackexp)]), pch=1, size=1.5) +
    geom_point(mapping=aes(x=dep_lon[1], y=dep_lat[1]), pch=15, color='blue') +
    geom_point(data=trackexp %>% filter(excreta > 0),
               mapping=aes(x=x,
                           y=y,
                           size=excreta), alpha=.2, color='firebrick') +
    theme_light() + xlab(NULL) + ylab(NULL) +
    labs(title = ptit,
         subtitle= psubtit,
         size = expression("Excreta\n(grams " *sec^-1* ")"))

  return(p)
}

(lf <- list.files('poop_tracks2/'))
poop_map(sample(1:length(lf), size=1))
lf[4303]
4264
poop_map(4264)

ggarrange(poop_map(3871),
          poop_map(3510),
          poop_map(4712),
          poop_map(4370),
          poop_map(5027),
          poop_map(4404),
          nrow=2, ncol=3,
          labels = 'auto')

ggsave('fig_poop_maps.png', width=12, height=5, bg='white')


################################################################################
################################################################################
# Splashdowns

(lf <- list.files('poop_tracks2/'))
i=1
splashdowns <- data.frame()
for(i in 1:length(lf)){
  message(i,' out of ', length(lf),' ...')
  (lfi <- paste0('poop_tracks2/',lf[i]))
  load(lfi)
  energetics %>% as.data.frame %>% head
  trackexp %>% as.data.frame %>% head
  splashdown_secs <- length(which(trackexp$status == 'splashdown'))
  tot_secs <- nrow(trackexp)
  maxdist <- trackexp$km %>% max
  (spi <- data.frame(trackexp %>% select(site:band, bird_id, tripID) %>% head(1),
                     tot_secs, maxdist, splashdown_secs))
  splashdowns <- rbind(splashdowns, spi)
}

splashdowns$splashdown_secs %>% table

################################################################################
################################################################################
################################################################################
# Gut status on return flight at 1 km

if(file.exists('guts.rds')){
  load('guts.rds')
}else{

  guts <- data.frame()
  (lf <- list.files('poop_tracks2/'))
  trips$colony_size %>% sort %>% table
  trips$maxdist_linear %>% max
  (kms <- data.frame(km=0:186))
  i=3
  for(i in 1:length(lf)){
    message('Preparing gut contents check at 1km for poop track ',i,' out of ', length(lf),' ...')
    (lfi <- paste0('poop_tracks2/',lf[i]))
    load(lfi)
    energetics %>% as.data.frame %>% head
    trackexp %>% as.data.frame %>% head
    trackexp$site %>% unique
    trackexp$bird_id %>% unique
    trackexp$tripID %>% unique
    trackexp$deployment %>% unique
    #tripsumm
    trackexp %>% as.data.frame %>% head
    #trackexp$coldist %>% plot
    (trackexp$cum_ingesta <- cumsum(trackexp$ingesta))
    #trackexp$cum_ingesta %>% plot
    (trackexp$cum_excreta <- cumsum(trackexp$excreta))
    #trackexp$cum_excreta %>% plot
    (trackis <- trackexp[trackexp$q > 0.9,])
    (matchi <- which.min(abs(trackis$km - 1)))
    (mindist <- trackis$km[matchi])
    grams_in_gut_at_min <- trackis$gut[matchi]
    quantile_at_min <- trackis$q[matchi]
    cum_ingesta <- trackis$cum_ingesta[matchi]
    cum_excreta <- trackis$cum_excreta[matchi]
    status_at_min <- trackis$status[matchi]
    (tripsumm <- data.frame(energetics,
                            min_dist_on_return = mindist,
                            quantile_at_min,
                            status_at_min,
                            grams_in_gut_at_min,
                            tot_ingested_by_min = cum_ingesta,
                            tot_excreted_by_min = cum_excreta))
    tripsumm %>% head
    guts <- rbind(guts, tripsumm)
  }
  save(guts, file='guts.rds')
}

load('guts.rds')
guts %>% head

guts$FMR %>% unique %>% length
guts$sst %>% mean
guts$sst %>% sd
guts$sst %>% min
guts$sst %>% max

guts$FMR %>% mean
guts$FMR %>% sd
guts$FMR %>% min
guts$FMR %>% max

guts$tot_flight_costs %>% mean
guts$tot_flight_costs %>% median
guts$tot_flight_costs %>% sd

ggplot(guts, aes(x=colony_size, y=tot_flight_costs)) +
  geom_point() + geom_smooth() + scale_x_continuous(trans='log')
ggplot(guts, aes(x=colony_size, y=kJ_per_day)) +
  geom_point() + geom_smooth() + scale_x_continuous(trans='log')
ggplot(guts, aes(x=colony_size, y=food_per_trip)) +
  geom_point() + geom_smooth() + scale_x_continuous(trans='log')
ggplot(guts, aes(x=colony_size, y=n_trips)) +
  geom_point() + geom_smooth() + scale_x_continuous(trans='log')

guts %>% filter(food_per_trip < 2000) %>% pull(food_per_trip) %>% hist(breaks=20)
guts$duration %>% sort %>% round %>% table
guts$n_trips %>% hist
ggplot(guts, aes(x=duration, y=n_trips)) + geom_point() + scale_x_continuous(trans='log')

guts <-
  guts %>%
  filter(status_at_min == 'fly') %>%
  filter(min_dist_on_return < 5) %>%
  mutate(gut_ing_prop = grams_in_gut_at_min / tot_ingested_by_min) %>%
  mutate(gut_exc_prop = grams_in_gut_at_min / tot_excreted_by_min) #%>%
#left_join(trip2join, by='colony_size') %>%

guts %>% nrow
guts %>% head

gutsumm <-
  guts %>%
  group_by(colony_size) %>%
  summarize(n = n(),
            daily_trips = mean(n_trips, na.rm=TRUE),
            grams_in_gut_mean = mean(grams_in_gut_at_min, na.rm=TRUE)*excreta_ingesta_ratio[1],
            grams_in_gut_median = median(grams_in_gut_at_min, na.rm=TRUE)*excreta_ingesta_ratio[1],
            grams_in_gut_sd = sd(grams_in_gut_at_min, na.rm=TRUE)*excreta_ingesta_ratio[1],
            gut_ing_prop_mean = mean(gut_ing_prop, na.rm=TRUE),
            gut_ing_prop_median = median(gut_ing_prop, na.rm=TRUE),
            gut_ing_prop_sd = sd(gut_ing_prop, na.rm=TRUE)) %>%
  mutate(gut_ing_prop_se = gut_ing_prop_sd / sqrt(n)) %>%
  mutate(mtons_in_gut_mean = (colony_size * daily_trips * grams_in_gut_mean)/(1000*1000*2),
         mtons_in_gut_median = (colony_size * daily_trips * grams_in_gut_median)/(1000*1000*2),
         mtons_in_gut_sd = (colony_size * daily_trips * grams_in_gut_sd)/(1000*1000*2)) %>%
  mutate(mtons_in_gut_se = mtons_in_gut_sd / sqrt(n)) %>%
  mutate(colony_size = round(colony_size / 100))

gutsumm %>% as.data.frame
gutsumm$gut_ing_prop_mean %>% range
gutsumm %>% arrange(colony_size) %>% pull(mtons_in_gut_mean)
{
  p1 <-
    ggplot(gutsumm,
           aes(x=colony_size,
               y=gut_ing_prop_mean)) +
    geom_point(size=2) +
    geom_segment(mapping=aes(x=colony_size,
                             xend=colony_size,
                             y=gut_ing_prop_mean - gut_ing_prop_se,
                             yend=gut_ing_prop_mean + gut_ing_prop_se),
                 lwd=.5, alpha=.3) +
    scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
    scale_y_continuous(breaks=seq(0, 1, by=.1), limits=c(-.075, 1.075)) +
    theme_light() +
    labs(title='Fraction of ingesta still in gut upon return') +
    xlab('Colony size (hundreds of pairs)') +
    ylab('Proportion')

  p2 <-
    ggplot(gutsumm,
           aes(x=colony_size,
               y=mtons_in_gut_mean)) +
    geom_point() +
    geom_segment(mapping=aes(x=colony_size,
                             xend=colony_size,
                             y=mtons_in_gut_mean - mtons_in_gut_se,
                             yend=mtons_in_gut_mean + mtons_in_gut_se),
                 #y=kg_in_gut_mean - kg_in_gut_sd,
                 #yend=kg_in_gut_mean + kg_in_gut_sd),
                 lwd=.5, alpha=.3) +
    scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
    scale_y_continuous(trans='log', breaks=c(.05, .25, 1, 4, 15, 50, 150)) +
    labs(title='Total ingesta excreted at colony per day') +
    xlab('Colony size (hundreds of pairs)') +
    ylab('log Excreta (metric tons)') +
    theme_light()

  ggarrange(p1, p2, nrow=1, labels='auto')
}

ggsave('fig_gut_at_return.png', width=10, height=4, bg='white')


################################################################################
################################################################################
################################################################################
# Halo & redistribution summaries

if(file.exists('halo_summaries.rds')){
  load('halo_summaries.rds')
}else{
  mrs <- data.frame()
  (lf <- list.files('poop_tracks2/'))
  trips$colony_size %>% sort %>% table
  trips$maxdist_linear %>% max
  (kms <- data.frame(km=0:186))
  i=3
  # troubleshooting
  #(is <- grep('20120712-Diabas-03-4181306', lf))
  #lf[is]
  #(i = is[1])
  for(i in 1:length(lf)){
    message('Preparing halo summary for poop track ',i,' out of ', length(lf),' ...')
    (lfi <- paste0('poop_tracks2/',lf[i]))
    load(lfi)
    energetics %>% as.data.frame %>% head
    trackexp %>% as.data.frame %>% head

    # Add cumulative and quantil diagnostics
    suppressMessages({
      tripsumm <-
        trackexp %>%
        arrange(km) %>%
        mutate(sum_ing = sum(ingesta)) %>%
        mutate(sum_exc = sum(excreta)) %>%
        mutate(cum_ing = cumsum(ingesta)) %>%
        mutate(cum_exc = cumsum(excreta)) %>%
        rowwise %>%
        mutate(p_ing = cum_ing / sum_ing) %>%
        mutate(p_exc = cum_exc / sum_exc) %>%
        ungroup %>%
        group_by(site, deployment, species, band, bird_id, year, tripID) %>%
        summarize(km_ing_20 = ifelse(max(abs(p_ing))>=0.2,
                                     head(km[which(abs(p_ing) >= .2)], 1),
                                     NA),
                  km_exc_20 = ifelse(max(abs(p_exc))>=0.2,
                                     head(km[which(abs(p_exc) >= .2)], 1),
                                     NA),
                  km_ing_50 = ifelse(max(abs(p_ing))>=0.5,
                                     head(km[which(abs(p_ing) >= .5)], 1),
                                     NA),
                  km_exc_50 = ifelse(max(abs(p_exc))>=0.5,
                                     head(km[which(abs(p_exc) >= .5)], 1),
                                     NA),
                  km_ing_80 = ifelse(max(abs(p_ing))>=0.8,
                                     head(km[which(abs(p_ing) >= .8)], 1),
                                     NA),
                  km_exc_80 = ifelse(max(abs(p_exc))>=0.8,
                                     head(km[which(abs(p_exc) >= .8)], 1),
                                     NA),
                  g_ing_20 = sum(ingesta[km <= km_ing_20], na.rm=TRUE),
                  g_ing_50 = sum(ingesta[km <= km_ing_50], na.rm=TRUE),
                  g_ing_80 = sum(ingesta[km <= km_ing_80], na.rm=TRUE),
                  g_exc_20 = sum(ingesta[km <= km_exc_20], na.rm=TRUE),
                  g_exc_50 = sum(ingesta[km <= km_exc_50], na.rm=TRUE),
                  g_exc_80 = sum(ingesta[km <= km_exc_80], na.rm=TRUE)) %>%
        mutate(g_ing_20_80 = g_ing_80 - g_ing_20) %>%
        mutate(g_exc_20_80 = g_exc_80 - g_exc_20) %>%
        ungroup %>%
        as.data.frame
    })
    tripsumm %>% head

    tripsumm <- data.frame(energetics %>% select(site:n_trips),
                           tripsumm %>% select(km_ing_20:g_exc_20_80))
    tripsumm %>% head
    tripsumm %>% tail
    mrs <- rbind(mrs, tripsumm)
  }
  # Save results
  save(mrs, file='halo_summaries.rds')
}

################################################################################
################################################################################
################################################################################
# Summarize those summaries

load('halo_summaries.rds') # mrs
mrs %>% head

mrs$km_ing_50 %>% hist
mrs$km_exc_50 %>% hist

mrs$colony_size %>% unique %>% sort

mrssum <-
  mrs %>%
  mutate(id = paste(site, species, bird_id, deployment, tripID, sep='_')) %>%
  mutate(rim_width = km_ing_80 - km_ing_20) %>%
  group_by(colony_size) %>%
  summarize(n=n(),
            dur_mean = duration %>% mean(na.rm=TRUE),
            dur_median = duration %>% median(na.rm=TRUE),
            dur_sd = duration %>% sd(na.rm=TRUE),

            trips_mean = n_trips %>% mean(na.rm=TRUE),
            trips_median = n_trips %>% median(na.rm=TRUE),
            trips_sd = n_trips %>% sd(na.rm=TRUE),

            ing20_mean = g_ing_20 %>% mean(na.rm=TRUE),
            ing20_median = g_ing_20 %>% median(na.rm=TRUE),
            ing20_sd = g_ing_20 %>% sd(na.rm=TRUE),

            ing50_mean = g_ing_50 %>% mean(na.rm=TRUE),
            ing50_median = g_ing_50 %>% median(na.rm=TRUE),
            ing50_sd = g_ing_50 %>% sd(na.rm=TRUE),

            ing80_mean = g_ing_80 %>% mean(na.rm=TRUE),
            ing80_median = g_ing_80 %>% median(na.rm=TRUE),
            ing80_sd = g_ing_80 %>% sd(na.rm=TRUE),

            exc20_mean = g_exc_20 %>% mean(na.rm=TRUE),
            exc20_median = g_exc_20 %>% median(na.rm=TRUE),
            exc20_sd = g_exc_20 %>% sd(na.rm=TRUE),

            exc50_mean = g_exc_50 %>% mean(na.rm=TRUE),
            exc50_median = g_exc_50 %>% median(na.rm=TRUE),
            exc50_sd = g_exc_50 %>% sd(na.rm=TRUE),

            exc80_mean = g_exc_80 %>% mean(na.rm=TRUE),
            exc80_median = g_exc_80 %>% median(na.rm=TRUE),
            exc80_sd = g_exc_80 %>% sd(na.rm=TRUE),

            km20_mean = km_ing_20 %>% mean(na.rm=TRUE),
            km20_median = km_ing_20 %>% median(na.rm=TRUE),
            km20_sd = km_ing_20 %>% sd(na.rm=TRUE),

            km50_mean = km_ing_50 %>% mean(na.rm=TRUE),
            km50_median = km_ing_50 %>% median(na.rm=TRUE),
            km50_sd = km_ing_50 %>% sd(na.rm=TRUE),

            km80_mean = km_ing_80 %>% mean(na.rm=TRUE),
            km80_median = km_ing_80 %>% median(na.rm=TRUE),
            km80_sd = km_ing_80 %>% sd(na.rm=TRUE),

            rim_width_mean = (km_ing_80 - km_ing_20) %>% mean(na.rm=TRUE),
            rim_width_median = (km_ing_80 - km_ing_20) %>% median(na.rm=TRUE),
            rim_width_sd = (km_ing_80 - km_ing_20) %>% sd(na.rm=TRUE),

            ekm20_mean = km_exc_20 %>% mean(na.rm=TRUE),
            ekm20_median = km_exc_20 %>% median(na.rm=TRUE),
            ekm20_sd = km_exc_20 %>% sd(na.rm=TRUE),

            ekm50_mean = km_exc_50 %>% mean(na.rm=TRUE),
            ekm50_median = km_exc_50 %>% median(na.rm=TRUE),
            ekm50_sd = km_exc_50 %>% sd(na.rm=TRUE),

            ekm80_mean = km_exc_80 %>% mean(na.rm=TRUE),
            ekm80_median = km_exc_80 %>% median(na.rm=TRUE),
            ekm80_sd = km_exc_80 %>% sd(na.rm=TRUE),

            erim_width_mean = (km_exc_80 - km_exc_20) %>% mean(na.rm=TRUE),
            erim_width_median = (km_exc_80 - km_exc_20) %>% median(na.rm=TRUE),
            erim_width_sd = (km_exc_80 - km_exc_20) %>% sd(na.rm=TRUE)) %>%
  rowwise %>%
  # redistribution
  mutate(exc_ing = ekm50_median - km50_median) %>%
  # proportionals
  mutate(rim_width_prop = rim_width_median / km50_median,
         erim_width_prop = erim_width_median / km50_median,
         exc_ing_prop = exc_ing / km50_median) %>%
  # SE's
  mutate(trips_se = trips_sd / sqrt(n),
         dur_se = dur_sd / sqrt(n),
         ing20_se = ing20_sd / sqrt(n),
         ing50_se = ing50_sd / sqrt(n),
         ing80_se = ing80_sd / sqrt(n),
         exc20_se = exc20_sd / sqrt(n),
         exc50_se = exc50_sd / sqrt(n),
         exc80_se = exc80_sd / sqrt(n),
         km20_se = km20_sd / sqrt(n),
         km50_se = km50_sd / sqrt(n),
         km80_se = km80_sd / sqrt(n),
         ekm20_se = ekm20_sd / sqrt(n),
         ekm50_se = ekm50_sd / sqrt(n),
         ekm80_se = ekm80_sd / sqrt(n),
         rim_width_se = rim_width_sd / sqrt(n),
         erim_width_se = erim_width_sd / sqrt(n),
         exc_ing_se = mean(c(ekm50_se, km50_se)),
         exc_ing_prop_se = exc_ing_prop / sqrt(n)) %>%
  mutate(exc_ing_prop_se_max = min(c(1, exc_ing_prop + exc_ing_prop_se))) %>%
  mutate(exc_ing_prop_se_min = max(c(0, exc_ing_prop - exc_ing_prop_se))) %>%
  # areas km2
  mutate(area20 = pi*(km20_median^2)) %>%
  mutate(area50 = pi*(km50_median^2)) %>%
  mutate(area80 = pi*(km80_median^2)) %>%
  mutate(area_rim = area80 - area20) %>%
  # density per trip (g / km2)
  mutate(d_ing_trip_halo = ing80_mean / area80) %>%
  mutate(d_ing_trip_rim = (ing80_mean - ing20_mean) / area80) %>%
  mutate(d_exc_trip_halo = exc80_mean / area80) %>%
  mutate(d_exc_trip_rim = (exc80_mean - exc20_mean) / area80) %>%
  # density per bird (g / km2)
  mutate(d_ing_bird_halo = d_ing_trip_halo * trips_mean) %>%
  mutate(d_ing_bird_rim = d_ing_trip_rim * trips_mean) %>%
  mutate(d_exc_bird_halo = d_exc_trip_halo * trips_mean) %>%
  mutate(d_exc_bird_rim = d_exc_trip_rim * trips_mean) %>%
  # density per colony (kg / km2)
  mutate(d_ing_col_halo = (d_ing_bird_halo * colony_size * 2)/(1000)) %>%
  mutate(d_ing_col_rim = (d_ing_bird_rim * colony_size * 2)/(1000)) %>%
  mutate(d_exc_col_halo = (d_exc_bird_halo * colony_size * 2)/(1000)) %>%
  mutate(d_exc_col_rim = (d_exc_bird_rim * colony_size * 2)/(1000)) %>%
  # totals metric tons
  mutate(ing_col_halo = (ing80_mean*colony_size*2)/(1000*1000)) %>%
  mutate(ing_col_rim = ((ing80_mean - ing20_mean)*colony_size*2)/(1000*1000)) %>%
  mutate(exc_col_halo = (exc80_mean*colony_size*2)/(1000*1000)) %>%
  mutate(exc_col_rim = ((exc80_mean - exc20_mean)*colony_size*2)/(1000*1000)) %>%
  # wrap up
  mutate(colony_size = colony_size / 100) %>%
  ungroup %>%
  as.data.frame

mrssum %>% head
mrssum %>% arrange(desc(d_exc_col_rim))

################################################################################
################################################################################
################################################################################
# Halos & Redistribution

mrssum %>% head

{
  dodge <- .06
  dodge_se <- 0.01
  ing_color <- 'darkblue'
  exc_color <- 'darkorange4'
  cex_subt = 8
  cex_t = 12
  cex_cap = 8

  # Foraging halos increase
  p1 <-
    ggplot(mrssum) +
    # Ingesta halo
    geom_segment(mapping=aes(x=colony_size - (dodge + dodge_se)*colony_size,
                             y=km20_mean - km20_se,
                             yend=km20_mean + km20_se),
                 alpha=.7, lwd=.2, color=ing_color) +
    geom_segment(mapping=aes(x=colony_size - (dodge - dodge_se)*colony_size,
                             y=km80_mean - km80_se,
                             yend=km80_mean + km80_se),
                 alpha=.7, lwd=.2, color=ing_color) +
    geom_segment(mapping=aes(x=colony_size - dodge*colony_size,
                             y=km20_mean,
                             yend=km80_mean),
                 lwd=3, alpha=.7, color=ing_color) +
    # Excreta halo
    geom_segment(mapping=aes(x=colony_size + (dodge + dodge_se)*colony_size,
                             y=ekm20_mean - ekm20_se,
                             yend=ekm20_mean + ekm20_se),
                 alpha=.7, lwd=.2, color=exc_color) +
    geom_segment(mapping=aes(x=colony_size + (dodge - dodge_se)*colony_size,
                             y=ekm80_mean - ekm80_se,
                             yend=ekm80_mean + ekm80_se),
                 alpha=.7, lwd=.2, color=exc_color) +
    geom_segment(mapping=aes(x=colony_size + dodge*colony_size,
                             y=ekm20_mean,
                             yend=ekm80_mean),
                 lwd=3, alpha=.7, color=exc_color) +
    scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
    ylab('Distance from colony (km)') +
    xlab('log Colony pairs x 100') +
    theme_light() +
    labs(title='Foraging halos increase with colony size',
         subtitle= 'Halo rims shown for <span style="color:darkblue;">**feeding**</span> and <span style="color:darkorange4;">**excretion**</span>') +
    theme(plot.title = element_text(size = cex_t),
          plot.subtitle = element_markdown(hjust=0, size=cex_subt),
          plot.caption = element_text(size=cex_cap))

  # EXCRETION - INGESTION anomaly
  p2 <-
    ggplot(mrssum) +
    geom_hline(yintercept=0, lty=2, alpha=.4) +
    geom_segment(mapping=aes(x=colony_size,
                             y=exc_ing + exc_ing_se,
                             yend=exc_ing - exc_ing_se),
                 alpha=.5) +
    geom_point(mapping=aes(x=colony_size, y=exc_ing)) +
    scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
    scale_y_continuous(breaks=c(-5, 0, 5, 10, 15, 20), limits=c(-6, 22.5)) +
    xlab('log Colony pairs x 100') +
    ylab('Difference (km)') +
    theme_light() +
    labs(title='Excretion vs. Ingestion anomaly',
         subtitle='Median Excretion Distance - Median Ingestion Distance') +
    theme(plot.title = element_text(size = cex_t),
          plot.subtitle = element_markdown(hjust=0, size=cex_subt),
          plot.caption = element_text(size=cex_cap))

  # Proportional anomaly
  p3 <-
    ggplot(mrssum) +
    geom_point(mapping=aes(x=colony_size, y=exc_ing_prop)) +
    geom_segment(mapping=aes(x=colony_size,
                             y=exc_ing_prop_se_min,
                             yend=exc_ing_prop_se_max),
                 alpha=.5) +
    geom_hline(yintercept=c(0), lty=2, alpha=.2) +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
    labs(title='Proportional anomaly',
         subtitle='Anamoly scaled by median ingestion distance') +
    xlab('log Colony pairs x 100') +
    ylab('Proportional Difference') +
    theme_light() +
    theme(plot.title = element_text(size = cex_t),
          plot.subtitle = element_markdown(hjust=0, size=cex_subt),
          plot.caption = element_text(size=cex_cap))

  ggarrange(p1, p2, p3, nrow=1, labels='auto')
}

ggsave('fig_redistribution.png', width=12, height=4, bg='white')

################################################################################
################################################################################
################################################################################
# Logistics - Range relationships

mrssum %>% as.data.frame %>% head
mrssumi <- mrssum
mrssumi$logcol <- log(mrssumi$colony_size)
cs <- mrssumi$colony_size
logcol <- mrssumi$logcol

# area
(lm2 <- lm(log(area80) ~ logcol, data=mrssumi))
lm2 %>% summary
predict(lm2, newdata=mrssumi, se.fit=TRUE)
p2 <- exp(predict(lm2, newdata=mrssumi))
plot(p2 ~ mrssumi$logcol)
ggplot() + geom_point(mapping=aes(x=cs, y=p2)) +
  scale_x_continuous(trans='log') + scale_y_continuous(trans='log')

# excreta density per bird g per km2
(lm4 <- lm(log(d_exc_bird_halo) ~ logcol, data=mrssumi))
lm4 %>% summary
p4 <- exp(predict(lm4, newdata=mrssumi))
mrssum$d_exc_bird_halo %>% plot(ylim=c(0,35))
p4 %>% plot(ylim=c(0,35))
plot(p4 ~ mrssumi$logcol)
ggplot() + geom_point(mapping=aes(x=cs, y=p4)) +
  scale_x_continuous(trans='log') +
  scale_y_continuous(trans='log', breaks=c(.1, 12))

# excreta density kg per km2
p5 <- (((p4*cs*100*2)/1000))
p5 %>% plot(ylim=c(0, 100))
mrssum$d_exc_col_halo %>% plot(ylim=c(0, 100))
ggplot() + geom_point(mapping=aes(x=cs, y=p5)) +
  scale_x_continuous(trans='log') + scale_y_continuous(trans='log')


{
  mrssum %>% head
  predx <- min(mrssum$colony_size):max(mrssum$colony_size)
  predx
  gammax <- 1 #.9 # choose the gamma that minimizes AIC

  # 1: Foraging range ============================================================
  {
    m1<- mgcv::gam(log(km80_median) ~ s(log(colony_size)), data=mrssum, gamma=gammax)
    p1 <- data.frame(predx,
                     predict(m1, se.fit=TRUE,
                             newdata=data.frame(colony_size = predx)))
    p1poly <- data.frame(x = c(predx, rev(predx)),
                         y= c(c(p1$fit + p1$se.fit),
                              c(rev(p1$fit) - rev(p1$se.fit))))
    (pv <- summary(m1)$s.pv %>% round(3))
    if(pv == 0){pv <- 'p < 0.001'}else{pv <- paste0('p = ', pv)} ; pv
    (de <- paste0('(DE = ', (summary(m1)$dev.expl * 100) %>% round(0), '%)'))
    (subt1 <- paste0(pv,' ',de))

    plot1 <-
      ggplot(mrssumi) +
      geom_point(mapping=aes(x=colony_size,
                             y=km80_median),
                 alpha=.8) +
      geom_polygon(data=p1poly,
                   mapping = aes(x = x, y = exp(y)), fill='black', alpha=.2) +
      geom_path(data = p1, mapping = aes(x = predx, y = exp(fit)),
                color= 'black', alpha=.7, lwd=1) +
      ylab('log Range (km)') +
      labs(title = 'Foraging Range',
           subtitle = subt1) + # mean across trips
      scale_y_continuous(trans='log', breaks=c(.5, 2.5, 5, 10, 25, 50), limits =) +
      scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
      xlab('log Colony pairs x 100') +
      theme_minimal() +
      theme(plot.title = element_text(size = cex_t),
            plot.subtitle = element_markdown(hjust=0, size=cex_subt))

    plot1 %>% print
  }

  # 2: Foraging Area ============================================================
  {
    m2<- mgcv::gam(log(area80) ~ s(log(colony_size)), data=mrssum, gamma=gammax)
    p2 <- data.frame(predx,
                     predict(m2, se.fit=TRUE,
                             newdata=data.frame(colony_size = predx)))
    p2poly <- data.frame(x = c(predx, rev(predx)),
                         y= c(c(p2$fit + p2$se.fit),
                              c(rev(p2$fit) - rev(p2$se.fit))))
    (pv <- summary(m2)$s.pv %>% round(3))
    if(pv == 0){pv <- 'p < 0.001'}else{pv <- paste0('p = ', pv)} ; pv
    (de <- paste0('(DE = ', (summary(m2)$dev.expl * 100) %>% round(0), '%)'))
    (subt2 <- paste0(pv,' ',de))

    plot2 <-
      ggplot(mrssum) +
      geom_point(mapping=aes(x=colony_size,
                             y=area80),
                 alpha=.8,
                 pch=halo_pch) +
      geom_polygon(data=p2poly,
                   mapping = aes(x = x, y = exp(y)), #fill=halo_color,
                   alpha=.2) +
      geom_path(data = p2,
                mapping = aes(x = predx,
                              y = exp(fit)),
                lty = halo_lty,
                alpha=.8, lwd=1) +
      scale_y_continuous(trans='log', breaks=c(5, 20, 100, 800, 5000, 20000)) +
      scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
      ylab(expression("log Area (" *km^2* ")")) +
      labs(title = 'Foraging Area',
           subtitle = subt2) + # mean across trips
      xlab('log Colony pairs x 100') +
      theme_minimal()+
      theme(plot.title = element_text(size = cex_t),
            plot.subtitle = element_markdown(hjust=0, size=cex_subt))

    plot2 %>% print
  }

  # 3: Trip duration ============================================================
  {
    m3<- mgcv::gam(log(dur_median) ~ s(log(colony_size)), data=mrssum, gamma=gammax)
    p3 <- data.frame(predx,
                     predict(m3, se.fit=TRUE,
                             newdata=data.frame(colony_size = predx)))
    p3poly <- data.frame(x = c(predx, rev(predx)),
                         y= c(c(p3$fit + p3$se.fit),
                              c(rev(p3$fit) - rev(p3$se.fit))))
    (pv <- summary(m3)$s.pv %>% round(3))
    if(pv == 0){pv <- 'p < 0.001'}else{pv <- paste0('p = ', pv)} ; pv
    (de <- paste0('(DE = ', (summary(m3)$dev.expl * 100) %>% round(0), '%)'))
    (subt3 <- paste0(pv,' ',de))

    plot3 <-
      ggplot(mrssum) +
      geom_point(mapping=aes(x=colony_size,
                             y=dur_median),
                 alpha=.8) +
      geom_polygon(data=p3poly,
                   mapping = aes(x = x, y = exp(y)), fill='black', alpha=.2) +
      geom_path(data = p3, mapping = aes(x = predx, y = exp(fit)),
                color= 'black', alpha=.7, lwd=1) +
      ylab('Trip duration (hr)') +
      labs(title = 'Trip Duration',
           subtitle = subt3) + # mean
      scale_y_continuous(breaks=c(0, 3, 6, 9, 12), limits=c(0,12)) +
      scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
      xlab('log Colony pairs x 100') +
      theme_minimal()+
      theme(plot.title = element_text(size = cex_t),
            plot.subtitle = element_markdown(hjust=0, size=cex_subt))

    plot3 %>% print
  }

  # 4: Daily trips ============================================================
  {
    m4<- mgcv::gam(trips_median ~ s(log(colony_size)), data=mrssum, gamma=gammax)
    p4 <- data.frame(predx,
                     predict(m4, se.fit=TRUE,
                             newdata=data.frame(colony_size = predx)))
    p4poly <- data.frame(x = c(predx, rev(predx)),
                         y= c(c(p4$fit + p4$se.fit),
                              c(rev(p4$fit) - rev(p4$se.fit))))
    (pv <- summary(m4)$s.pv %>% round(3))
    if(pv == 0){pv <- 'p < 0.001'}else{pv <- paste0('p = ', pv)} ; pv
    (de <- paste0('dev. exp. = ', (summary(m4)$dev.expl * 100) %>% round(0), '%'))
    (subt4 <- paste0(pv,' ',de))

    plot4 <-
      ggplot(mrssum) +
      geom_point(mapping=aes(x=colony_size,
                             y=trips_median),
                 alpha=.8) +
      geom_polygon(data=p4poly,
                   mapping = aes(x = x, y = (y)), fill='black', alpha=.2) +
      geom_path(data = p4, mapping = aes(x = predx, y = (fit)),
                color= 'black', alpha=.7, lwd=1) +
      ylab('Daily Trips') +
      labs(title = expression('Possible Daily Trips bird'^'-1'),
           subtitle = subt4) + # max possible
      scale_y_continuous(breaks=c(1:6), limits=c(1, 6.5)) +
      scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
      xlab('log Colony pairs x 100') +
      theme_minimal()+
      theme(plot.title = element_text(size = cex_t),
            plot.subtitle = element_markdown(hjust=0, size=cex_subt))

    plot4 %>% print
  }

  ggarrange(plot1, plot2, plot3, plot4, nrow=1, labels='auto')
  }

ggsave('fig_logistics.png', width=12, height=4, bg='white')


################################################################################
################################################################################
################################################################################
# Excreta only

{
  mrssum %>% head
  predx <- min(mrssum$colony_size):max(mrssum$colony_size)
  predx
  gammax <- 1.4 #.625 #.9 # choose the gamma that minimizes AIC
  ing_halo_color <- 'darkblue'
  ing_rim_color <- 'dodgerblue'
  exc_halo_color <- 'darkorange4'
  exc_rim_color <- 'goldenrod2'
  halo_pch <- 16
  rim_pch <- 1
  halo_lty <- 1
  rim_lty <- 2
  span <- .95
  span2 <- .7
  cex_subt = 7
  cex_t = 12
  cex_cap = 8

  # 11: Total excreta =============================================================
  {
    m11 <- mgcv::gam(log(exc_col_halo) ~ s(log(colony_size)),
                     data=mrssum,
                     gamma=gammax)
    AIC(m11)
    p11 <- data.frame(predx,
                      predict(m11, se.fit=TRUE,
                              newdata=data.frame(colony_size = predx)))
    p11poly <- data.frame(x = c(predx, rev(predx)),
                          y= c(c(p11$fit + p11$se.fit),
                               c(rev(p11$fit) - rev(p11$se.fit))))
    (pv <- summary(m11)$s.pv %>% round(3))
    if(pv == 0){pv <- 'p < 0.001'}else{pv <- paste0('p = ', pv)} ; pv
    (de <- paste0('(DE = ', (summary(m11)$dev.expl * 100) %>% round(0), '%)'))
    (subt11 <- paste0(pv,' ',de))

    plot11 <-
      ggplot(mrssum) +
      geom_point(mapping=aes(x=colony_size,
                             y=exc_col_halo),
                 #alpha=.8,
                 pch=halo_pch) +
      geom_polygon(data=p11poly,
                   mapping = aes(x = x, y = exp(y)), alpha=.2) +
      geom_path(data = p11,
                mapping = aes(x = predx,
                              y = exp(fit)),
                alpha=.8, lwd=1) +
      scale_y_continuous(trans='log', breaks=c(.2, .6, 2, 8, 30, 140))+
      ylab('log Excreta (metric tons)') +
      labs(title = expression('Total Excreta day'^'-1'),
           subtitle=subt11) +
      scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
      xlab('log Colony pairs x 100') +
      theme_minimal()+
      theme(plot.title = element_text(size = cex_t),
            plot.subtitle = element_markdown(hjust=0, size=cex_subt))

    plot11 %>% print
  }

  # 12: Excreta density per bird =========================================================
  {
    m12 <- mgcv::gam(log(d_exc_bird_halo) ~ s(log(colony_size)),
                     data=mrssum,
                     gamma=gammax)
    AIC(m12)
    p12 <- data.frame(predx,
                      predict(m12, se.fit=TRUE,
                              newdata=data.frame(colony_size = predx)))
    p12poly <- data.frame(x = c(predx, rev(predx)),
                          y= c(c(p12$fit + p12$se.fit),
                               c(rev(p12$fit) - rev(p12$se.fit))))
    (pv <- summary(m12)$s.pv %>% round(3))
    if(pv == 0){pv <- 'p < 0.001'}else{pv <- paste0('p = ', pv)} ; pv
    (de <- paste0('(DE = ', (summary(m12)$dev.expl * 100) %>% round(0), '%)'))
    (subt12 <- paste0(pv,' ',de))

    plot12 <-
      ggplot(mrssum) +
      geom_point(mapping=aes(x=colony_size,
                             y=d_exc_bird_halo),
                 #alpha=.8,
                 pch=halo_pch) +
      geom_polygon(data=p12poly,
                   mapping = aes(x = x, y = exp(y)), alpha=.2) +
      geom_path(data = p12,
                mapping = aes(x = predx,
                              y = exp(fit)),
                alpha=.8, lwd=1) +
      scale_y_continuous(trans='log', breaks=c(.03, .1, .4, 1.5, 8, 40))+
      ylab(expression("log Density (g " *km^-2* ")")) +
      labs(title = expression(paste("Excreta Density bird"^"-1","day"^"-1")),
           subtitle=subt12) +
      scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
      xlab('log Colony pairs x 100') +
      theme_minimal()+
      theme(plot.title = element_text(size = cex_t),
            plot.subtitle = element_markdown(hjust=0, size=cex_subt))

    plot12 %>% print
  }

  # 13: Excreta density whole colony =========================================================
  {
    m13 <- mgcv::gam(log(d_exc_col_halo) ~ s(log(colony_size)),
                     data=mrssum,
                     gamma=gammax)
    AIC(m13)
    p13 <- data.frame(predx,
                      predict(m13, se.fit=TRUE,
                              #newdata=data.frame(logcol = log(predx))))
                              newdata=data.frame(colony_size = predx)))
    p13poly <- data.frame(x = c(predx, rev(predx)),
                          y= c(c(p13$fit + p13$se.fit),
                               c(rev(p13$fit) - rev(p13$se.fit))))
    (pv <- summary(m13)$s.pv %>% round(3))
    if(pv == 0){pv <- 'p < 0.001'}else{pv <- paste0('p = ', pv)} ; pv
    (de <- paste0('(DE = ', (summary(m13)$dev.expl * 100) %>% round(0), '%)'))
    (subt13 <- paste0(pv,' ',de))

    plot13 <-
      ggplot(mrssum) +
      #ggplot(mrssum %>% filter(d_exc_col_halo < 200)) +
      geom_point(mapping=aes(x=colony_size,
                             y=d_exc_col_halo),
                 #alpha=.8,
                 pch=halo_pch) +
      geom_polygon(data=p13poly,
                   mapping = aes(x = x, y = exp(y)), alpha=.2) +
      geom_path(data = p13,
                mapping = aes(x = predx,
                              y = exp(fit)),
                alpha=.8, lwd=1) +
      scale_y_continuous(trans='log', breaks=c(2, 6, 25, 90, 300, 900))+
      ylab(expression("log Density (kg " *km^-2* ")")) +
      labs(title = expression(paste("Excreta Density colony"^"-1","day"^"-1")),
           subtitle = subt13) +
      scale_x_continuous(trans='log', breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
      xlab('log Colony pairs x 100') +
      theme_minimal()+
      theme(plot.title = element_text(size = cex_t),
            plot.subtitle = element_markdown(hjust=0, size=cex_subt))

    plot13 %>% print
  }

  ggarrange(plot11, plot12, plot13,
            ncol=3, nrow=1,
            labels = 'auto')
  #heights = c(.9, .9, .9, 1))
}

ggsave('fig_excreta_density.png', width=12, height=4, bg='white')


################################################################################
################################################################################
################################################################################
# KM summaries that combine trips

if(file.exists('km_summaries.rds')){
  load('km_summaries.rds')
}else{
  mr <- data.frame()
  (lf <- list.files('poop_tracks2/'))
  trips$colony_size %>% sort %>% table
  trips$maxdist_linear %>% max
  (kms <- data.frame(km=0:186))
  i=3
  # troubleshooting
  #(is <- grep('20120712-Diabas-03-4181306', lf))
  #lf[is]
  #(i = is[1])
  for(i in 1:length(lf)){
    message('Preparing KM summary for poop track ',i,' out of ', length(lf),' ...')
    (lfi <- paste0('poop_tracks2/',lf[i]))
    load(lfi)

    energetics %>% as.data.frame %>% head
    trackexp %>% as.data.frame %>% head

    tripi <-
      trackexp %>%
      mutate(km = round(km)) %>%
      select(km, secs = t, gut, excreta, ingesta)
    tripi
    tripkm <- full_join(kms, tripi, by='km')
    tripkm %>% head
    tripkm %>% filter(km == 1)
    tripsumm <-
      tripkm %>%
      group_by(km) %>%
      summarize(use = length(unique(secs[!is.na(secs)])),
                gut = sum(gut, na.rm=TRUE),
                ingesta = -1*sum(ingesta, na.rm=TRUE),
                excreta = sum(excreta, na.rm=TRUE)) %>%
      pivot_longer(cols=ingesta:excreta, names_to = 'direction', values_to = 'grams') %>%
      group_by(direction) %>%
      mutate(grams_cum = cumsum(grams)) %>%
      mutate(grams_tot = grams_cum[n()]) %>%
      mutate(grams_p = grams_cum / grams_tot) %>%
      ungroup() %>%
      mutate(grams_p = ifelse(direction == 'ingesta', grams_p*-1, grams_p))

    tripsumm %>% head
    tripsumm <- data.frame(energetics %>% select(site:n_trips), tripsumm)
    tripsumm %>% head
    #ggplot(tripsumm, aes(x=km, y=grams)) + geom_point()

    mr <- rbind(mr, tripsumm)
  }

  # Save results
  save(mr, file='km_summaries.rds')
}

################################################################################
################################################################################
################################################################################

load('km_summaries.rds')
mr %>% head
mr$colony_size %>% unique %>% sort

kmd %>% head
cs <- mr$colony_size %>% unique %>% sort
(csi <- cs[1])
drim <- data.frame()
for(csi in cs){
  message(csi)
  (mrcs <- mr %>% filter(colony_size == csi))
  kmi=15
  for(kmi in 10:max(mr$km)){
    message('---- rim ', kmi-10,' - ', kmi)
    (mri <-
        mr %>%
        filter(colony_size == csi) %>%
        filter(direction=='excreta') %>%
        group_by(colony_size) %>%
        summarize(g_outer = mean(grams_cum[km==kmi], na.rm=TRUE),
                  g_inner = mean(grams_cum[km==(kmi-10)], na.rm=TRUE)) %>%
        mutate(area_inner = pi*(kmi-10)^2,
               area_outer = pi*kmi^2,
               rim_mid = kmi-5) %>%
        mutate(area_rim = area_outer - area_inner) %>%
        mutate(d_trip = g_outer / area_outer, # g km2
               d_trip_rim = (g_outer - g_inner) / area_rim) %>% # g km2
        mutate(d_col = (d_trip * colony_size * 2)/1000, # kg km2
               d_col_rim = (d_trip_rim * colony_size * 2)/1000) %>% # kg km2
        as.data.frame)
    drim <- rbind(drim, mri)
  }
}

drim

options(scipen = 999)
plot_a <-
  ggplot(drim %>% mutate(colony_size = colony_size / 100),
         aes(x=rim_mid, y=d_col, group=colony_size, color=colony_size)) +
  geom_path() +
  viridis::scale_color_viridis(direction=-1, trans='log',
                               breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
  scale_x_continuous(trans='log', limits=c(5, 150), breaks=c(5, 10, 25, 50, 90, 150)) +
  scale_y_continuous(trans='log', limits=c(.003, 30),
                     breaks=c(0.003, .03, 0.3, 3, 30)) +
  #scale_x_continuous(trans='log', limits=c(5, 145), breaks=c(5, 15, 45, 90, 145)) +
  xlab('Distance from colony (km)') +
  ylab(expression("log Density (kg " *km^-2* ")")) +
  theme_light() +
  labs(title = 'Throughout foraging range',
       color='log Colony size') +
  theme(legend.position = "none")

options(scipen = 999)
plot_b <-
  ggplot(drim %>% mutate(colony_size = colony_size / 100),
         aes(x=rim_mid, y=d_col_rim, group=colony_size, color=colony_size)) +
  geom_path() +
  viridis::scale_color_viridis(direction=-1, trans='log',
                               breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
  scale_x_continuous(trans='log', limits=c(5, 150), breaks=c(5, 10, 25, 50, 90, 150)) +
  scale_y_continuous(trans='log', limits=c(.003, 30),
                     breaks=c(0.003, .03, 0.3, 3, 30)) +
  xlab('Distance from colony (km)') +
  ylab(expression("log Density (kg " *km^-2* ")")) +
  theme_light() +
  labs(title = '10-km running window',
       color='log Colony size')

ggarrange(plot_a, plot_b, nrow=1, labels='auto', widths=c(.85, 1))
ggsave('fig_excreta_km.png', width=12, height=4, bg='white')

################################################################################
################################################################################
################################################################################
# CDFs

(lf <- list.files('poop_tracks2/'))
i=1
cdfs <- data.frame()
for(i in 1:length(lf)){
  message(i,' out of ', length(lf),' ...')
  (lfi <- paste0('poop_tracks2/',lf[i]))
  load(lfi)
  energetics %>% as.data.frame %>% head
  trackexp %>% as.data.frame %>% head
  suppressMessages({
    cdfi <-
      trackexp %>%
      mutate(colony_size = energetics$colony_size[1]) %>%
      mutate(toting = sum(ingesta)) %>%
      mutate(totexc = sum(excreta)) %>%
      mutate(cuming = cumsum(ingesta)) %>%
      mutate(cumexc = cumsum(excreta)) %>%
      mutate(maxkm = max(km, na.rm=TRUE)) %>%
      mutate(tmax = max(t)) %>%
      mutate(tdur = max(t) - min(t)) %>%
      mutate(ping = cuming / toting) %>%
      mutate(pexc = cumexc / totexc) %>%
      mutate(pkm = km / maxkm) %>%
      mutate(pt = 1 - ((tmax - t) / tdur)) %>%
      mutate(ptr = round(pt, 2)) %>%
      group_by(colony_size, ptr) %>%
      summarize(pkm = mean(pkm),
                ing_active = diff(range(cuming))>0,
                exc_active = diff(range(cumexc))>0,
                ping = mean(ping),
                pexc = mean(pexc)) %>%
      as.data.frame #%>% head
  })
  cdfi$ptr
  cdfi$ing_active
  #ggplot(cdfi, aes(x=ptr, y=pkm)) + geom_point()
  #ggplot(cdfi %>% filter(ing_active == TRUE), aes(x=ptr, y=pkm)) + geom_point()
  cdfs <- rbind(cdfs, cdfi)
}

cdfs %>% head
cdfsumm <-
  cdfs %>%
  group_by(colony_size, ptr) %>%
  summarize(ing_active = length(which(ing_active)) / n(),
            exc_active = length(which(exc_active)) / n(),
            pkm = mean(pkm),
            ping = mean(ping),
            pexc = mean(pexc, na.rm=TRUE))
cdfsumm %>% head


options(scipen = 999)
ping <-
  ggplot() +
  geom_point(data =
              cdfsumm %>%
              mutate(colony_size = colony_size / 100),
            mapping = aes(x=ptr, y=ing_active, group=colony_size, color=colony_size),
            alpha=.6) +
  viridis::scale_color_viridis(direction=-1, trans='log',
                               breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
  xlab('Quantile of trip duration') +
  ylab('Fraction of trips with active feeding') +
  theme_light() +
  labs(title = 'Time distribution of ingestion',
       color='log Colony size') +
  theme(legend.position = "none")


pexc <-
  ggplot() +
  geom_point(data =
               cdfsumm %>%
               mutate(colony_size = colony_size / 100),
             mapping = aes(x=ptr, y=exc_active, group=colony_size, color=colony_size),
             alpha=.6) +
  viridis::scale_color_viridis(direction=-1, trans='log',
                               breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
  xlab('Quantile of trip duration') +
  ylab('Fraction of trips with active excretion') +
  theme_light() +
  labs(title = 'Time distribution of excretion',
       color='log Colony size') +
  theme(legend.position = "none")


options(scipen = 999)
pdist <-
  ggplot() +
  geom_path(data =
              cdfsumm %>%
              mutate(colony_size = colony_size / 100),
            mapping = aes(x=ptr, y=pkm, group=colony_size, color=colony_size),
            lwd=.5, alpha=.8) +
  viridis::scale_color_viridis(direction=-1, trans='log',
                               breaks=c(4, 10, 30, 100, 300, 1000, 4000)) +
  xlab('Quantile of trip duration') +
  ylab('Quantile distance from colony') +
  theme_light() +
  labs(title = 'Time distribution of distance',
       color='log Colony size')

ggarrange(ping, pexc, pdist, nrow=1, labels='auto', widths=c(.75, .75, 1))
ggsave('fig_quantiles.png', width=12, height=4, bg='white')


################################################################################
################################################################################
################################################################################
# Since foraging area increases with foraging range by a power of 2
# for each bird, the same amount of food will be dispersed over an increasingly larger area,
#
# how does nutrient cycling relate to the concept of a predation halo in central place foraging?
# what happens to nutrient cycling in ashmole's halo?
# self-fertilizing food webs

# the fact that nutrient fertilization involves thresholds


################################################################################
################################################################################
################################################################################

