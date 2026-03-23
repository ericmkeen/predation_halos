################################################################################
################################################################################
# The nutrient dimension of predation halos
# Eric K. Ezell, ekezell@sewanee.edu
#
# Step 4 Results
################################################################################
################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('halo_functions.R')
source('halo_data_prep.R')

################################################################################
################################################################################

(sites <-
  trips %>%
  group_by(site, colony_size) %>%
  tally %>%
  arrange(n, colony_size) %>%
  filter(n >= 5) %>%
  pull(site))


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

trips %>% nrow
trips$duration %>% hist
trips$duration %>% mean
trips$duration %>% median
trips$duration %>% sd
trips$duration %>% min
trips$duration %>% max

trips$totdist %>% hist
trips$totdist %>% mean(na.rm=TRUE)
trips$totdist %>% median(na.rm=TRUE)
trips$totdist %>% sd(na.rm=TRUE)
trips$totdist %>% min(na.rm=TRUE)
trips$totdist %>% max(na.rm=TRUE)

################################################################################
################################################################################
################################################################################
# TRACKS

#load('tracks.rds')

tracks %>% names
tracks <-
  tracks %>%
  filter(site %in% sites) %>%
  mutate(id = paste(site, species, bird_id, deployment, tripID, sep='_')) %>%
  left_join(trip2join, by='site') %>%
  mutate(colony_size = round(colony_size/100))

tracks %>% head

if(FALSE){
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
}

################################################################################
################################################################################
# Splashdowns

(lf <- list.files('poop_tracks/'))
i=1
splashdowns <- data.frame()
for(i in 1:length(lf)){
  message(i,' out of ', length(lf),' ...')
  (lfi <- paste0('poop_tracks/',lf[i]))
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
# Halo map

options(expressions = 500000)
(pa <- halo_map(3701,
                xlims = c(-14.7, -13.4),
                ylims = c(64.25, 64.7),
                image_xy = c(-14.46, 64.64),
                image_size = .5,
                col_xy = c(-14.185, 64.595),
                col_size = 1,
                scale_width = 0.4,
                labels = c('(a)', '(b)', '(c)')))
(pb <- halo_map(227000,
                xlims = c(16.5, 22.75),
                ylims = c(73.6, 75.2),
                image_xy = NULL,
                col_xy = c(19.05, 74.42),
                col_size = 1,
                scale_width = 0.4,
                labels = c('(d)', '(e)', '(f)')))
(pc <- halo_map(472409,
                xlims = c(-54.05, -51.95),
                ylims = c(49, 50.1),
                image_xy = NULL,
                col_xy = c(-53.195, 49.76),
                col_size = 1,
                scale_width = 0.3,
                labels = c('(g)', '(h)', '(i)')))

ggarrange(pa, pb, pc, nrow=3)
ggsave('fig_halos.png',
       width=14, height=12,
       bg='white')


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# NEW TRAX PLOTS

trax <- make_trax()

(sites <-
    trips %>%
    group_by(site, colony_size) %>%
    tally %>%
    arrange(n, colony_size) %>%
    filter(n >= 5) %>%
    pull(site))

trax <- trax %>% filter(site %in% sites)
trax %>% nrow


# Foraging range summary =====================================================

trax %>% as.data.frame %>% head
trax$n_trips %>% table

traxkm <-
  trax %>%
  filter(site %in% sites) %>%
  group_by(colony_size, kmfloor, kmceiling) %>%
  summarize(n = length(unique(track_id)),
            ing_mean = mean(ingesta),
            exc_mean = mean(excreta),
            exc_sd = sd(excreta),
            exc_se = exc_sd / sqrt(n)) %>%
  arrange(colony_size) %>%
  mutate(colony_size = round(colony_size/200))

traxkm %>% head

colevels <- traxkm$colony_size %>% unique %>% sort %>% prettyNum(big.mark=',')
ggplot(traxkm,
       aes(x=kmfloor,
                   #y=colony_size,
                   y=factor(prettyNum(colony_size, big.mark=','), levels=colevels),
                   size=exc_mean, color=exc_mean)) +
  geom_point(shape='|') +
  viridis::scale_color_viridis(direction=-1,
                               option='rocket') +
  scale_size_continuous(range=c(1,15)) +
  theme_minimal() +
  labs(color = expression(paste("Mean g Guano trip"^"-1")),
       title = 'Guano throughout foraging range') +
  xlab('Distance from colony (km)') +
  ylab('log Colony pairs x 100') +
  guides(size = 'none')

ggsave('fig_guano_range.png', width=9, height=6, bg='white')


# misc stats
traxkm$exc_mean %>% mean
traxkm$exc_mean %>% sd
traxkm$exc_mean %>% range %>% round(4)

# ingestion midpoint / max km
trax %>%
  group_by(colony_size) %>%
  summarize(ing = mean(km_ing50),
            kmax = mean(kmax)) %>%
  mutate(rat = ing/kmax) %>%
  summarize(mean(rat),
            sd(rat))

# excretion midpoint vs ingestion midpoint
trax %>%
  group_by(colony_size) %>%
  summarize(ing = mean(km_ing50),
            exc = mean(q50)) %>%
  mutate(rat = exc/ing) %>%
  summarize(mean(rat),
            sd(rat))

r <- 1:200
(halo <- (pi*((0.84*r)^2)) - (pi*1^2))
rim <- (pi*r^2) - halo

plot(halo, type='l', col='goldenrod')
lines(rim, col='firebrick')


################################################################################

traxtab <- trax_results(trax)
traxtab$n_trips
traxtab$col_exc_tot %>% range
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65,
          y_drop_col = c(0.07, .25, 1, 3, 9, 25))
ggsave('fig_punchline.png', width=10, height=12, bg='white')


# guano prod
traxtab %>% as.data.frame %>% head
traxtab$exc_tot_mean %>% head
traxtab$exc_tot_mean %>% tail
traxtab$exc_tot_mean %>% mean
traxtab$exc_tot_mean %>% sd

trax %>%
  group_by(colony_size) %>%
  summarize(ing = mean(km_ing50),
            exc = mean(q50)) %>%
  mutate(rat = exc/ing) %>%
  summarize(mean(rat),
            sd(rat))

traxtab$outer_sea_mean %>% mean
traxtab$outer_sea_mean %>% sd

traxtab %>% names

# per trip =====================================================================
(smallest <- traxtab %>%
  arrange(colony_size) %>%
  head(8) %>%
  summarize(in_mn = mean(d_inner),
            in_sd = sd(d_inner),
            out_mn = mean(d_outer),
            out_sd = sd(d_outer),
            all_mn = mean(d_all),
            all_sd = sd(d_all)))

# largest
(largest <- traxtab %>%
  arrange(colony_size) %>%
  tail(8) %>%
    summarize(in_mn = mean(d_inner),
              in_sd = sd(d_inner),
              out_mn = mean(d_outer),
              out_sd = sd(d_outer),
              all_mn = mean(d_all),
              all_sd = sd(d_all)))

smallest$all_mn / largest$all_mn

# colony =======================================================================
(smallest <- traxtab %>%
   arrange(colony_size) %>%
   head(8) %>%
   summarize(in_mn = mean(d_inner_col),
             in_sd = sd(d_inner_col),
             out_mn = mean(d_outer_col),
             out_sd = sd(d_outer_col),
             all_mn = mean(d_all_col),
             all_sd = sd(d_all_col)))

# largest
(largest <- traxtab %>%
    arrange(colony_size) %>%
    tail(8) %>%
    summarize(in_mn = mean(d_inner_col),
              in_sd = sd(d_inner_col),
              out_mn = mean(d_outer_col),
              out_sd = sd(d_outer_col),
              all_mn = mean(d_all_col),
              all_sd = sd(d_all_col)))


traxtab$col_exc_mean %>% mean
traxtab$col_exc_mean %>% sd

traxtab$col_exc_mean %>% head(8) %>% mean
traxtab$col_exc_mean %>% head(8) %>% sd
traxtab$col_exc_mean %>% tail(8) %>% mean
traxtab$col_exc_mean %>% tail(8) %>% sd

traxtab$col_exc_tot %>% head(1)
traxtab$col_exc_tot_se %>% head(1)
traxtab$col_exc_tot %>% max
traxtab$col_exc_tot_se %>% tail(1)
traxtab$colony_size %>% range

traxtab %>% select(km_ing50_median, q50) %>% filter(q50 < km_ing50_median)
nrow(traxtab)
16/23
################################################################################
################################################################################
################################################################################
# Trip quantiles

cdfs <-
  trip_quantiles(tracks_dir = 'poop_tracks/',
                 result_dir = '',
                 result_suffix = '')

qcdfs <-
  quantile_plot(cdfs,
              quantile_lims = c(0.5, 1),
              titles = c(' ', ' ', ' '))
qcdfs$plot
ggsave('fig_quantiles_v2.png', width=13, height=4, bg='white')

qcdfs$ingexc$pdiff %>% mean

################################################################################
################################################################################
################################################################################

