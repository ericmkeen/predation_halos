################################################################################
################################################################################
# The nutrient dimension of predation halos
# Eric K. Ezell, ekezell@sewanee.edu
#
# Step 5 Sensitivity analysis
################################################################################
################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('halo_functions.R')
source('halo_data_prep.R')
source('halo_parameters.R')

(sites <-
    trips %>%
    group_by(site, colony_size) %>%
    tally %>%
    arrange(n, colony_size) %>%
    filter(n >= 5) %>%
    pull(site))

################################################################################
################################################################################
# Subsampling

# if(file.exists('tripsub.rds')){
#   message('Subsampled trips file exists! Loading now.')
#   load('tripsub.rds')
# }else{
#   message('Subsampling trips ...')
#   trips %>% head
#   tripsub <- sample_trips(trips)
#   tripsub %>% head
#   tripsub %>% nrow # 1066
#   # for colonies w n < 25, keep all
#   # for colonies w 25 > n < 100, keep 25
#   # for colonies w n > 100, keep 25%
#   save(tripsub, file='tripsub.rds')
# }

################################################################################
################################################################################
# Simulating feeding distribution - centered

source('halo_parameters.R')
sim <- list(arrange = 'random',
            pt_mean = .5,
            pt_sd = 10,
            plot = FALSE)

trip_model(trips_to_process = 1579:nrow(trips),
           params = params,
           trips = trips,
           tracks = tracks,
           sim = sim,
           tracks_dir = 'sim/poop_tracks_centered/',
           diagnostics_dir = NULL,
           maps_dir = NULL,
           redo_analysis = TRUE,
           plots = FALSE)

################################################################################
################################################################################
# Simulating feeding distribution - left

source('halo_parameters.R')
sim <- list(arrange = 'frontload',
            pt_mean = .2,
            pt_sd = .2,
            plot = FALSE)

trip_model(trips_to_process = 1637:nrow(trips),
           params = params,
           trips = trips,
           tracks = tracks,
           sim = sim,
           tracks_dir = 'sim/poop_tracks_left/',
           diagnostics_dir = NULL,
           maps_dir = NULL,
           redo_analysis = TRUE,
           plots = FALSE)

################################################################################
################################################################################
# Simulating feeding distribution - right

source('halo_parameters.R')
sim <- list(arrange = 'backload',
            pt_mean = .8,
            pt_sd = .2,
            plot = FALSE)

trip_model(trips_to_process = 1:nrow(trips),
           params = params,
           trips = trips,
           tracks = tracks,
           sim = sim,
           tracks_dir = 'sim/poop_tracks_right/',
           diagnostics_dir = NULL,
           maps_dir = NULL,
           redo_analysis = TRUE,
           plots = FALSE)


################################################################################
################################################################################
# Excretion rate

# Half =========================================================================

source('halo_parameters.R')
params$excretion_rate <- 0.5*params$excretion_rate

trip_model(trips_to_process = 1:nrow(trips),
           params = params,
           trips = trips,
           tracks = tracks,
           tracks_dir = 'sim/excretion_half/',
           diagnostics_dir = NULL,
           maps_dir = NULL,
           redo_analysis = TRUE,
           plots = FALSE)


# Twice ========================================================================

source('halo_parameters.R')
params$excretion_rate <- 2*params$excretion_rate

trip_model(trips_to_process = 1:nrow(trips),
           params = params,
           trips = trips,
           tracks = tracks,
           tracks_dir = 'sim/excretion_twice/',
           diagnostics_dir = NULL,
           maps_dir = NULL,
           redo_analysis = TRUE,
           plots = FALSE)

################################################################################
################################################################################
# Transit time

# Half =========================================================================

source('halo_parameters.R')
params$ingesta_transit_time <- 0.5*params$ingesta_transit_time

trip_model(trips_to_process = 1:nrow(trips),
           params = params,
           trips = trips,
           tracks = tracks,
           tracks_dir = 'sim/transit_half/',
           diagnostics_dir = NULL,
           maps_dir = NULL,
           redo_analysis = TRUE,
           plots = FALSE)


# Twice ========================================================================

source('halo_parameters.R')
params$ingesta_transit_time <- 2*params$ingesta_transit_time

trip_model(trips_to_process = 1:nrow(trips),
           params = params,
           trips = trips,
           tracks = tracks,
           tracks_dir = 'sim/transit_twice/',
           diagnostics_dir = NULL,
           maps_dir = NULL,
           redo_analysis = TRUE,
           plots = FALSE)


################################################################################
################################################################################
# Ingestion rate

# Half =========================================================================

source('halo_parameters.R')
params$max_ingestion_rate <- 0.5*params$max_ingestion_rate

trip_model(trips_to_process = 1:nrow(trips),
           params = params,
           trips = trips,
           tracks = tracks,
           tracks_dir = 'sim/ingestion_half/',
           diagnostics_dir = NULL,
           maps_dir = NULL,
           redo_analysis = TRUE,
           plots = FALSE)


# Twice ========================================================================

source('halo_parameters.R')
params$max_ingestion_rate <- 2*params$max_ingestion_rate

trip_model(trips_to_process = 1:nrow(trips),
           params = params,
           trips = trips,
           tracks = tracks,
           tracks_dir = 'sim/ingestion_twice/',
           diagnostics_dir = NULL,
           maps_dir = NULL,
           redo_analysis = TRUE,
           plots = FALSE)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Traxtabs

# Main results
trax <- make_trax()
trax <- trax %>% filter(site %in% sites)

# Sensitivity results ==========================================================
trax_ctr <- make_trax(tracks_dir = 'sim/poop_tracks_centered/',
                      result_suffix = '_ctr',
                      result_dir = 'sim/')
trax_ctr <- trax_ctr %>% filter(site %in% sites)

trax_left <- make_trax(tracks_dir = 'sim/poop_tracks_left/',
                       result_suffix = '_left',
                       result_dir = 'sim/')
trax_left <- trax_left %>% filter(site %in% sites)

trax_right <- make_trax(tracks_dir = 'sim/poop_tracks_right/',
                        result_suffix = '_right',
                        result_dir = 'sim/')
trax_right <- trax_right %>% filter(site %in% sites)


trax_transit_half <- make_trax(tracks_dir = 'sim/transit_half/',
                       result_suffix = '_transit_half',
                       result_dir = 'sim/')
trax_transit_half <- trax_transit_half %>% filter(site %in% sites)

trax_transit_twice <- make_trax(tracks_dir = 'sim/transit_twice/',
                               result_suffix = '_transit_twice',
                               result_dir = 'sim/')
trax_transit_twice <- trax_transit_twice %>% filter(site %in% sites)


trax_excretion_half <- make_trax(tracks_dir = 'sim/excretion_half/',
                               result_suffix = '_excretion_half',
                               result_dir = 'sim/')
trax_excretion_half <- trax_excretion_half %>% filter(site %in% sites)

trax_excretion_twice <- make_trax(tracks_dir = 'sim/excretion_twice/',
                                result_suffix = '_excretion_twice',
                                result_dir = 'sim/')
trax_excretion_twice <- trax_excretion_twice %>% filter(site %in% sites)


trax_ingestion_half <- make_trax(tracks_dir = 'sim/ingestion_half/',
                                 result_suffix = '_ingestion_half',
                                 result_dir = 'sim/')
trax_ingestion_half <- trax_ingestion_half %>% filter(site %in% sites)

trax_ingestion_twice <- make_trax(tracks_dir = 'sim/ingestion_twice/',
                                  result_suffix = '_ingestion_twice',
                                  result_dir = 'sim/')
trax_ingestion_twice <- trax_ingestion_twice %>% filter(site %in% sites)



traxtab <- trax_results(trax_ctr)
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
ggsave('sim_ctr.png', width=10, height=12, bg='white')

traxtab <- trax_results(trax_left)
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
ggsave('sim_left.png', width=10, height=12, bg='white')

traxtab <- trax_results(trax_right)
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
ggsave('sim_right.png', width=10, height=12, bg='white')

traxtab <- trax_results(trax_transit_half)
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
ggsave('sim_transit_half.png', width=10, height=12, bg='white')

traxtab <- trax_results(trax_transit_twice)
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
ggsave('sim_transit_twice.png', width=10, height=12, bg='white')

traxtab <- trax_results(trax_excretion_half)
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
ggsave('sim_excretion_half.png', width=10, height=12, bg='white')

traxtab <- trax_results(trax_excretion_twice)
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
ggsave('sim_excretion_twice.png', width=10, height=12, bg='white')

traxtab <- trax_results(trax_ingestion_half)
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
ggsave('sim_ingestion_half.png', width=10, height=12, bg='white')

traxtab <- trax_results(trax_ingestion_twice)
trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
ggsave('sim_ingestion_twice.png', width=10, height=12, bg='white')

# # Foraging range
# traxtab <- trax_results(trax, col_lm = c(0.24178, 0.81420))
# trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
#
# traxtab <- trax_results(trax, col_lm = c(.5*0.24178, 1.5*0.81420))
# trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)
#
# traxtab <- trax_results(trax, col_lm = c(2*0.24178, -1.5*0.81420))
# trax_plot(traxtab, img_xy = c(12, 64), img_size=.65)

################################################################################
################################################################################
################################################################################
# Summarize

results_table <-
  rbind(
    trax_results(trax) %>% mutate(Scenario = 'Original predictions'),
    trax_results(trax_ctr) %>% mutate(Scenario = 'Simulated: centered feeding'),
    trax_results(trax_left) %>% mutate(Scenario = 'Simulated: early feeding'),
    trax_results(trax_right) %>% mutate(Scenario = 'Simulated: late feeding'),
    trax_results(trax_transit_half) %>% mutate(Scenario = 'Simulated: 0.5x transit time'),
    trax_results(trax_transit_twice) %>% mutate(Scenario = 'Simulated: 2.0x transit time'),
    trax_results(trax_excretion_half) %>% mutate(Scenario = 'Simulated: 0.5x excretion rate'),
    trax_results(trax_excretion_twice) %>% mutate(Scenario = 'Simulated: 2.0x excretion rate'),
    trax_results(trax_ingestion_half) %>% mutate(Scenario = 'Simulated: 0.5x ingestion rate'),
    trax_results(trax_ingestion_twice) %>% mutate(Scenario = 'Simulated: 2.0x ingestion rate')
  )

################################################################################
# Share exported to colony

results_table %>% head

rt <-
  results_table %>%
  group_by(colony_size) %>%
  mutate(frac_obs = col_exc_mean[Scenario=='Original predictions'][1]) %>%
  ungroup() %>%
  group_by(colony_size, Scenario) %>%
  summarize(frac = col_exc_mean,
            effect_frac = 100*(col_exc_mean / frac_obs) - 100) %>%
  ungroup %>%
  #as.data.frame
  group_by(Scenario) %>%
  summarize(frac_mn = mean(frac),
            frac_sd = sd(frac),
            effect_frac_mn = mean(effect_frac),
            effect_frac_sd = sd(effect_frac)) %>%
  mutate(`Share Exported to Colony` = paste0(round(frac_mn, 2), ' (', round(frac_sd, 2),')')) %>%
  mutate(`Effect` = paste0(round(effect_frac_mn, 1), '% (', round(effect_frac_sd, 1),'%)')) %>%
  select(Scenario, `Share Exported to Colony`, `Effect`) %>%
  slice(1, 8, 9, 10, 4, 7, 2, 5, 3, 6)

rt[1, 3] <- '-'
#rt <- rt[-1, ]

rt %>%
  tt() %>%
  style_tt(j = 1, align = "r") %>%
  style_tt(j = 2:3, align = "c") %>%
  style_tt(i = c(3, 8),
           j = 2:3,
           color = "red") %>%
  style_tt(i = 1,
           j = 2:3,
           color = "black") %>%
  style_tt(i = c(3, 4, 7:8),
           bold = TRUE)

################################################################################
# Guano allocation (halo v rim)

results_table %>% head

rt <-
  results_table %>%
  group_by(colony_size) %>%
  mutate(q50_obs = q50[Scenario=='Original predictions'][1],
         inner_obs = inner_sea_mean[Scenario=='Original predictions'][1],
         outer_obs = outer_sea_mean[Scenario=='Original predictions'][1],
         dinner_obs = d_inner[Scenario == 'Original predictions'][1],
         douter_obs = d_outer[Scenario == 'Original predictions'][1],
         dinner_col_obs = d_inner_col[Scenario == 'Original predictions'][1],
         douter_col_obs = d_outer_col[Scenario == 'Original predictions'][1]
  ) %>%
  ungroup() %>%
  group_by(colony_size, Scenario) %>%
  summarize(q50_eff = 100*(q50 / q50_obs) - 100,
            inner_eff = 100*(inner_sea_mean / inner_obs) - 100,
            outer_eff = 100*(outer_sea_mean / outer_obs) - 100,
            dinner_eff = 100*(d_inner / dinner_obs) - 100,
            douter_eff = 100*(d_outer / douter_obs) - 100,
            dinner_col_eff = 100*(d_inner_col / dinner_col_obs) - 100,
            douter_col_eff = 100*(d_outer_col / douter_col_obs) - 100) %>%
  ungroup %>%
  #as.data.frame %>% head
  group_by(Scenario) %>%
  summarize(q50_eff_mn = mean(q50_eff),
            q50_eff_sd = sd(q50_eff),
            inner_eff_mn = mean(inner_eff),
            inner_eff_sd = sd(inner_eff),
            outer_eff_mn = mean(outer_eff),
            outer_eff_sd = sd(outer_eff),
            dinner_eff_mn = mean(dinner_eff),
            dinner_eff_sd = sd(dinner_eff),
            douter_eff_mn = mean(douter_eff),
            douter_eff_sd = sd(douter_eff),
            dinner_eff_col_mn = mean(dinner_col_eff),
            dinner_eff_col_sd = sd(dinner_col_eff),
            douter_eff_col_mn = mean(douter_col_eff),
            douter_eff_col_sd = sd(douter_col_eff)) %>%
    mutate(`Excretion<br>midpoint` = paste0(round(q50_eff_mn, 1), '%<br>(', round(q50_eff_sd, 1),'%)')) %>%
    mutate(`Share<br>Inner<br>Halo` = paste0(round(inner_eff_mn, 1), '%<br>(', round(inner_eff_sd, 1),'%)')) %>%
    mutate(`Share<br>Outer<br>Rim` = paste0(round(outer_eff_mn, 1), '%<br>(', round(outer_eff_sd, 1),'%)')) %>%
  mutate(`Density<br>Inner<br>Halo<br>(per bird)` = paste0(round(dinner_eff_mn, 1), '%<br>(', round(dinner_eff_sd, 1),'%)')) %>%
  mutate(`Density<br>Outer<br>Rim<br>(per bird)` = paste0(round(douter_eff_mn, 1), '%<br>(', round(douter_eff_sd, 1),'%)')) %>%
  mutate(`Density<br>Inner<br>Halo<br>(colony)` = paste0(round(dinner_eff_col_mn, 1), '%<br>(', round(dinner_eff_col_sd, 1),'%)')) %>%
  mutate(`Density<br>Outer<br>Rim<br>(colony)` = paste0(round(douter_eff_col_mn, 1), '%<br>(', round(douter_eff_col_sd, 1),'%)')) %>%
  select(Scenario,
         `Excretion<br>midpoint`,
         `Share<br>Inner<br>Halo`,
         `Share<br>Outer<br>Rim`,
         `Density<br>Inner<br>Halo<br>(per bird)`,
         `Density<br>Outer<br>Rim<br>(per bird)`,
         `Density<br>Inner<br>Halo<br>(colony)`,
         `Density<br>Outer<br>Rim<br>(colony)`) %>%
    slice(1, 8, 9, 10, 4, 7, 2, 5, 3, 6) %>%
  mutate(Scenario = gsub('ted: ', 'ted:<br>', Scenario))

#rt[1, 2:4] <- '-'
rt <- rt[-1, ]

rt %>%
  tt() %>%
  style_tt(j = 1, align = "r") %>%
  style_tt(j = 2:3, align = "c") %>%
  style_tt(i = c(1:4, 7, 9),
           j = 2,
           color = "red") %>%
  style_tt(i = c(1, 2, 6, 9),
           j = 3,
           color = "red") %>%
  style_tt(i = c(3, 4, 7, 8),
           j = 4,
           color = "red") %>%
  style_tt(i = c(5, 6),
           j = 5,
           color = "red") %>%
  style_tt(i = c(1:3, 6, 8),
           j = 6,
           color = "red") %>%
  style_tt(i = c(5, 6),
           j = 7,
           color = "red") %>%
  style_tt(i = c(1:3, 6, 8),
           j = 8,
           color = "red") %>%
  style_tt(i = c(1:3),
           j=2,
           bold = TRUE) %>%
    style_tt(i = c(2),
             j=3,
             bold = TRUE) %>%
    style_tt(i = c(2),
             j=3,
             bold = TRUE) %>%
    style_tt(i = c(1:3),
             j=4,
             bold = TRUE) %>%
    style_tt(i = c(1:3, 6:7),
             j=5,
             bold = TRUE) %>%
  style_tt(i = c(1:3, 6:7),
           j=6,
           bold = TRUE) %>%
  style_tt(i = c(1:3, 6),
           j=7,
           bold = TRUE) %>%
  style_tt(i = c(1:3, 6:7),
           j=8,
           bold = TRUE)

################################################################################





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Re-distribution

# Main results
halo <- halo_summaries(tracks_dir = 'poop_tracks/',
                      result_dir = '',
                      result_suffix = '')
halo <- halo %>% filter(site %in% sites)
halos <- halo_summary_summary(halo)
plot_redistribution(halos)

# Sensitivity results ==========================================================

# KM summaries
halo_ctr <- halo_summaries(tracks_dir = 'sim/poop_tracks_centered/',
                       result_dir = 'sim/',
                       result_suffix = '_ctr')
names(halo_ctr) <- names(halo)
halo_ctr <- halo_ctr %>% filter(site %in% sites)

halo_left <- halo_summaries(tracks_dir = 'sim/poop_tracks_left/',
                           result_dir = 'sim/',
                           result_suffix = '_left')
halo_left <- halo_left %>% filter(site %in% sites)

halo_right <- halo_summaries(tracks_dir = 'sim/poop_tracks_right/',
                            result_dir = 'sim/',
                            result_suffix = '_right')
halo_right <- halo_right %>% filter(site %in% sites)

halo_transit_half <- halo_summaries(tracks_dir = 'sim/transit_half/',
                            result_dir = 'sim/',
                            result_suffix = '_transit_half')
halo_transit_half <- halo_transit_half %>% filter(site %in% sites)

halo_transit_twice <- halo_summaries(tracks_dir = 'sim/transit_twice/',
                                    result_dir = 'sim/',
                                    result_suffix = '_transit_twice')
halo_transit_twice <- halo_transit_twice %>% filter(site %in% sites)

halo_excretion_half <- halo_summaries(tracks_dir = 'sim/excretion_half/',
                                    result_dir = 'sim/',
                                    result_suffix = '_excretion_half')
halo_excretion_half <- halo_excretion_half %>% filter(site %in% sites)

halo_excretion_twice <- halo_summaries(tracks_dir = 'sim/excretion_twice/',
                                     result_dir = 'sim/',
                                     result_suffix = '_excretion_twice')
halo_excretion_twice <- halo_excretion_twice %>% filter(site %in% sites)

halo_ingestion_half <- halo_summaries(tracks_dir = 'sim/ingestion_half/',
                                    result_dir = 'sim/',
                                    result_suffix = '_ingestion_half')
halo_ingestion_half <- halo_ingestion_half %>% filter(site %in% sites)

halo_ingestion_twice <- halo_summaries(tracks_dir = 'sim/ingestion_twice/',
                                     result_dir = 'sim/',
                                     result_suffix = '_ingestion_twice')
halo_ingestion_twice <- halo_ingestion_twice %>% filter(site %in% sites)

# Summaries of the summaries
halos_ctr <- halo_summary_summary(halo_ctr)
halos_left <- halo_summary_summary(halo_left)
halos_right <- halo_summary_summary(halo_right)
halos_transit_half <- halo_summary_summary(halo_transit_half)
halos_transit_twice <- halo_summary_summary(halo_transit_twice)
halos_excretion_half <- halo_summary_summary(halo_excretion_half)
halos_excretion_twice <- halo_summary_summary(halo_excretion_twice)
halos_ingestion_half <- halo_summary_summary(halo_ingestion_half)
halos_ingestion_twice <- halo_summary_summary(halo_ingestion_twice)


# Plot  ========================================================================

subtitles <- c('Halo rims for <span style="color:darkblue;">**feeding**</span> and <span style="color:darkorange4;">**excretion**</span>',
               'Excretion - ingestion offset',
               'Proportional offset')

#==============================================

pobs <- plot_redistribution(halos,
                            titles = c('Original predictions', '', ''),
                            subtitles = subtitles, title_bold = TRUE,
                            labels = c('a', 'b', 'c'),
                            ylims3=c(-.2, .75))

#==============================================

ggarrange(
  pobs,
  plot_redistribution(halos_ctr,
                      titles = c('Simulated: centered feeding', '', ''),
                      subtitles = subtitles, title_bold = TRUE,
                      labels = c('d', 'e', 'f'),
                      ylims3=c(-.2,  .75)),
  plot_redistribution(halos_left,
                      titles = c('Simulated: early feeding', '', ''),
                      subtitles = subtitles, title_bold = TRUE,
                      labels = c('g', 'h', 'i'),
                      ylims3=c(-.2,  .75)),
  plot_redistribution(halos_right,
                      titles = c('Simulated: late feeding', '', ''),
                      subtitles = subtitles, title_bold = TRUE,
                      labels = c('j', 'k', 'l'),
                      ylims3=c(-.2,  .75)),
  nrow=4)
ggsave('sensitivity-halo-feeding-time-tall.png', height=11, width=8)
ggsave('sensitivity-halo-feeding-time-wide.png', height=10, width=13)

#==============================================

ggarrange(
  pobs,
  plot_redistribution(halos_transit_half,
                      titles = c('Simulated: 0.5x transit time', '', ''),
                      subtitles = subtitles, title_bold = TRUE,
                      labels = c('d', 'e', 'f'),
                      ylims3=c(-.2,  .75)),
  plot_redistribution(halos_transit_twice,
                      titles = c('Simulated: 2.0x transit time', '', ''),
                      subtitles = subtitles, title_bold = TRUE,
                      labels = c('g', 'h', 'i'),
                      ylims3=c(-.2,  .75)),
  nrow=3)
ggsave('sensitivity-halo-transit-time.png', height=8, width=10)

#==============================================

ggarrange(
  pobs,
  plot_redistribution(halos_excretion_half,
                      titles = c('Simulated: 0.5x excretion rate', '', ''),
                      subtitles = subtitles, title_bold = TRUE,
                      labels = c('d', 'e', 'f'),
                      ylims3=c(-.2,  .75)),
  plot_redistribution(halos_excretion_twice,
                      titles = c('Simulated: 2.0x excretion rate', '', ''),
                      subtitles = subtitles, title_bold = TRUE,
                      labels = c('g', 'h', 'i'),
                      ylims3=c(-.2,  .75)),
  nrow=3)
ggsave('sensitivity-halo-excretion-rate.png', height=8, width=10)

#==============================================

ggarrange(
  pobs,
  plot_redistribution(halos_ingestion_half,
                      titles = c('Simulated: 0.5x ingestion rate', '', ''),
                      subtitles = subtitles, title_bold = TRUE,
                      labels = c('d', 'e', 'f'),
                      ylims3=c(-.2, 1)),
  plot_redistribution(halos_ingestion_twice,
                      titles = c('Simulated: 2.0x ingestion rate', '', ''),
                      subtitles = subtitles, title_bold = TRUE,
                      labels = c('g', 'h', 'i'),
                      ylims3=c(-.2, 1)),
  nrow=3)
ggsave('sensitivity-halo-ingestion-rate.png', height=8, width=10)


# Summarize ====================================================================

results_table <-
  rbind(halos %>% mutate(Scenario = 'Original predictions'),

        halos_ctr %>% mutate(Scenario = 'Simulated: centered feeding'),
        halos_left %>% mutate(Scenario = 'Simulated: early feeding'),
        halos_right %>% mutate(Scenario = 'Simulated: late feeding'),

        halos_transit_half %>% mutate(Scenario = 'Simulated: 0.5x transit time'),
        halos_transit_twice %>% mutate(Scenario = 'Simulated: 2.0x transit time'),

        halos_excretion_half %>% mutate(Scenario = 'Simulated: 0.5x excretion rate'),
        halos_excretion_twice %>% mutate(Scenario = 'Simulated: 2.0x excretion rate'),

        halos_ingestion_half %>% mutate(Scenario = 'Simulated: 0.5x ingestion rate'),
        halos_ingestion_twice %>% mutate(Scenario = 'Simulated: 2.0x ingestion rate')
  ) %>%
  group_by(Scenario) %>%
  summarize(diff_mn = mean(exc_ing, na.rm=TRUE),
            diff_sd = sd(exc_ing, na.rm=TRUE),
            pdiff_mn = mean(exc_ing_prop, na.rm=TRUE),
            pdiff_sd = sd(exc_ing_prop, na.rm=TRUE)) %>%
  mutate(pdiff_obs = pdiff_mn[1]) %>%
  mutate(effect = 100*(pdiff_mn / pdiff_obs) - 100) %>%
  mutate(`Excretion - Ingestion Offset` = paste0(round(diff_mn, 3), ' (', round(diff_sd, 3),')')) %>%
  mutate(`Proportional Offset` = paste0(round(pdiff_mn, 3), ' (', round(pdiff_sd, 3),')')) %>%
  mutate(`Effect on Offset` = paste0(round(effect, 1), '%')) %>%
  select(Scenario, `Excretion - Ingestion Offset`, `Proportional Offset`, `Effect on Offset`) %>%
  slice(1, 8, 9, 10, 4, 7, 2, 5, 3, 6)

 results_table$`Effect on Offset`[1] <- '-'

 results_table %>%
   tt() %>%
  style_tt(j = 1, align = "r") %>%
  style_tt(j = 2:4, align = "c") %>%
  style_tt(
     i = which(results_table$`Effect on Offset` < 0),
     j = 4,
     color = "red"
   ) %>%
   style_tt(
     i = 1,
     j = 4,
     color = "black"
   ) %>%
   style_tt(i = c(2:4, 8),
            bold = TRUE)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Quantiles

# Main results
cdfs <- trip_quantiles(tracks_dir = 'poop_tracks/',
                       result_dir = '',
                       result_suffix = '')

# Sensitivity results ==========================================================

cdfs_ctr <- trip_quantiles(tracks_dir = 'sim/poop_tracks_centered/',
                       result_dir = 'sim/',
                       result_suffix = '_ctr')
cdfs_left <- trip_quantiles(tracks_dir = 'sim/poop_tracks_left/',
                           result_dir = 'sim/',
                           result_suffix = '_left')
cdfs_right <- trip_quantiles(tracks_dir = 'sim/poop_tracks_right/',
                           result_dir = 'sim/',
                           result_suffix = '_right')
cdfs_transit_half <- trip_quantiles(tracks_dir = 'sim/transit_half/',
                           result_dir = 'sim/',
                           result_suffix = '_transit_half')
cdfs_transit_twice <- trip_quantiles(tracks_dir = 'sim/transit_twice/',
                                    result_dir = 'sim/',
                                    result_suffix = '_transit_twice')
cdfs_excretion_half <- trip_quantiles(tracks_dir = 'sim/excretion_half/',
                                    result_dir = 'sim/',
                                    result_suffix = '_excretion_half')
cdfs_excretion_twice <- trip_quantiles(tracks_dir = 'sim/excretion_twice/',
                                     result_dir = 'sim/',
                                     result_suffix = '_excretion_twice')
cdfs_ingestion_half <- trip_quantiles(tracks_dir = 'sim/ingestion_half/',
                                    result_dir = 'sim/',
                                    result_suffix = '_ingestion_half')
cdfs_ingestion_twice <- trip_quantiles(tracks_dir = 'sim/ingestion_twice/',
                                     result_dir = 'sim/',
                                     result_suffix = '_ingestion_twice')


# Plot  ========================================================================

quantile_lims = c(0.5, 1)
qc <- quantile_plot(cdfs, quantile_lims,
                    titles = c('Original predictions', '', ''),
                    labels=c('a', 'b', 'c'))
qc$plot

#==============================================
qc_ctr <- quantile_plot(cdfs_ctr, quantile_lims,
                        titles = c('Simulated: centered feeding', '', ''),
                        labels=c('d', 'e', 'f'))
qc_ctr$plot
qc_left <- quantile_plot(cdfs_left, quantile_lims,
                         titles = c('Simulated: early feeding', '', ''),
                         labels=c('g', 'h', 'i'))
qc_left$plot
qc_right <- quantile_plot(cdfs_right, quantile_lims,
                          titles = c('Simulated: late feeding', '', ''),
                          labels=c('j', 'k', 'l'))
qc_right$plot

ggarrange(qc$plot,
          qc_ctr$plot,
          qc_left$plot,
          qc_right$plot,
          nrow=4)
ggsave('sensitivity-quantile-feeding-time-tall.png', height=16, width=14)
ggsave('sensitivity-quantile-feeding-time-wide.png', height=12, width=18)

#==============================================
qc_transit_half <- quantile_plot(cdfs_transit_half, quantile_lims,
                                 titles = c('Simulated: 0.5x transit time', '', ''),
                                 labels=c('d', 'e', 'f'))
qc_transit_half$plot

qc_transit_twice <- quantile_plot(cdfs_transit_twice, quantile_lims,
                                  titles = c('Simulated: 2.0x transit time', '', ''),
                                  labels=c('g', 'h', 'i'))
qc_transit_twice$plot

ggarrange(qc$plot,
          qc_transit_half$plot,
          qc_transit_twice$plot,
          nrow=3)
ggsave('sensitivity-quantile-transit-time.png', height=10, width=14)

#==============================================
qc_excretion_half <- quantile_plot(cdfs_excretion_half, quantile_lims,
                                   titles = c('Simulated: 0.5x excretion rate', '', ''),
                                   labels=c('d', 'e', 'f'))
qc_excretion_half$plot
qc_excretion_twice <- quantile_plot(cdfs_excretion_twice, quantile_lims,
                                    titles = c('Simulated: 2.0x excretion rate', '', ''),
                                    labels=c('g', 'h', 'i'))
qc_excretion_twice$plot

ggarrange(qc$plot,
          qc_excretion_half$plot,
          qc_excretion_twice$plot,
          nrow=3)
ggsave('sensitivity-quantile-excretion-rate.png', height=10, width=14)

#==============================================
qc_ingestion_half <- quantile_plot(cdfs_ingestion_half, quantile_lims,
                                   titles = c('Simulated: 0.5x ingestion rate', '', ''),
                                   labels=c('d', 'e', 'f'))
qc_ingestion_half$plot
qc_ingestion_twice <- quantile_plot(cdfs_ingestion_twice, quantile_lims,
                                    titles = c('Simulated: 2.0x ingestion rate', '', ''),
                                    labels=c('g', 'h', 'i'))
qc_ingestion_twice$plot

ggarrange(qc$plot,
          qc_ingestion_half$plot,
          qc_ingestion_twice$plot,
          nrow=3)
ggsave('sensitivity-quantile-ingestion-rate.png', height=10, width=14)



# Summarize ====================================================================

results_table <-
  rbind(qc$ingexc %>% mutate(Scenario = 'Original predictions'),

      qc_ctr$ingexc %>% mutate(Scenario = 'Simulated: centered feeding'),
      qc_left$ingexc %>% mutate(Scenario = 'Simulated: early feeding'),
      qc_right$ingexc %>% mutate(Scenario = 'Simulated: late feeding'),

      qc_transit_half$ingexc %>% mutate(Scenario = 'Simulated: 0.5x transit time'),
      qc_transit_twice$ingexc %>% mutate(Scenario = 'Simulated: 2.0x transit time'),

      qc_excretion_half$ingexc %>% mutate(Scenario = 'Simulated: 0.5x excretion rate'),
      qc_excretion_twice$ingexc %>% mutate(Scenario = 'Simulated: 2.0x excretion rate'),

      qc_ingestion_half$ingexc %>% mutate(Scenario = 'Simulated: 0.5x ingestion rate'),
      qc_ingestion_twice$ingexc %>% mutate(Scenario = 'Simulated: 2.0x ingestion rate')
      ) %>%
  group_by(Scenario) %>%
  summarize(ing = mean(pkm_ing50, na.rm=TRUE),
            ing_sd = sd(pkm_ing50, na.rm=TRUE),
            exc = mean(pkm_exc50, na.rm=TRUE),
            exc_sd = sd(pkm_exc50, na.rm=TRUE),
            diff = mean(pdiff, na.rm=TRUE),
            diff_sd = sd(pdiff, na.rm=TRUE)) %>%
  mutate(`Ingestion midpoint` = paste0(round(ing, 3), ' (', round(ing_sd, 3),')')) %>%
  mutate(`Excretion midpoint` = paste0(round(exc, 3), ' (', round(exc_sd, 3),')')) %>%
  mutate(`Ing. - Exc. Offset` = paste0(round(diff, 3), ' (', round(diff_sd, 3),')')) %>%
  mutate(pdiff_obs = diff[1]) %>%
  mutate(`Effect on Offset` = paste0(round((100*(diff / pdiff_obs)) - 100, 1), '%')) %>%
  select(-(ing:diff_sd), -pdiff_obs) %>%
  slice(1, 8, 9, 10, 4, 7, 2, 5, 3, 6)

results_table$`Effect on Offset`[1] <- '-'

results_table %>% tt() %>%
  style_tt(j = 1, align = "r") %>%
  style_tt(j = 2:5, align = "c") %>%
  style_tt(
    i = which(results_table$`Effect on Offset` < 0),
    j = 5,
    color = "red"
  ) %>%
  style_tt(
    i = 1,
    j = 5,
    color = "black"
  ) %>%
  style_tt(i = c(2:4),
           bold = TRUE)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
