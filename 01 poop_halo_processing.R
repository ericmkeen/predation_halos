# seabird poop halo

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(swfscMisc)
library(tidyr)

subsample <- FALSE
redo_analysis <- TRUE

################################################################################
################################################################################
################################################################################
# GPS track data from Allison Patterson et al 2022

#===============================================================================
# Trip summary key

(trips <- readr::read_delim('tripSummary.txt', delim=';'))
trips %>% as.data.frame %>% head
trips %>% nrow
trips_all <- trips

# Foraging range increases with colony size
ggplot(trips,
       aes(x=colony_size,
           y=maxdist_linear)) +
  scale_y_continuous(trans='log', breaks=c(1,2,5, 10, 20, 50, 100, 175)) +
  scale_x_continuous(trans='log', breaks=c(700, 2000, 5000, 15000, 50000, 150000, 500000)) +
  geom_jitter(alpha=.1, width=.08) +
  geom_smooth(method='lm') +
  theme_minimal() +
  xlab('log Colony Size') + ylab('log Max. Distance (km)')

# Proportion of area within 30km that is water
trips$prop_water_30km %>% hist

#===============================================================================
# GPS tracks

if(file.exists('tracks.rds')){
  load('tracks.rds')
}else{
  gps <- readr::read_delim('gpsTracks.txt', delim = ';')
  gps %>% as.data.frame %>% head
  gps$behaviour %>% table
  gps$band %>% table %>% sort

  # remove time at colony
  gps <- gps %>% filter(behaviour == 'inTrip')

  # Filter to trips
  gps <- gps %>% filter(tripID %in% trips$tripID,
                        bird_id %in% trips$bird_id)

  tracks <- gps %>%
    group_by(site, species, bird_id, year, tripID) %>%
    mutate(t = as.numeric(time)) %>%
    mutate(lon2 = lead(lon),
           lat2 = lead(lat),
           t2 = lead(t)) %>%
    mutate(secs = t2 - t) %>%
    rowwise() %>%
    mutate(meters = ifelse(any(is.na(c(lat, lon, lat2, lon2))),
                           NA,
                           1000*swfscMisc::distance(lat1=lat, lon1=lon, lat2=lat2, lon2=lon2, units='km'))) %>%
    ungroup() %>%
    mutate(ms = meters / secs) %>%
    mutate(kmh = ms*3.6) %>%
    mutate(kmh = ifelse(kmh > 120, NA, kmh)) %>%
    mutate(status = ifelse(kmh > 10, 'fly','sit')) %>%
    mutate(prev_status = lag(status),
           next_status = lead(status)) %>%
    as.data.frame

  tracks %>% head

  # Save dataset
  save(tracks, file='tracks.rds')
}
#
#
#
# exploratory plots  ===========================================================

#(tracks$t2 - tracks$t) %>% median(na.rm=TRUE)
tracks$kmh %>% hist
tracks %>% filter(kmh > 15) %>% pull(kmh) %>% mean
tracks %>% filter(kmh > 15) %>% pull(kmh) %>% median
tracks %>% filter(kmh > 15) %>% pull(kmh) %>% sd

ggplot(tracks, aes(x=kmh, fill=status)) + geom_histogram()

# test paths for birds with most data
(ids <- tracks$bird_id %>% table %>% sort %>% rev %>% head(8) %>% names)

ggplot(tracks %>%
         filter(band %in% ids),
       aes(x=lon, y=lat, color=status)) +
  geom_point() +
  facet_wrap(~bird_id, scales='free')

ggplot(tracks %>%
         filter(band %in% ids),
       aes(x=time, y=coldist, color=status)) +
  geom_point() +
  facet_wrap(~bird_id, scales='free_x')

  ggplot(trips, aes(x=colony_size, y=duration)) +
    geom_point(alpha=.2) +
    #(height=0, width=5000, alpha=.2) +
    scale_y_continuous(trans='log', n.breaks=10) +
    scale_x_continuous(trans='log') +
    geom_smooth(method='lm')

  trips %>% head
  ntrips <-
    trips %>%
    mutate(dt = lubridate::date(start)) %>%
    mutate(id = paste(site, species, bird_id, deployment, tripID, sep='_')) %>%
    group_by(colony_size, bird_id, dt) %>%
    tally()

  ggplot(ntrips, aes(x=colony_size, y=n)) +
    geom_jitter(height=.2, width=.05, alpha=.1) +
    scale_x_continuous(trans='log') +
    geom_smooth()


#===============================================================================
#===============================================================================
#===============================================================================
# set params
params <- list(excretion_rate = 0.00013, # % wet meal mass excreted per second (middle of the two values,  0.011 and 0.015, tried in Hilton)
               ingesta_transit_time = 70*60, # sec to first poop once foraging begins, hilton et al. 2000
               max_ingestion_rate = 0.037 , # grams per second, hilton pg 24 of 24
               max_feeding_time = 197*60, # minutes, mean of 186, 318, and 89, hilton et al pg 24 of pdf
               gut_capacity = 192,# g, or mL, hilton pg 24 of pdf
               murre_mass = 888, # g # from hilton et al. 2000
               FMR = 1789, # kJ, field metabolic rate
               RMR = 544, # kJ # basal energy needs
               Ebird = 1088, #kJ # daily energy needs of a resting bird
               Echick = 286, #kJ # daily energy needs of a chick
               Rfeed = 50, # J/s, energy during active feeding
               Rsit = 12.8, # J/s, energy during sitting on water
               E_density_diet = 6, # kJ/g; middle of 4-8 kJ from Hilton)
               trip_hours = 18, # based on the 18 hours per pair in Hilton
               dive_time = 77, # sec Evans et al 2013
               dives_per_bout = 7,
               post_dive_interval = 42, # s
               bout_duration = 833, # s 803 is observed from Evans et al; 833 makes sense based on params above
               bout_interval = 250, # s
               percent_feeding = 0.50, # percent of "sit" time spent feeding underwater; Hilton used 0.63; Evans et al. had 50%.
               assimilation_efficiency = 0.79, # hilton et al 2000b, used in Roth et al 2008
               excreta_ingesta_ratio = .427) # see below
params
# Proportion of ingesta that is excreted
# Hilton et al 2000b fed murres 10.88% of their body weight in fish
# Based on figure 3, the murre weight is 990g
# So they were fed 107g wet mass
# Using an equation from Larimer (1992), dry = 0.309*wet - 0.286, which is for sand lance,
# So this is 32.77g dry mass
# Based on figure 2, murres excreted 14g dry mass
# So they excreta was equal to 14/32.77 = 42.7% of their ingesta by dry weight


cost_of_flight_watts <- function(w){0.0018*(w^1.6)} # hilton et al 2000
# cost_of_flight_watts()
cost_of_flight_watts(888)
cost_of_flight_watts(888)


sit_session <- function(secs, toplot=FALSE){
  #secs <- 1*3600
  (bout_duration <- params$bout_duration)
  (dive <- params$dive_time)
  (total_feed_time <- round(secs * params$percent_feeding))
  (max_dives <- floor(total_feed_time / dive))

  if(secs < bout_duration){ bout_duration <- secs }
  bout_duration
  if(secs < dive | max_dives == 0){
    sit_seq <- rep('sit', times= secs)
  }else{

    # sequence for a single bout
    (post_dive_interval <- params$post_dive_interval)
    (dive_cycle <- dive + post_dive_interval)
    ((dives_per_bout <- floor(bout_duration / dive_cycle)))
    if(dives_per_bout > max_dives){dives_per_bout <- max_dives}
    (post_dive_interval <- floor((bout_duration - (dives_per_bout*dive)) / dives_per_bout))
    (dive_seq <- c(rep('feed', times = dive),
                   rep('sit', time = post_dive_interval)))
    (bout_seq <- rep(dive_seq, times= dives_per_bout))
    if(length(bout_seq) < bout_duration){
      bout_seq <- c(bout_seq,
                    rep('sit', times= (bout_duration - length(bout_seq))))
    }
    bout_seq
    length(bout_seq)
    bout_duration

    # sequence for a sit session
    (feed_time_per_bout <- length(which(bout_seq == 'feed')))
    (n_bouts <- floor(total_feed_time / feed_time_per_bout))
    #(n_bouts <- round(total_bout_time / bout_duration))
    (bout_interval <- floor((secs - (n_bouts * bout_duration)) / n_bouts))
    # sanity check
    bout_duration*n_bouts + bout_interval*n_bouts

    (bout_cycle_seq <- c(bout_seq,
                         rep('sit', times=bout_interval)))
    (sit_seq <- rep(bout_cycle_seq, times=n_bouts))

    if(length(sit_seq) < secs){
      sit_seq <- c(sit_seq,
                   rep('sit', times= (secs - length(sit_seq))))
    }
    length(sit_seq)
    secs
  }

  if(toplot){
    p <- ggplot2::ggplot(data.frame(t=1:length(sit_seq),
                                    sit_seq),
                         aes(x=t, y=sit_seq)) + geom_point(cex=.2)
    print(p)
  }
  return(sit_seq)
}

test <- sit_session(10000, toplot=TRUE)
test <- sit_session(7000, toplot=TRUE)
test <- sit_session(3600, toplot=TRUE)
test <- sit_session(900, toplot=TRUE)
test <- sit_session(700, toplot=TRUE)
test <- sit_session(500, toplot=TRUE)
test <- sit_session(200, toplot=TRUE)
test <- sit_session(154, toplot=TRUE)
test <- sit_session(100, toplot=TRUE)
test <- sit_session(78, toplot=TRUE)
test <- sit_session(50, toplot=TRUE)
test <- sit_session(3, toplot=TRUE)

#===============================================================================
#===============================================================================
#===============================================================================

# Loop settings for horizontal computing in Sewanee Ocean Lab
nrow(trips_all)
i_begin <- 1
i_end <- nrow(trips)
files_to_process <- i_begin:i_end
i=4303
(i=sample(1:5283, size=1))

for(i in files_to_process){
  (tripi <- trips[i, ]) %>% as.data.frame

  #if(tripi$duration < 48){ # only process tracks less than 2 days long
  if(TRUE){

    (tripid <- paste0(tripi$site,'_',tripi$year,'_',
                      tripi$bird_id,'_',tripi$tripID,
                      tripi$deployment))
    (fn <- paste0('poop_tracks2/',tripid,'.rds'))   # stage filename for result

    if(!file.exists(fn) | redo_analysis){
      message('======= Trip ', tripid,' (',i,' of ', nrow(trips),') =======')

      (tracki <- tracks %>% filter(bird_id == tripi$bird_id, tripID == tripi$tripID, deployment == tripi$deployment))
      tracki #%>% head
      nrow(tracki)

      # Handle anomalous file errors
      tracki$t2[nrow(tracki)] <- NA # Guarantee that final t2 is NA
      if(i == 1175){tracki <- tracki[1:214, ]}
      if(i == 1176){tracki <- tracki[1:11, ]}

      # Interpolate track
      message('--- interpolating track ...')
      trackexp <-
        tracki %>%
        filter(!is.na(t2)) %>%
        mutate(coldist2 = lead(coldist, 1)) %>%
        mutate(coldist2 = ifelse(is.na(coldist2), coldist, coldist2)) %>%
        arrange(t) %>%
        rowwise() %>%
        mutate(tnew = list(t:t2)) %>%
        mutate(x = list(seq(lon, lon2, length= (t2-t+1)))) %>%
        mutate(y = list(seq(lat, lat2, length= (t2-t+1)))) %>%
        mutate(km = list(seq(coldist, coldist2, length = (t2-t+1)))) %>%
        unnest(cols = c(tnew, x, y, km)) %>%
        ungroup() %>%
        select(-t, -t2, -lon, -lon2, -lat, -lat2, -dt, -dist, -secs, -time, -coldist, -coldist2, -meters) %>%
        rename(t = tnew) %>%
        select(-prev_status, -next_status) %>%
        mutate(q = dplyr::row_number() / dplyr::n()) %>%
        distinct(t, .keep_all=TRUE)

      trackexp %>% as.data.frame %>% head
      trackexp %>% nrow

      # fix NA status
      nas <- which(is.na(trackexp$status))
      if(length(nas)>0){
        trackexp$status[nas] <- 'error'
      }

      # splashdowns
      # https://www.int-res.com/articles/meps2013/475/m475p277.pdf
      # mead is 320m from colony, sd 126, 320 + 2 SDs = 572m
      # any sitting that happens within 0.572km of colony is considered a splashdown -- no foraging.
      message('--- flagging splashdown(s) ...')
      trackexp %>% as.data.frame %>% head
      #trackexp$km %>% plot
      trackexp
      (splash_zones <- which(trackexp$km < 0.572 & trackexp$q < 0.50))
      if(length(splash_zones)>0){
        trackexp$status[splash_zones] <- 'splashdown'
      }
      trackexp$status %>% table

      # Creating feeding sequence
      message('--- creating feeding sequence ...')
      trackexp$status
      if(FALSE){
        ggplot2::ggplot(trackexp, aes(x=t, y=status)) + geom_point(cex=.2)
      }
      trackexp1 <-
        trackexp %>%
        # Create group IDs for each status session
        mutate(status_last = lag(status, 1)) %>%
        mutate(status_last = ifelse(is.na(status_last), status, status_last)) %>%
        mutate(status_shift = ifelse(status_last != status,
                                     1, 0)) %>%
        mutate(status_id = cumsum(status_shift) + 1) %>%
        # group by these and work with each separately
        group_by(status_id) %>%
        group_split()

      length(trackexp1)
      trackexp2 <- data.frame()
      j=1
      for(j in 1:length(trackexp1)){
        (trackexpi <- trackexp1[[j]])
        stati <- trackexpi$status
        #message('--- --- status sequence ', j,' (',stati[1],', nrow=',nrow(trackexpi),')')
        trackexpi$status_detail <- stati
        if(stati[1] == 'sit'){
          trackexpi$t
          (secs <- trackexpi$t[nrow(trackexpi)] - trackexpi$t[1] + 1)
          trackexpi$status_detail <- sit_session(secs)
        }
        trackexp2 <- rbind(trackexp2, trackexpi)
      }

      trackexp2$status_detail %>% table
      if(FALSE){
        ggplot2::ggplot(trackexp2, aes(x=t, y=status_detail)) + geom_point(cex=.2)
      }

      # metabolic rate as function of temperature
      #63% of time sitting on water
      #MR = 17.39 - 0.60  x water temperature (Croll and McLaren, 1993)
      #(MR = 17.39 - (0.6*-2))
      #FMR + MR
      #(MR = 17.39 - (0.6*15))
      #1789 + MR
      #FMR = 1789 kJ / day # Cairns et al 1990

      # Get starting point energy needs
      message('--- calculating energy needs ...')

      # Calculate metabolic rate based on sea temperature
      (FMR <- params$FMR)
      percent_in_water <- 0.63 # Croll et al 1990
      (FMRwater = FMR*percent_in_water)
      (FMRdry = FMR*(1 - percent_in_water))
      (MR <- 17.39 - 0.6*tripi$sst)
      (MRnf = 17.39 - (0.6*5.67)) # croll et al. 1990
      (MRratio = MR/MRnf)
      (FMRwater = FMRwater*MRratio)
      (FMR = FMRwater + FMRdry)

      (duration <- tripi$duration)
      (n_trips <- ceiling(params$trip_hours / (duration*2)))
      (fly_per_trip <- trackexp2 %>% filter(status_detail == 'fly') %>% nrow)
      (sit_per_trip <- trackexp2 %>% filter(status_detail %in% c('sit', 'splashdown')) %>% nrow)
      (feed_per_trip <- trackexp2 %>% filter(status_detail == 'feed') %>% nrow)
      (error_per_trip <- trackexp2 %>% filter(status_detail == 'error') %>% nrow)
      (flap_cost_per_sec <- cost_of_flight_watts(params$murre_mass)/1000)  # kJ
      (tot_flight_costs <-  flap_cost_per_sec * fly_per_trip * n_trips) # kJ
      (rest_per_trip <- sit_per_trip - feed_per_trip) # sec
      (Ebird <- FMR)
      #Echick <- params$Echick / params$assimilation_efficiency)
      if(duration > 24){
        Ebird <- Ebird * (duration/24)
        #Echick <- Echick * (duration/24)
      }
      Ebird ; Echick
      (kJ_per_day <- Ebird + tot_flight_costs) #+
        # Echick + # built into FMR estimate
        # tot_feed_costs + tot_rest_costs + tot_error_costs) # kJ
      (food_per_day <- (kJ_per_day/ params$assimilation_efficiency) / params$E_density_diet) # grams
      (food_per_trip <- food_per_day / n_trips) # grams
      (ingest_rate <- food_per_trip / feed_per_trip) # grams/sec
      params$max_ingestion_rate
      # is this higher than the max possible? If so, replace with max possible.
      if(ingest_rate > params$max_ingestion_rate){ ingest_rate <- params$max_ingestion_rate }
      ingest_rate
      (energetics <- data.frame(tripi,
                               FMR, MR, MRnf, MRratio,
                               FMRwater, FMRdry,
                               AE = params$assimilation_efficiency,
                               Ebird,
                               #Echick,
                               n_trips,
                               fly_per_trip,
                               sit_per_trip,
                               feed_per_trip,
                               error_per_trip,
                               tot_flight_costs,
                               kJ_per_day,
                               food_per_day,
                               food_per_trip,
                               ingest_rate,
                               excreta_ingesta_ratio = params$excreta_ingesta_ratio))

      # figure out how much would be eaten if using the physiological maximum ingestion rate
      # convenience function ===================================================
      inout_model <- function(trackexp2,
                              params,
                              ingestion_scaling = 1){
        (ingesta_transit_time <- params$ingesta_transit_time)
        (gut_capacity <- params$gut_capacity)
        (excretion_rate <- params$excretion_rate)
        (max_ingestion_rate <- params$max_ingestion_rate * ingestion_scaling)
        #message(max_ingestion_rate)

        # dplyr attempt
        if(FALSE){
          trackexpi <-
            trackexp2 %>%
            mutate(transit_clock = 0,
                   gut = 0,
                   ingesta = 0,
                   excreta = 0) %>%
            rowwise() %>%
            mutate(pre_transit = lag(transit_clock, 1),
                   pre_gut = lag(gut, 1)) %>%
            mutate(pre_transit = ifelse(is.na(pre_transit), 0, pre_transit),
                   pre_gut = ifelse(is.na(pre_gut), 0, pre_gut)) %>%
            mutate(transit_clock = ifelse(pre_gut == 0, 0, transit_clock + 1)) %>%
            mutate(excreta = ifelse(transit_clock < ingesta_transit_time,
                                    0,
                                    pre_gut*excretion_rate)) %>%
            mutate(ingesta = ifelse(status_detail == 'feed',
                                    ifelse(pre_gut < gut_capacity,
                                           max_ingestion_rate,
                                           excreta),
                                    0)) %>%
            mutate(ingesta = ifelse(ingesta > (gut_capacity - pre_gut),
                                    (gut_capacity - pre_gut),
                                    ingesta)) %>%
            mutate(gut = pre_gut + ingesta - excreta) %>%
            ungroup()
          trackexpi %>% as.data.frame %>% head
        }

        if(TRUE){ # the loop way
          trackexpi <- trackexp2
          (stats <- trackexp2$status_detail)
          transit_clock <- gut <- ingesta <- excreta <- c(0)
          j=401
          #for(j in 1:400){
          for(j in 1:length(stats)){
            (this_status <- stats[j])
            #message(this_status)
            (this_transit_clock <- ifelse(tail(gut, 1) == 0,
                                          0,
                                          tail(transit_clock, 1) + 1))
            (this_excreta <- ifelse(this_transit_clock < ingesta_transit_time,
                                    0,
                                    tail(gut, 1)*excretion_rate))
            this_ingesta <- 0
            if(this_status == 'feed'){
              (this_ingesta <- ifelse(tail(gut, 1) < gut_capacity,
                                      min(c(max_ingestion_rate, (gut_capacity - tail(gut)[1]))),
                                      this_excreta))
            }
            (this_gut <- round(tail(gut, 1) - this_excreta + this_ingesta, 5))

            transit_clock[j] <- this_transit_clock
            gut[j] <- this_gut
            ingesta[j] <- this_ingesta
            excreta[j] <- this_excreta * params$excreta_ingesta_ratio
          }
          trackexpi <- data.frame(trackexpi, transit_clock, gut, ingesta, excreta)
        }
        return(trackexpi)
      }
      #=========================================================================

      message('--- applying poop model (first pass) ...')
      trackexp3 <- inout_model(trackexp2, params)
      if(FALSE){
        ggplot(trackexp3, aes(x=t, y=gut)) + geom_path()
        ggplot(trackexp3, aes(x=t, y=transit_clock)) + geom_path()
        ggplot(trackexp3, aes(x=t, y=ingesta)) + geom_path()
        ggplot(trackexp3, aes(x=t, y=excreta)) + geom_path()
      }

      # use that to scale down ingestion rate & recalculate ingesta/excreta distribution
      message('--- --- adjusting max. ingestion rate based on energy quotas etc ...')
      (tot_ingesta <- trackexp3$ingesta %>% sum)
      food_quota_per_trip
      (ingestion_rate_scaling <- food_quota_per_trip / tot_ingesta)
      message('--- applying poop model (second pass) ...')
      trackexp4 <- inout_model(trackexp2, params,
                               ingestion_scaling = ingestion_rate_scaling)
      if(FALSE){
        ggplot(trackexp4, aes(x=t, y=gut)) + geom_path()
        ggplot(trackexp4, aes(x=t, y=transit_clock)) + geom_path()
        ggplot(trackexp4, aes(x=t, y=ingesta)) + geom_path()
        ggplot(trackexp4, aes(x=t, y=excreta)) + geom_path()
      }
      (tot_ingesta <- trackexp4$ingesta %>% sum)
      food_quota_per_trip
      message('--- --- now total ingesta should match food quota = ',
              round(tot_ingesta, 2), ' g vs. ',round(food_quota_per_trip, 2))

      # revert name
      trackexp <- trackexp4
      trackexp %>% as.data.frame %>% head

      # Produce diagnostic plots
      message('--- producing diagnostics & maps...')
      p0 <- ggplot(trackexp, aes(x=x, y=y, color=status_detail)) + geom_point(size=.2) + xlab(NULL) + theme_minimal()
      p1 <- ggplot(trackexp, aes(x=t, y=status_detail)) + geom_point(size=.1, alpha=.5) + xlab(NULL)
      p2 <- ggplot(trackexp, aes(x=t, y=gut)) + geom_line() + xlab(NULL)
      p3 <- ggplot(trackexp, aes(x=t, y=transit_clock)) + geom_line() + xlab(NULL) +
        geom_hline(yintercept = params$ingesta_transit_time, color='firebrick', lty=2, alpha=.5)
      p4 <- ggplot(trackexp, aes(x=t, y=ingesta)) + geom_line() + xlab(NULL)
      p5 <- ggplot(trackexp, aes(x=t, y=excreta)) + geom_line() + xlab(NULL)
      ggpubr::ggarrange(p0, p1, p2, p3, p4, p5, nrow=2, ncol=3) #%>% print
      (fn_diag <- paste0('diagnostics2/',tripid,'.png'))
      ggsave(file=fn_diag, height=10, width=10, bg='white')

      # Save track map to file
      map1 <-
        ggplot(trackexp, aes(x=x, y=y, color=status_detail)) + geom_point(size=.1) +
        geom_point(data=tripi, mapping=aes(x=col_lon, y=col_lat), size=5, color='black') +
        geom_point(data=trackexp[1,], mapping=aes(x=dep_lon, y=dep_lat), size=5, color='black', alpha=.5, pch=1) +
        labs(title=tripid) + theme_minimal()

      #trackexp$excreta %>% hist
      poops <- trackexp %>% filter(excreta > 0)
      map2 <-
        ggplot(trackexp, aes(x=x, y=y)) + geom_path() +
        geom_point(data=tripi, mapping=aes(x=col_lon, y=col_lat), size=5, color='black') +
        geom_point(data=trackexp[1,], mapping=aes(x=dep_lon, y=dep_lat), size=5, color='black', alpha=.5, pch=1) +
        scale_size(breaks=seq(0.001, 0.020, length=6), range=c(.1, 6)) +
        labs(title='Predicted excreta rates:') + theme_minimal() +
        geom_point(data=poops, mapping=aes(x=x, y=y, size=excreta), color='firebrick', alpha=.02)

      ggpubr::ggarrange(map1, map2, ncol=2, widths=c(.55, .45))
      (fn_plot <- paste0('track_maps2/',tripid,'.png'))
      ggsave(file=fn_plot, height=5, width=10, bg='white')

      # Save interpolated poop tracks to file
      save(trackexp, energetics, file=fn)
      message('')

      #message('\n\nFinished!\n\n')
    } # redo analysis?
  }# only process tracks less than 2 days long
}
