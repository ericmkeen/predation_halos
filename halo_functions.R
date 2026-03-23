################################################################################
################################################################################
# The nutrient dimension of predation halos
# Eric K. Ezell, ekezell@sewanee.edu
#
# Functions & dependencies
################################################################################
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(swfscMisc)
library(tidyr)
library(ggtext)
library(truncnorm)
library(tinytable)
library(sf)
library(ggspatial)
library(ggimage)

load('gutcurve.rds')

################################################################################
################################################################################

quick_track_plot <- function(trackexp){
  p <-
    ggarrange(
      ggplot() +
        geom_path(data=trackexp,
                  mapping=aes(x=x, y=y), alpha=.6) +
        geom_point(data=trackexp %>% filter(status=='sit'),
                   mapping=aes(x=x, y=y), color='dodgerblue3') +
        theme_light(),
      ggplot() +
        geom_path(data=trackexp,
                  mapping=aes(x=t, y=km), alpha=.6) +
        geom_point(data=trackexp %>% filter(status=='sit'),
                   mapping=aes(x=t, y=km), color='dodgerblue3') +
        theme_light(),
      ggplot() +
        geom_path(data=trackexp,
                  mapping=aes(x=t, y=kmh), alpha=.6) +
        geom_point(data=trackexp %>% filter(status=='sit'),
                   mapping=aes(x=t, y=kmh), color='dodgerblue3') +
        theme_light(),
      ncol=1)
  print(p)
  return(p)
}

################################################################################
################################################################################
################################################################################
################################################################################
# DIVE SEQUENCE INSIDE A BOUT

if(FALSE){
  source('halo_parameters.R')
  params <<- params
  bouti <- 900
  dive_seq(300)
}

dive_seq <- function(bouti){
  (total_feed_time <- round(bouti * params$percent_feeding))
  (dives <- rtruncnorm(mean = params$dive_time,
                       sd = params$dive_time_sd,
                       a = 30,
                       b = 180,
                       n=10000) %>% round)
  # Determine max number of dives
  (keeps <- which(cumsum(dives) <= total_feed_time))
  max_dives <- 0
  if(length(keeps) > 0){
    dives <- dives[keeps]
    max_dives <- length(dives)
  }
  max_dives
  #message(max_dives)

  if(max_dives == 0){
    sit_seq <- rep('sit', times= bouti)
  }else{
    (post_dive_interval <- params$post_dive_interval)
    (pdi_sum <- bouti - sum(dives))
    (pdi_mean <- ceiling(pdi_sum / length(dives)))
    (dive_cycles <- dives + pdi_mean)
    sum(dive_cycles) # should now equal...
    bouti
    dive_seqs <- sapply(dives, function(x){rep('feed', times= x)})
    if(! 'list' %in% class(dive_seqs)){ dive_seqs <- list(dive_seqs)}
    dive_seqs
    (sit_seq <-
        lapply(dive_seqs,
               function(x){
                 c(x, rep('sit', times=pdi_mean))
               }
        ) %>% unlist)
    sit_seq %>% length
    if(length(sit_seq) < bouti){
      sit_seq <- c(sit_seq,
                   rep('sit', times= (bouti - length(sit_seq))))
    }
    if(length(sit_seq) > bouti){
      sit_seq <- sit_seq[1:bouti]
    }
    sit_seq %>% length
  }
  sit_seq
  return(sit_seq)
}

################################################################################
################################################################################
################################################################################
################################################################################
# BOUT SEQUENCE during a WATER SESSION

if(FALSE){
  secs <- 900
  secs <- 3600
  secs <- 6000
  toplot=TRUE
}

water_session <- function(secs, params, toplot=FALSE){
  # This function relies on the dive_seq function above

  params <<- params
  if(FALSE){ # test dive_seq()
    (dive_seq(300))
  }

  (bouts_all <- rtruncnorm(mean = params$bout_duration,
                           sd = params$bout_duration_sd,
                           a=params$dive_time + 1,
                           b=params$bout_duration + 3*params$bout_duration_sd,
                           n=10000) %>% round)
  (keeps <- which(cumsum(bouts_all) <= secs))
  n_bouts <- 1
  if(length(keeps) > 0){
    bouts <- bouts_all[keeps]
    n_bouts <- length(bouts)
  }else{
    bouts <- secs
  }
  bouts
  n_bouts
  bouts_all %>% head(10)

  (bout_interval <- params$post_dive_interval)
  (boutint_sum <- secs - sum(bouts))
  (boutint_mean <- ceiling(boutint_sum / length(bouts)))
  (bout_cycles <- bouts + boutint_mean)
  sum(bout_cycles) # should now equal...
  secs

  # Get a dive sequence for every bout
  (surf_seqs <- sapply(bouts, function(x){dive_seq(x)}))
  if(! 'list' %in% class(surf_seqs)){ surf_seqs <- list(surf_seqs)}

  (surf_seq <-
      lapply(surf_seqs,
             function(x){
               c(x, rep('sit', times=boutint_mean))
             }) %>% unlist)

  # Handle any length discrepancies
  surf_seq %>% length
  secs
  if(length(surf_seq) < secs){
    surf_seq <- c(surf_seq,
                  rep('sit', times= (secs - length(surf_seq))))
  }
  if(length(surf_seq) > secs){
    surf_seq <- surf_seq[1:secs]
  }
  surf_seq %>% length

  #return(surf_seq)

  if(toplot){
    p <- ggplot2::ggplot(data.frame(t=1:length(surf_seq),
                                    surf_seq),
                         aes(x=t, y=surf_seq)) + geom_point(cex=.2)
    print(p)
  }
  return(surf_seq)
}

# Testing
if(FALSE){
  source('halo_parameters.R')
  test <- water_session(10000, params=params, toplot=TRUE)
  test <- water_session(7000, params=params,  toplot=TRUE)
  test <- water_session(3600, params=params,  toplot=TRUE)
  test <- water_session(900, params=params,  toplot=TRUE)
  test <- water_session(700, params=params,  toplot=TRUE)
  test <- water_session(500, params=params,  toplot=TRUE)
  test <- water_session(200, params=params,  toplot=TRUE)
  test <- water_session(154, params=params,  toplot=TRUE)
  test <- water_session(100, params=params,  toplot=TRUE)
  test <- water_session(78, params=params,  toplot=TRUE)
  test <- water_session(50, params=params,  toplot=TRUE)
  test <- water_session(3,params=params,  toplot=TRUE)
}

################################################################################
################################################################################
################################################################################
################################################################################
# Gut predictor

if(FALSE){
  g_init <- NULL
  g_init <- 47
  (ti <- 60*60*6) # 6 hours -- an overnight
}

gut_predictor <- function(ti,
                          g_init = NULL,
                          gutcurve){
  gutcurvi <- gutcurve
  t_offset <- 0
  if(!is.null(g_init)){
    diffs <- abs(gutcurvi$gut - g_init)
    #plot(diffs)
    (mini <- which.min(diffs))
    (t_offset <- gutcurvi$t[mini])
  }

  guti <- 0
  (ti_adj <- ti + t_offset)
  if(ti_adj < max(gutcurvi$t)){
    diffs <- abs(gutcurvi$t - ti_adj)
    #plot(diffs)
    (mini <- which.min(diffs))
    (guti <- gutcurvi$guts[mini])
  }
  return(guti)
}

################################################################################
################################################################################
################################################################################
################################################################################
# Ingesta-Excreta model

if(FALSE){
  testi <- inout_model(trackexp2, params,
                       ingestion_scaling = 1,
                       gut0 = gut0)
  ggplot(testi, aes(x=t, y=gut)) + geom_path() + ylim(0, NA)
}

inout_model <- function(trackexp2,
                        params,
                        ingestion_scaling = 1,
                        gut0 = 0){ # allow for stomach contents to begin non-zero
  (ingesta_transit_time <- params$ingesta_transit_time)
  (gut_capacity <- params$gut_capacity)
  (excretion_rate <- params$excretion_rate)
  (max_ingestion_rate <- params$max_ingestion_rate * ingestion_scaling)
  #message(max_ingestion_rate)

  trackexpi <- trackexp2
  (stats <- trackexp2$status_detail)
  (gut <- gut0)
  (leftovers <- ifelse(gut0 > 0, TRUE, FALSE)) # beginning with gut contents
  transit_clock <- ingesta <- excreta <- c(0)
  j=1
  #for(j in 1:400){
  for(j in 1:length(stats)){
    (this_status <- stats[j])
    (pre_stats <- stats[1:j])

    # Determine transit clock
    this_transit_clock <- 0
    if(any(pre_stats == 'feed')){ # there has been some feeding on this trip
      (this_transit_clock <- ifelse(tail(gut, 1) == 0,
                                    0,
                                    tail(transit_clock, 1) + 1))
    } # this conditional zeroes transit time if starting gut contents are non-zero
    this_transit_clock

    # Calculate excreta
    this_excreta <- 0
    if(tail(gut,1) > 0){
      # gut is not empty!
      if(this_transit_clock < ingesta_transit_time){
        # but transit time is too low
        if(leftovers){
          # if the bird began trip with leftover gut contents, excrete them!
          # should not effect current transit time
          this_excreta <- tail(gut, 1)*excretion_rate
        }
      }else{
        # good to excrete -- do it!
        this_excreta <- tail(gut, 1)*excretion_rate
      }
    }

    # Calculate ingestion
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

  return(trackexpi)
}

################################################################################
################################################################################
################################################################################
################################################################################

track_interpolate <- function(tracki){
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

  return(trackexp)
}

################################################################################
################################################################################
################################################################################
################################################################################

if(FALSE){
  debug = TRUE
  i=1636
  (tripi <- trips[i, ]) %>% as.data.frame
  (tracki <- tracks %>% filter(bird_id == tripi$bird_id, tripID == tripi$tripID, deployment == tripi$deployment))
  tracki %>% tail
  nrow(tracki)
  tracki$t2[nrow(tracki)] <- NA # Guarantee that final t2 is NA
  trackexp <- track_interpolate(tracki)
  nrow(trackexp)

  sim <- list(arrange = 'random',
              pt_mean = .5,
              pt_sd = 10,
              plot = FALSE)

  tracksim <- track_simulate(trackexp, sim, debug=FALSE)
}

# ==============================================================================

track_simulate <- function(trackexp,
                           sim=NULL,
                           debug=FALSE){
  tracksim <- trackexp
  if(!is.null(sim)){
    message('--- simulating a new sit/fly sequence ...')

    trackexp %>% as.data.frame %>% head
    trackexp$ms[!is.finite(trackexp$ms)] <- NA
    (kmtot <- (trackexp$ms %>% sum(na.rm=TRUE))/1000)
    (ttot <- trackexp$t %>% range %>% diff)
    if(debug){
      quick_track_plot(trackexp)
    }

    message('    --- summarizing each stage of the trip ...')
    trackmod <-
      trackexp %>%
      # Create group IDs for each status session
      filter(is.finite(ms)) %>%
      mutate(t_trip = 0:(n()-1)) %>%
      mutate(t_prop = t_trip / max(t_trip)) %>%
      mutate(kmcum = cumsum(ms)/1000) %>%
      mutate(status_last = lag(status, 1)) %>%
      mutate(status_last = ifelse(is.na(status_last), status, status_last)) %>%
      mutate(status_shift = ifelse(status_last != status,
                                   1, 0)) %>%
      mutate(status_id = cumsum(status_shift) + 1) %>%
      as.data.frame

    trackmod %>% head
    if(debug){
      trackmod$status_id %>% plot
      trackmod$kmcum %>% plot
    }

    stati <-
      trackmod %>%
      group_by(status_id, status) %>%
      summarize(n=n(),
                t1 = t_trip[1],
                tn = t_trip[n()],
                km1 = kmcum[1],
                kmn = kmcum[n()],
                ms = mean(ms, na.rm=TRUE),
                .groups='drop') %>%
      mutate(tdur = tn - t1) %>%
      mutate(pdur = tdur / ttot) %>%
      mutate(kmdur = kmn - km1) %>%
      mutate(pkm = kmdur / kmtot)
    stati

    (sits <- stati %>% filter(status == 'sit'))
    sits
    (nsit <- sits %>% nrow)
    (tsit <- sits$tdur %>% sum)
    (psit <- tsit / ttot)
    (tfly <- ttot - tsit)
    (pfly <- 1 - psit)
    # mean flight speed
    (ms_mean <- trackmod %>% filter(status=='fly') %>% pull(ms) %>% mean(na.rm=TRUE))

    if(any(stati$status == 'fly') & nsit > 0){

      # Re-order  =========================================================
      message('    --- reordering the sit sessions ...')
      if(sim$arrange == 'random'){
        (siti <- sits[sample(nrow(sits)), ])
      }
      if(sim$arrange == 'frontload'){
        (siti <- sits %>% tibble %>% arrange(desc(tdur)))
      }
      if(sim$arrange == 'backload'){
        (siti <- sits %>% tibble %>% arrange(tdur))
      }
      siti

      # Decide how to space them in time ==================================
      # message('    --- re-spacing those sit sessions in time ...')
      # siti
      # (nflys <- nrow(siti) + 1)
      # pfly
      # (pfly_left <- pfly)

      # Decide how to space them in time ==================================
      message('    --- re-spacing those sit sessions in time ...')
      (next_min_t <- 0)
      (pt_wiggle <- pfly)
      pt1_news <- ptn_news <- c()
      si=2
      for(si in 1:nrow(siti)){
        (sitii <- siti[si, ])
        (pduri <- sitii$pdur)
        (pt_min <- next_min_t)
        # Give enough wiggle room to account for all the wiggles to come
        (pt_max <- (pt_min + (pt_wiggle/(nrow(siti) + 1 - si))))
        (pt_max <- min(c(pt_max, 1)))
        # Randomy place new start time based on a rnorm
        if(pt_min == pt_max){
          pt1_new <- pt_min
        }else{
          (pt1_new <- rtruncnorm(n=1,
                                 a=pt_min,
                                 b=pt_max,
                                 mean=sim$pt_mean,
                                 sd=sim$pt_sd))
        }
        pt1_new
        (pt_wiggle <- pt_wiggle - (pt1_new - pt_min))
        if(pt_wiggle <= 0){pt_wiggle <- 0}
        (ptn_new <- pt1_new + pduri)
        next_min_t <- ptn_new
        pt1_news <- c(pt1_news, pt1_new)
        ptn_news <- c(ptn_news, ptn_new)
        #message(pt1_news)
      }
      pt1_news
      ptn_news
      siti$pt1_new <- pt1_news
      siti$ptn_new <- ptn_news
      siti

      (sitsim <-
          siti %>%
          tibble %>%
          select(-t1, -tn, -km1, -kmn) %>%
          #select(-ms) %>%
          rename(pt1 = pt1_new, ptn = ptn_new) %>%
          as.data.frame)

      (ttot_sit <- sitsim$tdur %>% sum)
      ttot
      (kmtot_sit <- sitsim$kmdur %>% sum)
      kmtot
      (ttot_fly <- ttot - ttot_sit)
      (kmtot_fly <- kmtot - kmtot_sit)
      (kmh <- kmtot_fly / (ttot_fly/3600))
      (ms_mean_new <- kmh * 0.277778)
      ms_mean

      # Add flights in between  ===========================================
      message('    --- adding flight stages in between ...')
      sitsim

      (flitsim <-
          sitsim %>%
          tibble %>%
          mutate(ptn_lag = lag(ptn)) %>%
          mutate(ptn_lag = ifelse(is.na(ptn_lag), 0, ptn_lag)) %>%
          mutate(pt1_fly = ptn_lag,
                 ptn_fly = pt1) %>%
          mutate(status = 'fly') %>%
          mutate(pt1 = pt1_fly,
                 ptn = ptn_fly) %>%
          select(-n, -ptn_lag, - pt1_fly, - ptn_fly) %>%
          mutate(tdur = (ptn - pt1)*ttot) %>%
          mutate(mdur = tdur * ms_mean_new)  %>%
          mutate(kmdur = mdur / 1000) %>%
          mutate(ms = ms_mean_new)
      )

      sitsim
      flitsim
      (newsim <-
          rbind(sitsim %>% select(-n),
                flitsim %>% select(-mdur)
          ) %>%
          arrange(pt1) %>%
          as.data.frame)
      newsim

      # make sure this agenda goes all the way to 1
      if(is.na(tail(newsim, 1)$ptn)){newsim$ptn[nrow(newsim)] <- 1}
      newsim
      if(tail(newsim, 1)$ptn < 1){
        (simn <- data.frame(status_id = max(newsim$status_id) + 1,
                            status = 'fly',
                            ms = ms_mean,
                            pt1 = tail(newsim, 1)$ptn,
                            ptn = 1) %>%
           mutate(pdur = (ptn - pt1)) %>%
           mutate(tdur = ttot*pdur) %>%
           mutate(ms = ms_mean_new) %>%
           mutate(kmdur = tdur*ms/1000) %>%
           mutate(pkm = kmdur / kmtot) %>%
           select(status_id, status, ms, tdur, pdur, kmdur, pkm, pt1, ptn)
        )
        newsim <- rbind(newsim, simn)
      }
      newsim
      if(any(newsim$ptn > 1)){
        (max_ptn <- max(newsim$ptn, na.rm=TRUE))
        newsim$ptn <- newsim$ptn / max_ptn
        newsim$pt1 <- newsim$pt1 / max_ptn
      }
      newsim

      # figure out distance traveled along route based on average speed of stage
      message('    --- determing distance traveled in each stage ...')
      newsim$kmdur %>% sum
      kmtot
      (newsim <-
          newsim %>%
          mutate(ptdur = ptn - pt1) %>%
          mutate(tdur = ptdur*ttot) %>%
          mutate(mdur = kmdur*1000) %>%
          #mutate(mdur = tdur*ms) %>%
          mutate(mn = cumsum(mdur)) %>%
          mutate(m1 = lag(mn, default=0)) %>%
          mutate(status_id = 1:n()) #%>%
        #select(status_id, status:mdur, m1, mn)
      )

      # Re-time the track =================================================
      message('    --- re-timing the track...')
      trackmod$newid <- 1:nrow(trackmod)
      tracknew <- data.frame()
      si=1
      si=nrow(newsim)
      for(si in 1:nrow(newsim)){
        (simi <- newsim[si,])
        (stati <- simi$status)
        (m1 <- simi$m1)
        (mn <- simi$mn)
        (mprops <- trackmod$ms %>% cumsum)
        mprops %>% tail

        # find index for t1
        (diffs <- abs(mprops - m1))
        #diffs %>% plot
        (m1i <- which.min(diffs))
        # find index for tn
        (diffs <- abs(mprops - mn))
        #diffs %>% plot
        (mni <- which.min(diffs))

        # new time intervals between track records
        (tracki <- trackmod[m1i:mni, ]) %>% head
        tracki %>% nrow
        message('        --- time scaling factor: ', round((simi$tdur / nrow(tracki)), 5))
        (tseg <- seq(0, simi$tdur, length=nrow(tracki)))
        tracki$tseg <- tseg
        tracki$segid <- si

        # replace other values in the track
        tracki$ms <- simi$ms
        tracki$kmh <- simi$ms * 3.6
        tracki$status <- simi$status
        tracki$status_id <- simi$status_id
        tracki <- tracki %>% tibble %>%
          select(ms, kmh, status, x, y, km, tseg, status_id, newid)
        tracki %>% head

        tracknew <- rbind(tracknew, tracki)
        #
      }

      tracknew %>% head

      # handle duplicated trackmod rows
      tracknew %>% nrow
      tracknew %>% group_by(newid) %>% tally() %>% pull(n) %>% table
      tracknew <-
        tracknew %>%
        group_by(newid) %>%
        slice(1)
      tracknew %>% nrow

      # ensure proportional time allocation is in tact
      tracknew %>%
        group_by(status_id, status) %>%
        summarize(dur = max(tseg),
                  .groups='drop') %>%
        group_by(status) %>%
        summarize(dur = sum(dur)/ttot, .groups='drop')

      # setting a master time marker
      tracknew %>% head
      tracknew <-
        tracknew %>%
        group_by(status_id) %>%
        arrange(tseg) %>%
        mutate(tlag = lag(tseg, default=0)) %>%
        mutate(tdiff = tseg - tlag)  %>%
        ungroup() %>%
        arrange(status_id, tseg) %>%
        mutate(t_trip = cumsum(tdiff))

      # ensure proportional time allocation is in tact
      tracknew %>% head
      tracknew %>%
        group_by(status_id, status) %>%
        summarize(dur = t_trip[n()] - t_trip[1],
                  .groups='drop') %>%
        group_by(status) %>%
        summarize(dur = sum(dur)/ttot,
                  .groups='drop')

      if(debug){
        tracknew$t_trip %>% head
        tracknew$t_trip %>% plot
        tracknew$status_id %>% plot
      }
      nrow(trackmod)
      nrow(tracknew)

      # Simplify to one row per second  ===================================
      message('    --- reducing to one record per second ...')
      tracknew %>% head
      tracknew <-
        tracknew %>%
        mutate(tr = round(t_trip)) %>%
        arrange(tr) %>%
        group_by(tr, status_id, status) %>%
        summarize(ms = mean(ms, na.rm=TRUE),
                  kmh = mean(kmh, na.rm=TRUE),
                  km = mean(km),
                  x = mean(x),
                  y = mean(y),
                  .groups='drop') %>%
        ungroup()

      tracknew %>% head
      nrow(tracknew)
      tracknew %>% group_by(status_id, status) %>% tally

      # ensure proportional time allocation is in tact
      tracknew %>%
        group_by(status_id, status) %>%
        summarize(dur = tr[n()] - tr[1],
                  .groups='drop') %>%
        group_by(status) %>%
        summarize(dur = sum(dur)/ttot,
                  .groups='drop')

      # setup leading values for interpolation,
      # remove redundant timestamps
      tracknew <-
        tracknew %>%
        arrange(tr) %>%
        mutate(t2 = lead(tr)) %>%
        mutate(secs = t2 - tr) %>%
        mutate(x2 = lead(x),
               y2 = lead(y)) %>%
        mutate(x2 = ifelse(is.na(x2), x, x2)) %>%
        mutate(y2 = ifelse(is.na(y2), y, y2)) %>%
        mutate(km2 = lead(km, 1)) %>%
        mutate(km2 = ifelse(is.na(km2), km, km2)) %>%
        filter(!is.na(t2)) %>%
        filter(secs != 0)

      #tracknew$tr %>% plot
      tracknew$secs %>% table
      tracknew$secs %>% sum

      # Re-interpolate times  =============================================
      message('--- reinterpolating in case any seconds were missed ...')
      tracknew %>% head
      tracknewi <-
        tracknew %>%
        rowwise() %>%
        mutate(tnew = ifelse(secs > 1,
                             list(tr:t2),
                             list(tr))) %>%
        mutate(xnew = ifelse(secs > 1,
                             list(seq(x, x2, length= (secs + 1))),
                             list(x))) %>%
        mutate(ynew = ifelse(secs > 1,
                             list(seq(y, y2, length= (secs + 1))),
                             list(y))) %>%
        mutate(kmnew = ifelse(secs > 1,
                              list(seq(km, km2, length = (secs + 1))),
                              list(km))) %>%
        unnest(cols = c(tnew, xnew, ynew, kmnew)) %>%
        ungroup()

      tracknewi %>% nrow

      # finalize, format and setup columns
      tracknewi <-
        tracknewi %>%
        distinct(tnew, .keep_all=TRUE) %>%
        arrange(tr) %>%
        mutate(t_trip = 1:n()) %>%
        mutate(t = trackmod$t[1] + t_trip) %>%
        mutate(t_prop = t_trip / t_trip[n()]) %>%
        mutate(q = cumsum(ms) / sum(ms)) %>%
        mutate(status_last = lag(status, 1)) %>%
        mutate(status_shift = ifelse(status_last != status, 1, 0)) %>%
        select(ms, kmh, status, t, x=xnew, y=ynew, km=kmnew,
               q, t_prop, status_last, status_shift,
               status_id)

      tracknewi %>% head
      tracknewi %>% group_by(status_id, status) %>% tally

      # Re-constitute dataframe  ==========================================
      message('--- reconstituting column names ...')
      suppressWarnings({
        tracksim <-
          data.frame(trackmod %>% select(site:tripID) %>% head(1),
                     tracknewi)
      })
      tracksim %>% head

      tracksim %>% group_by(status_id, status) %>%  tally()

      # Compare old to new
      if(sim$plot | debug){
        ggarrange(quick_track_plot(trackexp),
                  quick_track_plot(tracksim),
                  ncol=2) %>% print
      }

      (kmtot_sim <- (tracksim$ms %>% sum(na.rm=TRUE))/1000)
      (ttot_sim <- tracksim$t %>% range %>% diff)
      (psit_sim <- (tracksim %>% filter(status=='sit') %>% nrow()) / ttot_sim)

      message('--- Finished!')
      message('    --- Original track: ',
              round(kmtot, 1), ' km, ',
              round(ttot/3600, 1), ' hrs, ',
              round(psit*100, 1), '% sitting')
      message('    --- Simulated track: ',
              round(kmtot_sim, 1), ' km, ',
              round(ttot_sim/3600, 1), ' hrs, ',
              round(psit_sim*100, 1), '% sitting')

    } # end if at least 1 sit
  } # end if null sim
  return(tracksim)
}

################################################################################
################################################################################
################################################################################
################################################################################

trip_model <- function(trips_to_process,
                       params,
                       trips,
                       tracks,
                       sim = NULL,
                       tracks_dir = 'poop_tracks/',
                       diagnostics_dir = 'diagnostics/',
                       maps_dir = 'track_maps/',
                       redo_analysis = FALSE,
                       plots=TRUE){

  if(FALSE){ # debugging
    trips_to_process = 1:100
    i=1
    i=100
    i=3000

    (trips_to_process = sample(1:nrow(trips), size=10))
    tracks_dir = 'sim/poop_tracks_centered/'
    diagnostics_dir = 'diagnostics/'
    maps_dir = 'track_maps/'
    redo_analysis = TRUE

    i=1
    i=3000
    sim <- list(arrange = 'backload',
                pt_mean = .5,
                pt_sd = .5,
                plot = TRUE)
  }

  load('gutcurve.rds')

  # core ===========
  # interpolate
  # feeding_sequence
  # energy needs
  # inout model

  for(i in trips_to_process){
    (tripi <- trips[i, ]) %>% as.data.frame

    if(TRUE){
      (tripid <- paste0(tripi$site,'_',tripi$year,'_',
                        tripi$bird_id,'_',tripi$tripID,
                        tripi$deployment))
      (fn <- paste0(tracks_dir,tripid,'.rds'))   # stage filename for result

      if(!file.exists(fn) | redo_analysis){
        message('======= Trip ', tripid,' (',i,' of ', nrow(trips),') =======')

        (tracki <- tracks %>% filter(bird_id == tripi$bird_id, tripID == tripi$tripID, deployment == tripi$deployment))
        tracki #%>% head
        tracki %>% tail
        nrow(tracki)
        #tracki$secs %>% plot

        # Handle anomalous file errors
        tracki$t2[nrow(tracki)] <- NA # Guarantee that final t2 is NA
        if(i == 1175){tracki <- tracki[1:214, ]}
        if(i == 1176){tracki <- tracki[1:11, ]}

        # Interpolate track ======================================================

        message('--- interpolating track ...')
        trackexp <- track_interpolate(tracki)

        if(FALSE){
          quick_track_plot(trackexp)
        }
        trackexp %>% as.data.frame %>% head
        trackexp %>% as.data.frame %>% tail(50)
        trackexp %>% nrow

        # fix NA status
        nas <- which(is.na(trackexp$status))
        if(length(nas)>0){
          trackexp$status[nas] <- 'error'
        }


        # Add feeding simulation? ==============================================
        if(!is.null(sim)){
          tracksim <- track_simulate(trackexp = trackexp,
                                     sim=sim,
                                     debug=FALSE)
          trackexp <- tracksim
        }


        # splashdowns ============================================================

        # https://www.int-res.com/articles/meps2013/475/m475p277.pdf
        # mean is 320m from colony, sd 126, 320 + 2 SDs = 572m
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

        # Creating feeding sequence ==============================================

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
            trackexpi$status_detail <- water_session(secs, params)
          }
          trackexp2 <- rbind(trackexp2, trackexpi)
        }

        trackexp2$status_detail %>% table
        if(FALSE){
          ggplot2::ggplot(trackexp2, aes(x=t, y=status_detail)) + geom_point(cex=.1)
        }

        # Get starting point energy needs ========================================
        message('--- calculating energy needs ...')

        (fly_per_trip <- trackexp2 %>% filter(status_detail == 'fly') %>% nrow)
        (sit_per_trip <- trackexp2 %>% filter(status_detail %in% c('sit', 'splashdown')) %>% nrow)
        (feed_per_trip <- trackexp2 %>% filter(status_detail == 'feed') %>% nrow)
        (error_per_trip <- trackexp2 %>% filter(status_detail == 'error') %>% nrow)

        (duration <- tripi$duration)
        (n_trips_exact <- params$trip_hours / (duration*2)) # exact
        (n_trips <- ceiling(n_trips_exact)) # rounded

        # Play with new approach based on Burke & Montevecchi (2018)
        # which is based on Elliot & Gaston 2014 (TBMU) and applied to COMU
        # DEE = 508*Tf + 1.01 SUM [ 1 - e^(Td / 1.23)] + (113 - 2.75*SST)*Ts + (72.2 - 2.75*SST)*Tw + 33Tc
        # Tf = hours flying per day
        # Td = Duration = dive duration in **minutes**, sum of all dives in a day
        # Ts = hours spent active on the water in a day
        # Tw = hours spent inactive on the water
        # Fw = fraction of sit time spent inactive on the water
        # Tc = hours spent at colony

        # Per-trip time allocation =============================================

        (SST <- tripi$sst)
        (Fw = params$percent_inactive)
        (Tf <- fly_per_trip / 3600) # hours
        (Td <- feed_per_trip / 60) # minutes
        Td / 60 # in hours fyi
        (Ts <- (sit_per_trip/3600) * (1 - Fw))
        (Tw <- (sit_per_trip/3600) * Fw)
        Ts + Tw + (Td/60) + Tf # should be same as trip length

        # Trip Energy Expenditure (TEE) ========================================

        # Compute energy spent across dives
        trackexp2 %>% as.data.frame %>% head
        trackexp2$status_detail
        trackexp2$status_detail %>% table
        dives <- trackexp2
        dives %>% as.data.frame %>% head
        dives <- dives %>%
          mutate(status_last = lag(status_detail)) %>%
          mutate(status_shift = ifelse(status_last != status_detail & status_detail == 'feed', 1, 0)) %>%
          mutate(status_shift = ifelse(is.na(status_shift) & status_detail == 'feed', 1, status_shift)) %>%
          mutate(dive_id = cumsum(status_shift)) %>%
          filter(status_detail == 'feed') %>%
          group_by(dive_id) %>%
          summarize(Td = n()) %>%
          mutate(Td = Td/60) %>%
          mutate(exp_term = exp((-1*Td)/1.23)) %>%
          mutate(one_minus_exp = 1 - exp_term)
        dives %>% nrow
        dives$Td %>% sum
        Td # should match
        (Ed <- 1.01 * (sum(dives$one_minus_exp)))

        # energy spent on this trip
        (TEE <- 508*Tf + Ed + (113 - 2.75*SST)*Ts + (72.2 - 2.75*SST)*Tw)

        # energy spent for same duration during next brood shift
        (BEE <- 33*min(c(18 - duration, duration)))

        # energy spent during overnight rest at colony
        (CEE_tot <- 33*(24 - params$trip_hours))
        (n_trips <- ceiling(params$trip_hours / (duration*2)))
        (CEE_trip = CEE_tot / n_trips)

        (kJ_per_trip <- TEE + BEE + CEE_trip)

        # Translate to ingestion quota and mean ingestion rate ===================
        (food_per_trip <- (kJ_per_trip / params$assimilation_efficiency) / params$E_density_diet) # grams
        (ingest_rate <- food_per_trip / feed_per_trip) # grams/sec
        ingest_rate_raw <- ingest_rate

        # is this higher than the max possible? If so, replace with max possible.
        params$max_ingestion_rate
        if(ingest_rate > params$max_ingestion_rate){ ingest_rate <- params$max_ingestion_rate }
        ingest_rate

        #=========================================================================
        # figure out how much would be eaten if using the physiological maximum ingestion rate
        # here assume initial gut contents are 0 (to maximize ingestion rate further)
        message('--- applying poop model (first pass) ...')
        trackexp3 <- inout_model(trackexp2, params,
                                 ingestion_scaling = 1,
                                 gut0 = 0)
        if(FALSE){
          ggplot(trackexp3, aes(x=t, y=gut)) + geom_path() + ylim(0, NA)
          ggplot(trackexp3, aes(x=t, y=transit_clock)) + geom_path()
          ggplot(trackexp3, aes(x=t, y=ingesta)) + geom_path()
          ggplot(trackexp3, aes(x=t, y=excreta)) + geom_path()
        }

        #=========================================================================
        # use that to scale down ingestion rate & recalculate ingesta/excreta distribution

        message('--- --- adjusting max. ingestion rate based on energy quotas etc ...')
        (tot_ingesta <- trackexp3$ingesta %>% sum)
        food_per_trip
        (ingestion_rate_scaling <- food_per_trip / tot_ingesta)
        message('--- applying poop model (second pass) ...')
        trackexp4 <- inout_model(trackexp2, params,
                                 ingestion_scaling = ingestion_rate_scaling,
                                 gut0 = 0)
        if(FALSE){
          ggplot(trackexp4, aes(x=t, y=gut)) + geom_path()
          ggplot(trackexp4, aes(x=t, y=transit_clock)) + geom_path()
          ggplot(trackexp4, aes(x=t, y=ingesta)) + geom_path()
          ggplot(trackexp4, aes(x=t, y=excreta)) + geom_path()
        }
        (tot_ingesta <- trackexp4$ingesta %>% sum)
        food_per_trip
        message('--- --- now total ingesta should match food quota = ',
                round(tot_ingesta, 2), ' g vs. ',round(food_per_trip, 2))

        #=========================================================================
        # Use this revised model to include initial gut contents

        message('--- --- re-running, now  adding initial gut contents ...')
        # grams left in stomach at end of trip
        (gut_from_last <- trackexp4$gut %>% tail(1))

        # get last sit record
        (sits <- which(trackexp4$status == 'sit'))
        if(length(sits) > 0){
          (last_sit <- sits %>% tail(1))
        }else{
          last_sit <- nrow(trackexp4)
        }
        (km_last_sit <- trackexp4$km[last_sit])

        # Get average flight speed
        (flits <- which(trackexp4$status == 'fly'))
        if(length(flits) > 0){
          (mean_flight_kmh <- trackexp4$kmh[flits] %>% mean(na.rm=TRUE))
        }else{
          (mean_flight_kmh <- 60.12) # dataset-wide mean
        }
        mean_flight_kmh

        # estimate flight time back to colony from last sit
        (flyback_hrs <- km_last_sit / mean_flight_kmh)

        # to this, add the expectation that
        # the bird's mate will make a foraging trip of the same duration
        (hrs_since_last_eat <- flyback_hrs + tripi$duration)

        # now add overnight brood rest,
        # in the event that this is the first flight of the day for the pair.
        # use a random test to determine if it is,
        # based on number of trips expected during day for the mates combined
        #(pair_trips <- 2*n_trips)
        (rando_threshold <- (1 - (1/n_trips)))
        #(rando_threshold <- (1 - 1/pair_trips))
        (rando <- runif(1))
        (first_trip <- ifelse(rando > rando_threshold, TRUE, FALSE))

        # if this is first trip of day, add overnight gap to time since last eat
        if(first_trip){
          (hrs_overnight <- 24 - params$trip_hours)
          (hrs_since_last_eat <- hrs_since_last_eat + hrs_overnight)
        }
        hrs_since_last_eat
        (t_since_last_eat <- hrs_since_last_eat * 3600)

        # Now estimate contents left in gut at start of foraging trip
        (gut0 <- gut_predictor(t_since_last_eat,
                               g_init = gut_from_last,
                               gutcurve = gutcurve))

        message('--- applying poop model (third pass) ...')
        trackexp4 <- inout_model(trackexp2, params,
                                 ingestion_scaling = ingestion_rate_scaling,
                                 gut0 = gut0)
        if(FALSE){
          ggplot(trackexp4, aes(x=t, y=gut)) + geom_path()
          ggplot(trackexp4, aes(x=t, y=transit_clock)) + geom_path()
          ggplot(trackexp4, aes(x=t, y=ingesta)) + geom_path()
          ggplot(trackexp4, aes(x=t, y=excreta)) + geom_path()
        }

        #=========================================================================

        (energetics <- data.frame(tripi,
                                  fly_per_trip,
                                  sit_per_trip,
                                  feed_per_trip,
                                  error_per_trip,
                                  n_trips,
                                  TEE,
                                  Tf, Td, Tw, Ts, BEE, CEE_trip,
                                  kJ_per_trip,
                                  food_per_trip,
                                  ingest_rate_raw,
                                  ingest_rate_pass1 = ingest_rate,
                                  ingestion_rate_scaling,
                                  ingest_rate_pass2 = ingest_rate*ingestion_rate_scaling,
                                  gut0,
                                  data.frame(params)))

        # revert name
        trackexp <- trackexp4
        trackexp %>% as.data.frame %>% head

        #=========================================================================
        # Produce diagnostic plots

        if(plots){

          message('--- producing diagnostics & maps...')
          p0 <- ggplot(trackexp, aes(x=x, y=y, color=status_detail)) + geom_point(size=.2) + xlab(NULL) + theme_minimal()
          p1 <- ggplot(trackexp, aes(x=t, y=status_detail)) + geom_point(size=.1, alpha=.5) + xlab(NULL)
          p2 <- ggplot(trackexp, aes(x=t, y=gut)) + geom_line() + xlab(NULL)
          p3 <- ggplot(trackexp, aes(x=t, y=transit_clock)) + geom_line() + xlab(NULL) +
            geom_hline(yintercept = params$ingesta_transit_time, color='firebrick', lty=2, alpha=.5)
          p4 <- ggplot(trackexp, aes(x=t, y=ingesta)) + geom_line() + xlab(NULL)
          p5 <- ggplot(trackexp, aes(x=t, y=excreta)) + geom_line() + xlab(NULL)
          ggpubr::ggarrange(p0, p1, p2, p3, p4, p5, nrow=2, ncol=3) #%>% print
          (fn_diag <- paste0(diagnostics_dir,tripid,'.png'))
          ggsave(file=fn_diag, height=10, width=10, bg='white')

          #=======================================================================
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
          (fn_plot <- paste0(maps_dir,tripid,'.png'))
          ggsave(file=fn_plot, height=5, width=10, bg='white')
        } # end if PLOTS

        # Save interpolated poop tracks to file
        save(trackexp, energetics, file=fn)
        message('')
      } # redo analysis?
    }
  } # end of loop
  message('\n\nFinished!\n\n')
}


################################################################################
################################################################################
################################################################################
################################################################################
# Ingestion-excretion halo maps

if(FALSE){

  trips %>% filter(colony_size == 3701) %>% pull(site)
  trak <- tracks %>% filter(site == 'Papey')
  ggplot(trak, aes(x=lon, y=lat, group=tripID)) +
    geom_point(size=.1, alpha=.7) +
    geom_path(alpha=.5, linewidth=.4) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave('track_base.png', width=8, height=8)

  # confirm each site has only one colony size
  trips %>% group_by(site, colony_size) %>% tally() %>% group_by(site) %>% tally() %>% as.data.frame
  trips %>% group_by(site, colony_size) %>% tally() %>% arrange(colony_size)
  trips %>% group_by(colony_size) %>% tally %>% arrange(desc(n)) %>% as.data.frame

  colsize = 2672
  scale_pos = 'br'
  labels='auto'
  image_xy <- c(-1.8, 57.36)
  xlims = c(NA, NA)
  ylims = c(NA, NA)

  # try it
  halo_map(2295)
  halo_map(2672, image_xy = c(-1.8, 57.36))

  # see sample sizes for all colony sizes
  colsizes <-
    trips %>%
    group_by(colony_size) %>%
    tally %>% arrange(n) %>%
    pull(colony_size)

  # Loop through all options and save previews
  #sizi = colsizes[1]
  #for(sizi in colsizes){
  #  halo_map(sizi) %>% print
  #  ggsave(paste0('track_map_',sizi, '.png'),
  #         width=14, height=4, bg='white')
  #  message(sizi)
  #  #readline()
  #}
}

#===============================================================================

halo_map <- function(colsize,
                     scale_pos = 'br',
                     xlims = NULL,
                     ylims = NULL,
                     image_xy = NULL,
                     image_size = .3,
                     col_xy = NULL,
                     col_color = 'black',
                     col_size = 3,
                     scale_width = .4,
                     labels = 'auto'){
  (tripsub <- trips %>% filter(colony_size == colsize))
  tripsub$site %>% unique
  (sites <- tripsub$site %>% unique %>% head(1))

  coarse_trax <- tracks %>% filter(site %in% sites)
  coarse_trax %>% head

  (lf <- list.files('poop_tracks/'))
  (keeps <- grep(sites, lf))
  trax <- data.frame()
  i=2250
  for(i in keeps){
    (lfi <- lf[i])
    message(lfi)
    load(paste0('poop_tracks/', lfi))
    trax <- rbind(trax, trackexp)
  }

  trax %>% nrow
  trax %>% head

  traxi <-
    trax %>%
    mutate(id = paste0(bird_id,'_', deployment ,'_', tripID)) %>%
    group_by(id) %>%
    mutate(ing_tot = ingesta %>% sum(na.rm=TRUE),
           exc_tot = excreta %>% sum(na.rm=TRUE)) %>%
    arrange(t) %>%
    mutate(ing_cum = ingesta %>% cumsum,
           exc_cum = excreta %>% cumsum) %>%
    rowwise %>%
    mutate(qing = ing_cum / ing_tot,
           qexc = exc_cum / exc_tot) %>%
    mutate(qing_diff = abs(.5 - qing),
           qexc_diff = abs(.5 - qexc)) %>%
    ungroup


  ids <- traxi %>%
    group_by(id) %>%
    summarize(maxkm = max(km, na.rm=TRUE)) %>%
    arrange(desc(maxkm)) %>%
    pull(id)

  # Tracks =====================================================================

  traxi %>% names

  traxisf <- sf::st_as_sf(coarse_trax, coords = c("lon", "lat"), crs = 4326)
  traxisf %>% as.data.frame %>% head

  (lines <- traxisf %>%
      mutate(id = paste0(bird_id, '_', tripID)) %>%
      group_by(id) %>%
      summarize(do_union = FALSE) %>%
      st_cast("LINESTRING"))

  traxi %>% nrow
  traxi %>% as.data.frame %>% head

  if(FALSE){
    (traxip <-
       traxi %>%
       group_by(id) %>%
       mutate(t1 = t[1]) %>%
       mutate(maxkm = max(km, na.rm=TRUE)) %>%
       mutate(ttrip = t - t1) %>%
       mutate(qkm = km / maxkm) %>%
       ungroup)

    ggplot(traxip,
           aes(x=ttrip, y=km, group=id)) +
      geom_path()

    ggplot(traxip,
           aes(x=q, y=qkm, group=id)) +
      geom_path()

    ggplot() +
      geom_path(data=traxip,
                mapping=aes(x=q, y=qing, group=id),
                color='blue') +
      geom_path(data=traxip,
                mapping=aes(x=q, y=qexc, group=id),
                color='firebrick')

    ggplot() +
      geom_path(data=traxip,
                mapping=aes(y=qkm, x=qing, group=id),
                color='blue') +
      geom_path(data=traxip,
                mapping=aes(y=qkm, x=qexc, group=id),
                color='firebrick')
  }


  p1 <- ggplot() + theme_light() +
    geom_sf(data=traxisf, mapping=aes(group=tripID),
            alpha=.3, size=.3)

  if(!is.null(image_xy)){
    imgdf <- data.frame(lon=image_xy[1], lat = image_xy[2])
    imgsf <- sf::st_as_sf(imgdf, coords = c("lon", "lat"), crs = 4326)
    p1 <- p1 +
      geom_image(data=imgdf,
                 aes(x=lon,
                     y=lat,
                     image = 'murre-edit copy.png'),
                 size = image_size) +
      xlab(NULL) + ylab(NULL)
  }

  p1 <- p1 +
    geom_sf(data=lines, alpha=.3) +
    coord_sf() +
    annotation_scale(location = scale_pos, width_hint = scale_width) +
    labs(title=paste0('Site: ', traxi$site[1],
                      '  |  Colony: ', prettyNum(colsize, big.mark = ",")),
         subtitle = paste0('n = ', nrow(tripsub),' trips (Patterson et al. 2022)')) +
    theme(plot.title = element_text(face = "bold"))

  if(!is.null(xlims)){
    p1 <- p1 + xlim(xlims[1], xlims[2])
  }
  if(!is.null(ylims)){
    p1 <- p1 + ylim(ylims[1], ylims[2])
  }
  #if(!is.null(col_xy)){
  #  p1 <- p1 + geom_point(mapping=aes(x=col_xy[1], y=col_xy[2]),
  #                        pch=8, color=col_color, size=col_size) +
  #    xlab(NULL) + ylab(NULL)
  #}


  # Ingestion ==================================================================
  # Add ids from farthest out in
  ids
  p2 <- ggplot() + theme_light() + coord_sf()
  for(idi in ids){
    ti <- traxi %>% filter(id == idi)
    message(idi)

    # subsample
    subsampler <- rep(c(TRUE, rep(FALSE, times=10)), times=nrow(ti))
    #subsampler
    subsampler <- subsampler[1:nrow(ti)]
    ti <- ti[subsampler, ]
    tisf <- sf::st_as_sf(ti, coords = c("x", "y"), crs = 4326)

    # plot
    p2 <-
      p2 +
      geom_sf(data = tisf,
              mapping = aes(geometry = geometry,
                            group = id,
                            color = qing_diff,
                            alpha = .1*(.5 - qing_diff)),
              size=.3)
  }

  p2 <-
    p2 +
    viridis::scale_color_viridis(direction=1, option='mako', breaks=c(.5, .4, .3, .2, .1, 0)) +
    annotation_scale(location = scale_pos, width_hint = scale_width) +
    guides(color = guide_colourbar(reverse = TRUE)) +
    guides(alpha = "none") +
    labs(color='Distance\nfrom\ningestion\nmidpoint',
         title = ' ',
         subtitle = 'Ingestion model')

  if(!is.null(xlims)){
    p2 <- p2 + xlim(xlims[1], xlims[2])
  }
  if(!is.null(ylims)){
    p2 <- p2 + ylim(ylims[1], ylims[2])
  }
  if(!is.null(col_xy)){
    p2 <- p2 + geom_point(mapping=aes(x=col_xy[1], y=col_xy[2]),
                          pch=8, color=col_color, size=col_size) +
      xlab(NULL) + ylab(NULL)
  }

  #p2

  # Excretion =========================================================
  # Add ids from farthest out in
  ids
  p3 <- ggplot() + theme_light()
  for(idi in ids){
    ti <- traxi %>% filter(id == idi)
    message(idi)

    # subsample
    subsampler <- rep(c(TRUE, rep(FALSE, times=10)), times=nrow(ti))
    #subsampler
    subsampler <- subsampler[1:nrow(ti)]
    ti <- ti[subsampler, ]
    tisf <- sf::st_as_sf(ti, coords = c("x", "y"), crs = 4326)
    p3 <- p3 +
      geom_sf(data = tisf,
              mapping = aes(geometry = geometry,
                            group = id,
                            color = qexc_diff,
                            alpha = .1*(.5 - qexc_diff)),
              size=.3)
  }

  p3 <-
    p3 +
    viridis::scale_color_viridis(direction=1, option='rocket', breaks=c(.5, .4, .3, .2, .1, 0)) +
    annotation_scale(location = scale_pos, width_hint = scale_width) +
    guides(color = guide_colourbar(reverse = TRUE)) +
    guides(alpha = "none") +
    labs(color='Distance\nfrom\nexcretion\nmidpoint',
         title = ' ',
         subtitle = 'Excretion model')

  if(!is.null(xlims)){
    p3 <- p3 + xlim(xlims[1], xlims[2])
  }
  if(!is.null(ylims)){
    p3 <- p3 + ylim(ylims[1], ylims[2])
  }
  if(!is.null(col_xy)){
    p3 <- p3 + geom_point(mapping=aes(x=col_xy[1], y=col_xy[2]),
                          pch=8, color=col_color, size=col_size) +
      xlab(NULL) + ylab(NULL)
  }

  ggarrange(p1, p2, p3,
            nrow=1,
            widths=c(.29, .355, .355),
            labels=labels)

}


################################################################################
################################################################################
################################################################################
################################################################################

make_trax <- function(tracks_dir = 'poop_tracks/',
                      result_dir = '',
                      result_suffix = ''){

  (fn <- paste0(result_dir, 'trax', result_suffix, '.rds'))
  message(fn)

  if(file.exists(fn)){
    message('--- File exists! Loading now.')
    load(fn)
  }else{
    message('--- File not found. Creating now... \n')
    trax <- data.frame()
    (lf <- list.files(tracks_dir))
    #trips$colony_size %>% sort %>% table
    #trips$maxdist_linear %>% max
    (kms <- data.frame(km=0:186))
    i=3
    for(i in 1:length(lf)){
      message('Preparing trip stage summary for poop track ',i,' out of ', length(lf),' ...')
      (lfi <- paste0(tracks_dir,lf[i]))
      load(lfi)
      energetics %>% as.data.frame %>% head
      trackexp %>% as.data.frame %>% head

      if(sum(trackexp$ingesta, na.rm=TRUE)>0){
        suppressMessages({
          traxi <-
            trackexp %>%

            # basic stats
            mutate(track_id = i) %>%
            mutate(excreta_ingesta_ratio = energetics$excreta_ingesta_ratio[1]) %>%
            mutate(colony_size = energetics$colony_size[1]) %>%
            mutate(n_trips = energetics$n_trips[1]) %>%
            mutate(kmax = max(km, na.rm=TRUE)) %>%
            mutate(qkm = km/kmax) %>%
            mutate(kmfloor = floor(km)) %>%
            mutate(kmceiling = ceiling(km)) %>%
            #select(-(deployment:waterdist), -(status_last:transit_clock)) %>%

            # inventory gut contents at beginning and end of trip
            mutate(gutn = gut[n()]) %>%
            mutate(gut0 = energetics$gut0) %>%
            mutate(gutcol = gutn - gut0) %>%

            # ingestion / excretion totals
            mutate(tot_ing = sum(ingesta, na.rm=TRUE)) %>%
            mutate(exc_trip = sum(excreta, na.rm=TRUE)) %>%
            mutate(exc0 = gut0 * excreta_ingesta_ratio) %>% # this is included in exc_trip
            mutate(excn = gutn * excreta_ingesta_ratio) %>%
            mutate(tot_exc = exc_trip[1] + excn[1]) %>%
            mutate(sea_exc = exc_trip[1]) %>%
            mutate(col_exc = gutcol * excreta_ingesta_ratio) %>%
            select(-(gutn:gutcol), -(exc_trip:excn), -excreta_ingesta_ratio, -gut) %>%

            # find t and km for ingestion median
            arrange(km) %>%
            #arrange(t) %>%
            mutate(cum_ing = cumsum(ingesta)) %>%
            mutate(qing = cum_ing / tot_ing) %>%
            mutate(km_ing50 = km[which.min(abs(0.50 - qing))]) %>%
            #mutate(t_ing50 = tail(t[qing <= 0.5], 1)) %>%
            #mutate(tdiffs = t - t_ing50) %>%
            #mutate(tdmini = which.min(abs(tdiffs))) %>%
            #mutate(km_ing50 = km[tdmini]) %>%
            #select(-cum_ing, -tdiffs, -tdmini, -t_ing50) %>%

            # define trip stages based on ingestion median
            mutate(stage = ifelse(km < km_ing50, 'inner', 'outer')) %>%
            # anything within 1km of colony is colony
            mutate(stage = ifelse(km < 1, 'colony', stage)) %>%

            # same for excretion
            mutate(cum_exc = cumsum(excreta)) %>%
            mutate(qexc = cum_exc / sea_exc) %>%
            mutate(q01 = km[which.min(abs(0.01 - qexc))]) %>%
            mutate(q05 = km[which.min(abs(0.05 - qexc))]) %>%
            mutate(q10 = km[which.min(abs(0.10 - qexc))]) %>%
            mutate(q20 = km[which.min(abs(0.20 - qexc))]) %>%
            mutate(q25 = km[which.min(abs(0.25 - qexc))]) %>%
            mutate(q50 = km[which.min(abs(0.50 - qexc))]) %>%
            mutate(q75 = km[which.min(abs(0.75 - qexc))]) %>%
            mutate(q80 = km[which.min(abs(0.80 - qexc))]) %>%
            mutate(q90 = km[which.min(abs(0.90 - qexc))]) %>%
            mutate(q95 = km[which.min(abs(0.95 - qexc))]) %>%
            mutate(q99 = km[which.min(abs(0.99 - qexc))]) %>%

            # rearrange columns
            select(colony_size, site, track_id, kmax, n_trips, x, y, status, ms, kmh,
                   t, q, km, kmfloor, kmceiling, qkm,
                   tot_ing, ingesta, qing, km_ing50, qexc,
                   q01, q05, q10, q20, q25, q50, q75, q80, q90, q95, q99,
                   stage, tot_exc, sea_exc, col_exc, excreta)

          traxi %>% as.data.frame %>% tail

          # Summarize by km
          traxkm <-
            traxi %>%
            group_by(colony_size, site, track_id, kmax, n_trips, tot_ing, km_ing50,
                     q01, q05, q10, q20, q25, q50, q75, q80, q90, q95, q99,
                     tot_exc, sea_exc, col_exc, stage, kmfloor, kmceiling) %>%
            summarize(qkm = mean(qkm, na.rm=TRUE),
                      ingesta = sum(ingesta, na.rm=TRUE),
                      excreta = sum(excreta, na.rm=TRUE))
        })
        trax <- rbind(trax, traxkm)
      }
    }
    save(trax, file=fn)
  }
  return(trax)
}


################################################################################
################################################################################
################################################################################
# The big multiplot ==========================================================

  trax_results <- function(trax,
                           col_lm = NULL){ # slope # int
    trax %>% as.data.frame %>% head

    traxtab <-
      trax %>%
      filter(site %in% sites) %>%
      # excretion quantiles
      #group_by(colony_size) %>%
      #mutate(cum_exc = cumsum(excreta)) %>%
      #mutate(q_exc = cum_exc / sea_exc) %>%
      #mutate(q01 = kmfloor[which.min(abs(0.01 - q_exc))]) %>%
      #mutate(q05 = kmfloor[which.min(abs(0.05 - q_exc))]) %>%
      #mutate(q10 = kmfloor[which.min(abs(0.10 - q_exc))]) %>%
      #mutate(q20 = kmfloor[which.min(abs(0.20 - q_exc))]) %>%
      #mutate(q25 = kmfloor[which.min(abs(0.25 - q_exc))]) %>%
      #mutate(q50 = kmfloor[which.min(abs(0.5 - q_exc))]) %>%
      #mutate(q75 = kmfloor[which.min(abs(0.75 - q_exc))]) %>%
      #mutate(q80 = kmfloor[which.min(abs(0.8 - q_exc))]) %>%
      #mutate(q90 = kmfloor[which.min(abs(0.9 - q_exc))]) %>%
      #mutate(q95 = kmfloor[which.min(abs(0.95 - q_exc))]) %>%
      #mutate(q99 = kmfloor[which.min(abs(0.99 - q_exc))]) %>%
      group_by(colony_size, track_id) %>%
      mutate(inner_exc = sum(excreta[stage == 'inner'])) %>%
      mutate(outer_exc = sum(excreta[stage == 'outer'])) %>%
      slice(1) %>%
      group_by(colony_size) %>%
      #group_by(colony_size, q01, q05, q10, q20, q25, q50, q75, q80, q90, q95, q99) %>%
      summarize(n = length(unique(track_id)),
                n_trips = mean(n_trips, na.rm=TRUE),
                kmax_mean = median(kmax, na.rm=TRUE),
                kmax_sd = sd(kmax, na.rm=TRUE),
                km_ing50_mean = mean(km_ing50, na.rm=TRUE),
                km_ing50_median = median(km_ing50, na.rm=TRUE),
                km_ing50_sd = sd(km_ing50, na.rm=TRUE),
                q01 = median(q01),
                q05 = median(q05),
                q10 = median(q10),
                q20 = median(q20),
                q25 = median(q25),
                q50 = median(q50),
                q75 = median(q75),
                q80 = median(q80),
                q90 = median(q90),
                q95 = median(q95),
                q99 = median(q99),
                exc_tot_mean = mean(tot_exc),
                exc_tot_sd = sd(tot_exc),
                inner_raw = mean(inner_exc),
                inner_raw_sd = sd(inner_exc),
                outer_raw = mean(outer_exc),
                outer_raw_sd = sd(outer_exc),
                col_exc_mean = mean(col_exc / tot_exc),
                col_exc_sd = sd(col_exc / tot_exc),
                sea_exc_mean = mean(sea_exc / tot_exc),
                sea_exc_sd = sd(sea_exc / tot_exc),
                inner_tot_mean = mean(inner_exc / tot_exc),
                inner_tot_sd = sd(inner_exc / tot_exc),
                outer_tot_mean = mean(outer_exc / tot_exc),
                outer_tot_sd = sd(outer_exc / tot_exc),
                inner_sea_mean = mean(inner_exc / sea_exc),
                inner_sea_sd = sd(inner_exc / sea_exc),
                outer_sea_mean = mean(outer_exc / sea_exc),
                outer_sea_sd = sd(outer_exc / sea_exc)) %>%
      # corrections
      mutate(q50 = ifelse(q50 > q75, q75, q50)) %>%
      mutate(q50 = ifelse(q50 < q25, q25, q50)) %>%
      #mutate(q99 = ifelse(q99 < q75, q75, q99)) %>%
      #mutate(q95 = ifelse(q95 < q75, q75, q95)) %>%
      #mutate(q90 = ifelse(q90 < q75, q75, q90)) %>%
      #mutate(q80 = ifelse(q80 < q75, q75, q80)) %>%
      # SEs
      mutate(kmax_se = kmax_sd / sqrt(n),
             km_ing50_se = km_ing50_sd / sqrt(n),
             exc_tot_se = exc_tot_sd / sqrt(n),
             col_exc_se = col_exc_sd / sqrt(n),
             sea_exc_se = sea_exc_sd / sqrt(n),
             inner_tot_se = inner_tot_sd / sqrt(n),
             outer_tot_se = outer_tot_sd / sqrt(n),
             inner_sea_se = inner_sea_sd / sqrt(n),
             outer_sea_se = outer_sea_sd / sqrt(n)) %>%
      # colony-wide values
      mutate(col_exc_tot = col_exc_mean * exc_tot_mean * n_trips* colony_size,
             col_exc_tot_se = col_exc_se * exc_tot_mean * n_trips * colony_size) %>%
      # change to metric tons
      mutate(col_exc_tot = col_exc_tot / (1000*1000),
             col_exc_tot_se = col_exc_tot_se / (1000*1000))

    if(!is.null(col_lm)){
      traxtab %>% head
      # lm(log(km_ing50_mean) ~ log(colony_size), data=traxtab) %>% summary
      # (x <- traxtab$colony_size %>% sort)
      # (y <- 0.24178*log(x) + 0.81420)
      # (y <- exp(y))

      (ing50_scale <- traxtab$km_ing50_mean / traxtab$kmax_mean)
      y <- (col_lm[1]*log(traxtab$colony_size)) + col_lm[2]
      y <- exp(y)
      traxtab$km_ing50_mean <- y
      traxtab$kmax_mean <- y/ing50_scale

    }

    traxtab <- traxtab %>%
      # densities
      #mutate(area_inner = (km_ing50_mean^2)*pi - (1^2)*pi) %>%
      #mutate(area_outer = (kmax_mean^2)*pi - (km_ing50_mean^2)*pi) %>%
      mutate(area_inner = (km_ing50_median^2)*pi - (1^2)*pi) %>%
      mutate(area_outer = (kmax_mean^2)*pi - (km_ing50_median^2)*pi) %>%
      mutate(area_all = (kmax_mean^2)*pi - (1^2)*pi) %>%
      mutate(d_all = (sea_exc_mean*exc_tot_mean) / area_all,
             d_all_sd = mean(c(sea_exc_sd, exc_tot_sd)) / area_all) %>%
      mutate(d_inner = inner_raw / area_inner,
             d_inner_sd = inner_raw_sd / area_inner,
             d_outer = outer_raw / area_outer,
             d_outer_sd = outer_raw_sd / area_outer) %>%
      mutate(d_inner_se = d_inner_sd / sqrt(n),
             d_outer_se = d_outer_sd / sqrt(n)) %>%
      mutate(d_all_col = d_all * n_trips * colony_size,
             d_all_col_sd = d_all_sd * n_trips * colony_size) %>%
      mutate(d_inner_col = d_inner * n_trips * colony_size,
             d_inner_sd_col = d_inner_sd * colony_size,
             d_outer_col = d_outer * n_trips * colony_size,
             d_outer_sd_col = d_outer_sd * colony_size) %>%
      mutate(d_all_se = d_all_sd / sqrt(n),
             d_inner_se = d_inner_sd / sqrt(n),
             d_outer_se = d_outer_sd / sqrt(n),
             d_all_col_se = d_all_col_sd / sqrt(n),
             d_inner_se_col = d_inner_sd_col / sqrt(n),
             d_outer_se_col = d_outer_sd_col / sqrt(n)) %>%
      mutate(outer_km = kmax_mean - km_ing50_mean) %>%
      rowwise %>%
      mutate(outer_km_se = mean(c(kmax_se, km_ing50_se))) %>%
      arrange(colony_size) %>%
      mutate(colony_size  = colony_size / 200) %>%
      as.data.frame

    return(traxtab)
  }


################################################################################
################################################################################
################################################################################

  img_xy = c(7, 60)
  img_size = .1
  y_km <- c(0, 80)
  y_trip <- c(70, 220)
  y_share <- c(0.2, 0.8)
  y_conc_col <- c(90, 300, 900, 3000, 9000, 30000)
  y_drop_col <- c(300, 900, 3000, 9000, 30000, 90000)
  y_export <- c(0.1, 0.5)

trax_plot <- function(traxtab,
                      img_xy = c(7, 60),
                      img_size = .3,
                      y_km = c(0, 90),
                      y_trip = c(70, 220),
                      y_share = c(0.2, 0.8),
                      y_conc_bird = c(0, .6),
                      y_conc_col = c(90, 300, 900, 3000, 9000, 40000),
                      y_drop_col = c(0.07, .25, 1, 3, 9, 25),
                      y_export = c(0.0, 0.5)){

    viridis::rocket(n=10, alpha=1)
    #outer_col <- '#841E5AFF'
    outer_col <- 'firebrick'
    outer_col_smooth <- "#841E5A80"
    inner_col <- '#F5936AFF'
    inner_col_smooth <- "#F5936A80"

    #=============================================================================

    imgdf <- data.frame(x=img_xy[1], y = img_xy[2])
    poly_outer <- data.frame(x=c(traxtab$colony_size, rev(traxtab$colony_size)),
                         y=c(traxtab$kmax_mean, rev(traxtab$km_ing50_median)))
    poly_inner <- data.frame(x=c(traxtab$colony_size, rev(traxtab$colony_size)),
                             y=c(rep(1, times=nrow(traxtab)), rev(traxtab$km_ing50_median)))
    poly_col <- data.frame(x=c(traxtab$colony_size, rev(traxtab$colony_size)),
                             y=c(rep(1, times=nrow(traxtab)), rep(0, times=nrow(traxtab))))
    (p_halo <-
        ggplot() +
        geom_image(data=imgdf, mapping=aes(x=x, y=y, image = 'fig_inset.png'),
                   size = img_size)+
        geom_polygon(data=poly_outer, aes(x=x, y=y), col=NA, fill=outer_col, alpha=.6) +
        geom_polygon(data=poly_inner, aes(x=x, y=y), col=NA, fill=inner_col, alpha=.2) +
        geom_polygon(data=poly_col, aes(x=x, y=y), col=NA, fill='darkslategray4', alpha=.4) +
        geom_point(data=traxtab,
                  mapping=aes(x=colony_size,
                              y=kmax_mean),
                  color='grey70', pch=8) +
       geom_linerange(data=traxtab,
                      mapping=aes(x=colony_size,
                                  ymin=kmax_mean - kmax_se, ymax=kmax_mean + kmax_se),
                      color='grey70') +
       geom_linerange(data=traxtab,
                      mapping=aes(x=colony_size,
                                  ymin=km_ing50_median - km_ing50_se,
                                  ymax=km_ing50_median + km_ing50_se),
                      color='black') +
       geom_point(data=traxtab,
                  mapping=aes(x=colony_size,
                              y=km_ing50_median),
                  color='white', size=3) +
       geom_point(data=traxtab,
                  mapping=aes(x=colony_size,
                              y=km_ing50_median),
                  color='black', size=2) +
       scale_x_continuous(trans='log', breaks=c(5, 10, 30, 100, 300, 1000, 3000), limits=c(5, 3000)) +
       scale_y_continuous(limits=y_km) +
       xlab('log Colony pairs x 100') +
       ylab('Colony distance (km)') +
       labs(title = 'Predation halo') +
       theme_minimal()
    )

    #=============================================================================

    traxtab %>% head
    qcol <- 'grey60'
    (p_qs <- ggplot() +
      geom_linerange(data=traxtab, mapping=aes(x=colony_size, ymin=q01, ymax=q99),
                     color=qcol, alpha=1, linewidth=.5) +
      geom_linerange(data=traxtab, mapping=aes(x=colony_size, ymin=q05, ymax=q95),
                     color=qcol, alpha=1, linewidth=1) +
      geom_linerange(data=traxtab, mapping=aes(x=colony_size, ymin=q10, ymax=q90),
                     color=qcol, alpha=1, linewidth=1.5) +
      geom_linerange(data=traxtab, mapping=aes(x=colony_size, ymin=q20, ymax=q80),
                     color=qcol, alpha=1, linewidth=2.5) +
      geom_linerange(data=traxtab, mapping=aes(x=colony_size, ymin=q25, ymax=q75),
                     color=qcol, alpha=1, linewidth=5) +
      geom_linerange(data=traxtab, mapping=aes(x=colony_size, ymin=q50-.1, ymax=q50+.1),
                     color='black', alpha=.8, linewidth=7) +
      #geom_linerange(data=traxtab, mapping=aes(x=colony_size, ymin=km_ing50_median-.1, ymax=km_ing50_median+.1),
      #               color='black', alpha=.5, linewidth=6) +
      geom_point(data=traxtab, mapping=aes(x=colony_size, y=km_ing50_median),
                 color='white', pch=16, size=2.2) +
      geom_point(data=traxtab, mapping=aes(x=colony_size, y=km_ing50_median),
                 color='black', pch=16, size=1.5) +
      geom_point(data=traxtab, mapping=aes(x=colony_size, y=kmax_mean),
                 color='grey70', pch=8, size=2) +
      scale_x_continuous(trans='log', breaks=c(5, 10, 30, 100, 300, 1000, 3000), limits=c(5, 3000)) +
        scale_y_continuous(limits=y_km) +
        xlab('log Colony pairs x 100') +
      ylab('Colony distance (km)') +
      labs(title = 'Offshore guano distribution') +
      theme_minimal())

    #=============================================================================

    (p_drop_trip <-
       ggplot(traxtab, aes(x=colony_size,
                           y=exc_tot_mean,
                           ymin=exc_tot_mean - exc_tot_se, ymax=exc_tot_mean + exc_tot_se)) +
       geom_point() +
       geom_linerange() +
       scale_x_continuous(trans='log', breaks=c(5, 10, 30, 100, 300, 1000, 3000), limits=c(5, 3000)) +
       scale_y_continuous(limits=y_trip) +
       xlab('log Colony pairs x 100') +
       ylab('g') +
       labs(title = expression('Guano production trip'^'-1')) +
       theme_minimal()
    )

    #=============================================================================

    (p_share <-
       ggplot(traxtab, aes(x=colony_size,
                           y=col_exc_mean,
                           ymin=col_exc_mean - col_exc_se, ymax=col_exc_mean + col_exc_se)) +
       geom_point() +
       geom_linerange() +
       scale_x_continuous(trans='log', breaks=c(5, 10, 30, 100, 300, 1000, 3000), limits=c(5, 3000)) +
       scale_y_continuous(limits=y_export) +
       xlab('log Colony pairs x 100') +
       ylab('Proportion') +
       labs(title = 'Share exported to colony') +
       theme_minimal()
    )

    #=============================================================================

    (p_drop_col <-
       ggplot(traxtab,
              aes(x=colony_size,
                  y=col_exc_tot,
                  ymin=col_exc_tot - col_exc_tot_se, ymax=col_exc_tot + col_exc_tot_se)) +
       geom_point() +
       geom_linerange() +
       scale_x_continuous(trans='log',
                          breaks=c(5, 10, 30, 100, 300, 1000, 3000),
                          limits=c(5, 3000)) +
       scale_y_continuous(trans='log',
                          labels=scales::label_comma(),
                          breaks=y_drop_col) +
       xlab('log Colony pairs x 100') +
       ylab('mton') +
       labs(title = expression('Daily guano deposition at colony')) +
       theme_minimal()
    )

    #=============================================================================

    (p_zone <-
       ggplot() +
       geom_point(data=traxtab,
                  mapping = aes(x=colony_size,
                                y=inner_sea_mean),
                  color = inner_col, pch=8) +
       geom_linerange(data=traxtab,
                      mapping = aes(x=colony_size,
                                    ymin=inner_sea_mean - inner_sea_se, ymax=inner_sea_mean + inner_sea_se),
                      color = inner_col, lty=1) +
       geom_smooth(data=traxtab,
                   mapping = aes(x=colony_size,
                                 y=inner_sea_mean),
                   color = inner_col_smooth, se=FALSE, linewidth=.5) +
       geom_point(data=traxtab,
                  mapping = aes(x=colony_size,
                                y=outer_sea_mean),
                  color = outer_col) +
       geom_linerange(data=traxtab,
                      mapping = aes(x=colony_size,
                                    ymin=outer_sea_mean - outer_sea_se, ymax=outer_sea_mean + outer_sea_se),
                      color = outer_col) +
       geom_smooth(data=traxtab,
                   mapping = aes(x=colony_size,
                                 y=outer_sea_mean),
                   color = outer_col_smooth, se=FALSE, linewidth=.5) +
       geom_hline(yintercept = .5, lty=2, alpha=.5) +
       scale_x_continuous(trans='log', breaks=c(5, 10, 30, 100, 300, 1000, 3000), limits=c(5, 3000)) +
       scale_y_continuous(limits=y_share) +
       xlab('log Colony pairs x 100') +
       ylab('Proportion') +
       labs(title = 'Offshore guano placement') +
       theme_minimal()
    )

    #=============================================================================

    (p_d_bird <-
       ggplot() +
       geom_point(data=traxtab,
                  mapping = aes(x=colony_size,
                                y=d_inner),
                  color = inner_col, pch=8) +
       geom_linerange(data=traxtab,
                      mapping = aes(x=colony_size,
                                    ymin=d_inner - d_inner_se,
                                    ymax=d_inner + d_inner_se),
                      color = inner_col, lty=1) +
       geom_smooth(data=traxtab,
                   mapping = aes(x=colony_size,
                                 y=d_inner),
                   color = inner_col_smooth, se=FALSE, linewidth=.5) +
       geom_point(data=traxtab,
                  mapping = aes(x=colony_size,
                                y=d_outer),
                  color = outer_col) +
       geom_linerange(data=traxtab,
                      mapping = aes(x=colony_size,
                                    ymin=d_outer - d_outer_se, ymax=d_outer + d_outer_se),
                      color = outer_col) +
       geom_smooth(data=traxtab,
                   mapping = aes(x=colony_size,
                                 y=d_outer),
                   color = outer_col_smooth, se=FALSE, linewidth=.5) +
       scale_x_continuous(trans='log', breaks=c(5, 10, 30, 100, 300, 1000, 3000), limits=c(5, 3000)) +
       scale_y_continuous(limits=y_conc_bird) +
       xlab('log Colony pairs x 100') +
       ylab(expression(paste("g km"^"-2"))) +
       labs(title = expression(paste("Offshore guano concentration bird"^"-1"))) +
       theme_minimal()
    )

    #=============================================================================

    (p_d_col <-
       ggplot() +
       geom_point(data=traxtab,
                  mapping = aes(x=colony_size,
                                y=d_inner_col),
                  color = inner_col, pch=8) +
       geom_linerange(data=traxtab,
                      mapping = aes(x=colony_size,
                                    ymin=d_inner_col - d_inner_se_col,
                                    ymax=d_inner_col + d_inner_se_col),
                      color = inner_col, lty=1) +
       geom_smooth(data=traxtab,
                   mapping = aes(x=colony_size,
                                 y=d_inner_col),
                   color = inner_col_smooth, se=FALSE, linewidth=.5) +
       geom_point(data=traxtab,
                  mapping = aes(x=colony_size,
                                y=d_outer_col),
                  color = outer_col) +
       geom_linerange(data=traxtab,
                      mapping = aes(x=colony_size,
                                    ymin=d_outer_col - d_outer_se_col, ymax=d_outer_col + d_outer_se_col),
                      color = outer_col) +
       geom_smooth(data=traxtab,
                   mapping = aes(x=colony_size,
                                 y=d_outer_col),
                   color = outer_col_smooth, se=FALSE, linewidth=.5) +
       scale_x_continuous(trans='log',
                          breaks=c(5, 10, 30, 100, 300, 1000, 3000), limits=c(5, 3000)) +
       scale_y_continuous(trans='log',
                          labels = scales::label_comma(),
                          breaks=y_conc_col) +
       xlab('log Colony pairs x 100') +
       ylab(expression(paste("log g km"^"-2"))) +
       labs(title = expression(paste("Offshore guano concentration colony"^"-1"))) +
       theme_minimal()
    )

    ggarrange(p_halo, p_drop_trip,
              p_qs, p_zone,
              p_d_bird, p_d_col,
              p_share, p_drop_col,
              nrow=4, ncol=2, labels='auto')
  }


################################################################################
################################################################################
################################################################################
################################################################################
# Trip quantiles

trip_quantiles <- function(tracks_dir = 'poop_tracks2/',
                           result_dir = '',
                           result_suffix = ''){

  if(FALSE){
    tracks_dir <- 'poop_tracks2/'
    result_suffix <- ''
    result_dir <- ''
  }

  (fn <- paste0(result_dir, 'trip_quantiles', result_suffix, '.rds'))
  message(fn)

  if(file.exists(fn)){
    message('--- File exists! Loading now.')
    load(fn)
  }else{
    message('--- File not found. Creating now... \n')

    (lf <- list.files(tracks_dir))
    i=1
    cdfs <- data.frame()
    for(i in 1:length(lf)){
      message(i,' out of ', length(lf),' ...')
      (lfi <- paste0(tracks_dir,lf[i]))
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
                    sing = sum(ingesta, na.rm=TRUE),
                    sexc = sum(excreta, na.rm=TRUE),
                    ming = mean(ingesta, na.rm=TRUE),
                    mexc = mean(excreta, na.rm=TRUE),
                    ping = mean(ping),
                    pexc = mean(pexc)) %>%
          as.data.frame #%>% head
      })
      cdfi %>% head
      cdfi$ptr
      cdfi$ing_active
      #ggplot(cdfi, aes(x=ptr, y=pkm)) + geom_point()
      #ggplot(cdfi %>% filter(ing_active == TRUE), aes(x=ptr, y=pkm)) + geom_point()
      cdfs <- rbind(cdfs, cdfi)
    }
    save(cdfs, file=fn)
  }

  return(cdfs)
}

################################################################################
################################################################################
################################################################################
################################################################################
# Quantile comparison plot

quantile_plot <- function(cdfs,
                          quantile_lims = c(0, 1),
                          titles = c('Observed results', '', ''),
                          subtitles = c('Ingestion', 'Excretion', 'Ingestion-excretion offset'),
                          labels='auto'){

  cdfs %>% head
  cdfsumm <-
    cdfs %>%
    group_by(colony_size, ptr) %>%
    summarize(ing_active = length(which(ing_active)) / n(),
              exc_active = length(which(exc_active)) / n(),
              sing = mean(sing, na.rm=TRUE),
              sexc = mean(sexc, na.rm=TRUE),
              ming = mean(ming, na.rm=TRUE),
              mexc = mean(mexc, na.rm=TRUE),
              pkm = mean(pkm, na.rm=TRUE),
              ping = mean(ping, na.rm=TRUE),
              pexc = mean(pexc, na.rm=TRUE))
  cdfsumm %>% head

  # Median ingestion and excretion quantile-distances
  ingexc <-
    cdfsumm %>%
    group_by(colony_size) %>%
    summarize(pkm_ing50 = pkm[which(ping >= .5)[1]],
              pkm_exc50 = pkm[which(pexc >= .5)[1]]) %>%
    mutate(pdiff = pkm_exc50 - pkm_ing50)

  ingexc %>% as.data.frame
  ingexc$pdiff %>% mean(na.rm=TRUE)

  p <-
    ggarrange(
      ggplot(cdfsumm %>% mutate(colony_size = colony_size / 100),
             aes(x=ptr, y=pkm, group=colony_size,
                 color=ming)) +
        geom_path(lwd=.8, alpha=.7) +
        #geom_hline(data=ing_medians, mapping=aes(yintercept=pkm),
        #           alpha=.2, color='black') +
        viridis::scale_color_viridis(direction=-1) +
        xlab('Quantile of trip duration') +
        ylab('Quantile distance from colony') +
        theme_light() +
        labs(title = titles[1],
             subtitle = subtitles[1],
             color='Mean\ningestion') +
        theme(plot.title = element_text(face = "bold")),

      ggplot(cdfsumm %>% mutate(colony_size = colony_size / 100),
             aes(x=ptr, y=pkm, group=colony_size,
                 color=mexc)) +
        geom_path(lwd=.8, alpha=.7) +
        #geom_hline(data=exc_medians, mapping=aes(yintercept=pkm),
        #           alpha=.2, color='black') +
        viridis::scale_color_viridis(direction=-1) +
        xlab('Quantile of trip duration') +
        ylab('Quantile distance from colony') +
        theme_light() +
        labs(title = titles[2],
             subtitle = subtitles[2],
             color='Mean\nexcretion') +
        theme(plot.title = element_text(face = "bold")),

      ggplot(ingexc %>% mutate(colony_size = colony_size / 100),
             aes(x=pkm_ing50,
                 y=pkm_exc50,
                 size=colony_size)) +
        geom_point(alpha=.6) +
        geom_abline(lty=2) +
        xlim(quantile_lims[1], quantile_lims[2]) +
        ylim(quantile_lims[1], quantile_lims[2]) +
        theme_light() +
        xlab('Quantile-distance of ingestion midpoint') +
        ylab('Quantile-distance of excretion midpoint') +
        labs(title = titles[3],
             subtitle = subtitles[3],
             size = 'Colony size\n(hundreds\nof pairs)') +
        theme(plot.title = element_text(face = "bold")),

      nrow=1,
      labels=labels
    )

  return(list(cdfsumm = cdfsumm,
              plot = p,
              ingexc = ingexc))
}

################################################################################
################################################################################
################################################################################
################################################################################

sample_trips <- function(trips){
  # get lay of land
  trips %>%
    mutate(slide_id = 1:n()) %>%
    group_by(colony_size) %>%
    mutate(n=n()) %>%
    mutate(colcat = case_when(n < 25 ~ '0-25',
                              n >= 25 & n < 100 ~ '25-100',
                              .default = '>100')) %>%
    group_by(colcat, colony_size) %>%
    tally() %>%
    as.data.frame

  subtrips <-
    trips %>%
    mutate(slide_id = 1:n()) %>%
    group_by(colony_size) %>%
    mutate(n=n()) %>%
    mutate(colcat = case_when(n < 25 ~ '0-25',
                              n >= 25 & n < 100 ~ '25-100',
                              .default = '>100')) %>%
    ungroup %>%
    group_by(colcat) %>%
    group_split

  sliced_trips <-
    rbind(subtrips[[1]],
          subtrips[[2]] %>% group_by(colony_size) %>% slice_sample(n=25),
          subtrips[[3]] %>% group_by(colony_size) %>% slice_sample(prop=.25))

  return(sliced_trips)
}
