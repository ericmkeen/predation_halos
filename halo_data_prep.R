################################################################################
################################################################################
# The nutrient dimension of predation halos
# Eric K. Ezell, ekezell@sewanee.edu
#
# Step 1 Data prep
################################################################################
################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
################################################################################
################################################################################
# GPS track data from Allison Patterson et al 2022

#===============================================================================
# Trip summary key

suppressMessages({
  (trips <- readr::read_delim('tripSummary.txt', delim=';'))
  trips %>% as.data.frame %>% head
  trips %>% nrow
  trips_all <- trips
})

if(FALSE){
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
}

#===============================================================================
# GPS tracks

if(file.exists('tracks_all.rds')){
  load('tracks_all.rds')
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
    mutate(kmh = ifelse(kmh > 108, NA, kmh)) %>%
    # after Patterson 2022 thesis, PDF page 58 -- used 30 m/s as cutoff, translates to 108 kmh
    mutate(status = ifelse(kmh > 7.2, 'fly','sit')) %>%
    # after Patterson 2022 thesis, PDF page 58 -- used 2 m/s as cutoff; translates to 7.2 kmh
    #mutate(status = ifelse(kmh > 10, 'fly','sit')) %>%
    mutate(prev_status = lag(status),
           next_status = lead(status)) %>%
    as.data.frame
  
  tracks %>% head
  length(which(is.na(tracks$kmh))) / nrow(tracks)
  
  # Save dataset
  save(tracks, file='tracks_all.rds')
}


################################################################################
################################################################################
# Filter based on sample size

if(file.exists('trips_used.rds')){
  load('trips_used.rds')
}else{
  # Assess sample size of each trip's track
  mr <- data.frame()
  trips_to_process <- 1:nrow(trips)
  i=1
  for(i in trips_to_process){
    message(i)
    (tripi <- trips[i, ]) %>% as.data.frame
    (tripid <- paste0(tripi$site,'_',tripi$year,'_',
                      tripi$bird_id,'_',tripi$tripID,
                      tripi$deployment))
    (tracki <- 
        tracks %>% 
        mutate(i = 1:n()) %>% 
        filter(bird_id == tripi$bird_id, tripID == tripi$tripID, deployment == tripi$deployment) %>% 
        mutate(maxkm = max(coldist, na.rm=TRUE)) %>% 
        mutate(qkm = coldist / maxkm) %>% 
        mutate(trip_sec = cumsum(secs)) %>% 
        mutate(trip_dur = max(trip_sec, na.rm=TRUE)) %>% 
        mutate(qt = trip_sec / trip_dur))
    
    #tracki$coldist %>% plot
    #ggplot(tracki, aes(x=qt, y=qkm)) + geom_point()
    #(tracki$coldist / max(tracki$coldist, na.rm=TRUE)) %>% plot
    
    tracki #%>% head
    nrow(tracki)
    mri <- data.frame(tripi,
                      tripid,
                      track_n=nrow(tracki), 
                      max_gap = max(tracki$secs, na.rm=TRUE),
                      t2k3 = tracki %>% filter(qt <= .2, qkm <= .3) %>% nrow(),
                      t8k3 = tracki %>% filter(qt >= .8, qkm <= .3) %>% nrow()
                      )
    mri
    mr <- rbind(mr, mri)
  }
  
  mr %>% head
  
  # Try filtering
  mrf <- 
    mr %>% 
    filter(track_n >= 50) %>% 
    filter(t2k3 > 0) %>% 
    filter(t8k3 > 0)
  
  mrf %>% nrow
  trips <- mrf
  save(trips, file='trips_used.rds')
}


#===============================================================================

if(file.exists('tracks_used.rds')){
  load('tracks_used.rds')
}else{
  # Assess sample size of each trip's track
  mr <- data.frame()
  trips_to_process <- 1:nrow(trips)
  i=1
  for(i in trips_to_process){
    message(i)
    (tripi <- trips[i, ]) %>% as.data.frame
    (tracki <- 
        tracks %>% 
        filter(bird_id == tripi$bird_id, tripID == tripi$tripID, deployment == tripi$deployment))
    tracki #%>% head
    nrow(tracki)
    mr <- rbind(mr, tracki)
  }
  
  mr %>% nrow
  tracks <- mr
  save(tracks, file='tracks_used.rds')
}


trips %>% nrow
trips$bird_id %>% unique %>% length

# Review sample size
trips %>% 
  group_by(colony_size) %>% 
  tally() %>% 
  as.data.frame





################################################################################
################################################################################
# exploratory plots

if(FALSE){
  
  tracks %>% 
    group_by(bird_id, tripID, deployment) %>% 
    tally() %>% 
    arrange(n) %>% 
    pull(n) %>% 
    tail(3363)
  
  #(tracks$t2 - tracks$t) %>% median(na.rm=TRUE)
  tracks$kmh %>% hist
  tracks %>% filter(kmh > 15) %>% pull(kmh) %>% mean
  tracks %>% filter(kmh > 15) %>% pull(kmh) %>% median
  tracks %>% filter(kmh > 15) %>% pull(kmh) %>% sd
  
  ggplot(tracks, aes(x=kmh, fill=status)) + geom_histogram()
  ggplot(tracks %>% filter(kmh < 15), aes(x=kmh, fill=status)) + geom_histogram() +
    geom_vline(xintercept=7.2, lty=3)
  
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
  
}
