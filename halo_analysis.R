################################################################################
################################################################################
# The nutrient dimension of predation halos
# Eric K. Ezell, ekezell@sewanee.edu
#
# Step 3 Main text analysis
################################################################################
################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('halo_functions.R')
source('halo_data_prep.R')
source('halo_parameters.R')

################################################################################
################################################################################

# Test it
trip_model(trips_to_process = 1,
           params,
           trips,
           tracks,
           tracks_dir = 'poop_tracks/',
           diagnostics_dir = 'diagnostics/',
           maps_dir = 'track_maps/',
           redo_analysis = TRUE)

################################################################################
# Loop settings for horizontal computing in Sewanee Ocean Lab

# default
nrow(trips)
i_begin <- 1
i_end <- nrow(trips)

if(Sys.info()[7] == "ezelllab"){
  i_begin <- 1
  i_end <- 779
}
if(Sys.info()[7] == "ezelllab2"){
  i_begin <- 780
  i_end <- 1557
}
if(Sys.info()[7] == "ezelllab3"){
  i_begin <- 1558
  i_end <- 2337
}

trips_to_process <- i_begin:i_end
(i=sample(1:2337, size=1))


################################################################################
# Run

trip_model(trips_to_process,
           params,
           trips,
           tracks,
           tracks_dir = 'poop_tracks/',
           diagnostics_dir = 'diagnostics/',
           maps_dir = 'track_maps/',
           redo_analysis = TRUE)

