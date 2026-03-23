################################################################################
################################################################################
# The nutrient dimension of predation halos
# Eric K. Ezell, ekezell@sewanee.edu
#
# Step 2 Parameters
################################################################################
################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
################################################################################

params <-
  list(
    murre_mass = 888,
    # g
    # COMU - Hilton et al. 2000

    # Ingestion-Excretion ======================================================
    excretion_rate = 0.00013,
    # % wet meal mass excreted per second
    # COMU (middle of the two values,  0.011 and 0.015, Hilton et al 2000)
    ingesta_transit_time = 70*60,
    # sec to first poop once foraging begins
    # COMU - Hilton et al. 2000
    max_ingestion_rate = 0.037 ,
    # grams per second
    # COMU - Hilton et al 2000, pg 24 of 24
    max_feeding_time = 197*60,
    # second
    # COMU - mean of 186, 318, 89; Hilton et al 2000, pg 24 of PDF
    gut_capacity = 192,
    # g or mL,
    # COMU - Hilton et al 2000 pg 24 of pdf
    assimilation_efficiency = 0.79,
    # COMU - Hilton et al 2000b, also used in Roth et al 2008
    excreta_ingesta_ratio = .427,
    # Proportion of ingesta that is excreted
    # COMU
    # Hilton et al 2000b fed murres 10.88% of their body weight in fish
    # Based on figure 3, the murre weight is 990g -- so they were fed 107g wet mass
    # Using an equation from Larimer (1992), dry = 0.309*wet - 0.286, which is for sand lance, that would be 32.77g dry mass
    # Based on figure 2, murres excreted 14g dry mass -- so their excreta was equal to 14/32.77 = 42.7% of their ingesta by dry weight
    E_density_diet = 6,
    # kJ/g; middle of 4-8 kJ from Hilton et al 2000

    # Trip / diving ============================================================
    trip_hours = 18,
    # hours
    # COMU - based on the 18 hours per pair in Hilton et al 2000
    dive_time = 77,
    dive_time_sd = 28.6,
    # sec
    # COMU - Evans et al 2013, Table 3
    dives_per_bout = 7,
    dives_per_bout_sd = 5.74,
    # COMU - Evans et al 2013, Table 3 6.9
    post_dive_interval = 42,
    post_dive_interval_sd = 30,
    # s
    # COMU - Evans et al 2013
    bout_duration = 833,
    bout_duration_sd = 712,
    # s
    # COMU - 803 is observed from Evans et al; 833 makes sense based on params above
    bout_interval = 250,
    # s
    # COMU - Evans et al 2013
    percent_feeding = 0.50,
    # percent of "sit" time spent feeding underwater;
    # COMU - Evans et al. 2013 had 50%. (Note: Hilton used 0.63)
    percent_inactive = 0.31
    # percent of "sit" non-feed time inactive (legs tucked)
    # TBMU - Elliott & Gaston 2014 - reported fraction for Sep - October
    #  less diel variation (most within 20% - 40%) then later months (5% - 50%) but still substantial
    # simplify it by taking mean across hours of day: mean(c(47, 48, 46, 22, 15, 18, 20, 23, 24, 27, 38, 45)) = .31
    )

params



#===============================================================================
#===============================================================================
#===============================================================================
# Re-run gut curve based on parameters

if(!file.exists('gutcurve.rds')){
  message('Reproducing the gut contents decay curve ...')
  # Create a decay curve for gut contents
  t = 1:86400 # one day
  (er <- params$excretion_rate)
  (guts <- params$gut_capacity)
  for(i in t){
    last_gut <- guts[length(guts)]
    guti <- last_gut - (last_gut*er)
    #message(guti)
    guts <- c(guts, guti)
  }
  plot(guts)
  save(gutcurve, file='gutcurve.rds')
}




