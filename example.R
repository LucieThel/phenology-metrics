####################################################################################################################################################################
#### TEST THE METRICS ON SIMULATED PHENOLOGY OF BIRTHS #############################################################################################################
####################################################################################################################################################################

# This code is designed to test all the functions coding the different phenology metrics (timing, synchrony, rhythmicity, regularity) on simulated phenology of
# births.

# SOURCE R CODES ###################################################################################################################################################

# source("simulation_birth_phenology.R")
# source("phenology_metrics.R")

# SIMULATE A PHENOLOGY OF BIRTHS ###################################################################################################################################

## Simulate a phenology of births over 10 years of 365 days, characterized by one birth peak, a variable mean birth date and a variable standard deviation of the
# distribution of births
birth_pattern <- simulation.pattern(format="row", graph=T, graph_format="ind",
                                    births_distrib="peaks", nb_tu_per_cycle=365, nb_draws_per_tu=10, nb_cycles=10,
                                    delta_prob=0, delta_mean=30, delta_sd=10, delta_off_period_start=0, delta_off_period_end=0, variability=0,
                                    nb_peaks=1, probs=1, mean_peaks=182, sd_peaks=30,
                                    off_period_start=0, off_period_end=0)

## Extract the first year of the phenology of births
birth_pattern_1 <- birth_pattern$data_supracyc[birth_pattern$data_supracyc$cycle%in%1, c("tu", "nb_births")]

## Extract the two first years of the phenology of births
birth_pattern_2 <- birth_pattern$data_supracyc[birth_pattern$data_supracyc$cycle%in%c(1, 2), c("cycle", "tu", "nb_births")]

## Extract the complete phenology of birth
birth_pattern_3 <- birth_pattern$data_supracyc

# TEST THE METRICS ON THE SIMULATED PHENOLOGY OF BIRTHS ############################################################################################################

## Arguments for the metrics
opt.graph <- T
opt.interval <- "check" 
opt.corr <- T
opt.percent_min <- 1
opt.percent_min_end <- 1
opt.consecutive <- F
opt.first_quantile <- 25
opt.last_quantile <- 75
opt.percent <- 80
opt.nsim <- 1000
opt.confidence <- F
opt.CI <- 95
opt.nboots <- 1000
opt.origin <- "corr"
opt.transformation <- "none"
opt.post_hoc <- F
opt.resolution <- 0.1
opt.period <- 90 # equivalent to three months

## Apply the functions associated to the metrics extracted from the literature about phenology of births in large ungulates

# When one year of data is available
pheno.findlaylambin(pattern=birth_pattern_1, graph=opt.graph)

pheno.johnson(pattern=birth_pattern_1, interval=opt.interval, corr=opt.corr)

pheno.bunnellyomtov(pattern=birth_pattern_1, graph=opt.graph)

pheno.moe1(pattern=birth_pattern_1, percent_min=opt.percent_min, consecutive=opt.consecutive, graph=opt.graph)

pheno.sinclair(pattern=birth_pattern_1)

pheno.gaillardmajluf(pattern=birth_pattern_1, first_quantile=opt.first_quantile, last_quantile=opt.last_quantile, percent=NA, percent_min=NA, percent_min_end=NA,
                     graph=opt.graph)

pheno.gaillardmajluf(pattern=birth_pattern_1, first_quantile=NA, last_quantile=NA, percent=opt.percent, percent_min=NA, percent_min_end=NA, graph=opt.graph)

pheno.gaillardmajluf(pattern=birth_pattern_1, first_quantile=NA, last_quantile=NA, percent=NA, percent_min=opt.percent_min, percent_min_end=NA, graph=opt.graph)

pheno.gaillardmajluf(pattern=birth_pattern_1, first_quantile=NA, last_quantile=NA, percent=NA, percent_min=NA, percent_min_end=opt.percent_min_end, graph=opt.graph)

pheno.sigouin(pattern=birth_pattern_1, graph=opt.graph)

pheno.jemison(pattern=birth_pattern_1, graph=opt.graph)

pheno.adamsjemison(pattern=birth_pattern_1, reference="median", period=opt.period, graph=opt.graph)

pheno.adamsjemison(pattern=birth_pattern_1, reference="mode", period=opt.period, graph=opt.graph)

pheno.calabrese(pattern=birth_pattern_1, percent=opt.percent, graph=opt.graph)

pheno.owensmith1(pattern=birth_pattern_1, period=opt.period, object="max", graph=opt.graph)

pheno.owensmith1(pattern=birth_pattern_1, period=opt.period, object="min", graph=opt.graph)

pheno.owensmith2(pattern=birth_pattern_1, period=opt.period, graph=opt.graph)

pheno.campos(pattern=birth_pattern_1, graph=opt.graph)

pheno.skinner(pattern=birth_pattern_1, period=opt.period, percent=opt.percent, graph=opt.graph)

pheno.meng(pattern=birth_pattern_1, percent=opt.percent, graph=opt.graph)

pheno.zerbe(pattern=birth_pattern_1, percent=opt.percent, graph=opt.graph)

pheno.rutberg(pattern=birth_pattern_1, percent=opt.percent, consecutive=opt.consecutive, graph=opt.graph)

pheno.caughley(pattern=birth_pattern_1, graph=opt.graph)

pheno.paoli(pattern=birth_pattern_1, graph=opt.graph)

pheno.dibitetti(pattern=birth_pattern_1, graph=opt.graph)

pheno.riedmanmeng(pattern=birth_pattern_1, ref="uniform", origin=opt.origin, graph=opt.graph)

pheno.riedmanmeng(pattern=birth_pattern_1, ref="gaussian", origin=opt.origin, graph=opt.graph)

pheno.moe2(pattern=birth_pattern_1, nb_cycles=1, CI=opt.CI, resolution=opt.resolution, graph=opt.graph)

# When two years of data are available

pheno.whiting(pattern=birth_pattern_2, CI=opt.CI, graph=opt.graph)

pheno.pare(pattern=birth_pattern_2, graph=opt.graph)

pheno.berger(pattern=birth_pattern_2, graph=opt.graph)

pheno.diff(pattern=birth_pattern_2, graph=opt.graph)

pheno.schaik(pattern=birth_pattern_2, graph=opt.graph)

# When more than two years of data are available

pheno.millar(pattern=birth_pattern_3, graph=opt.graph)

pheno.linnell(pattern=birth_pattern_3, post_hoc=opt.post_hoc, graph=opt.graph)

pheno.loe(pattern=birth_pattern_3, graph=opt.graph)

pheno.hass(pattern=birth_pattern_3, graph=opt.graph)

pheno.bowyer(pattern=birth_pattern_3, transformation=opt.transformation, graph=opt.graph)

pheno.adams(pattern=birth_pattern_3, period=opt.period, graph=opt.graph)

pheno.paoli2(pattern=birth_pattern_3, graph=opt.graph)

####################################################################################################################################################################
####################################################################################################################################################################
