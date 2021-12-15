####################################################################################################################################################################
#### TEST THE METRICS ON SIMULATED PHENOLOGY OF BIRTHS #############################################################################################################
####################################################################################################################################################################

# This code is designed to test all the functions coding for the different phenology metrics (timing, synchrony, rhythmicity, regularity) individually on simulated
# phenology of births. Each metric is assorted with several outputs.
# A generic function also returns the outputs used in the associated article, namely one output for each metric.

# SOURCE R CODES ###################################################################################################################################################

source("simulation_birth_phenology_norm.R")
source("phenology_metrics.R")
source("phenology_metrics_top4.R")

# SIMULATE A PHENOLOGY OF BIRTHS ###################################################################################################################################

## Simulate a phenology of births over 10 years of 365 days, characterized by one birth peak, a variable mean birth date and a variable standard deviation of the
# distribution of births
birth_pattern <- simulation.pattern(format="row", graph=T, graph_format="ind",
                                    births_distrib="peaks", nb_tu_per_cycle=365, nb_draws_per_tu=10, nb_cycles=10,
                                    delta_prob=0, delta_mean=30, delta_sd=10, variability=0,
                                    nb_peaks=1, probs=1, mean_peaks=182, sd_peaks=30)

## Extract the first year of the phenology of births
birth_pattern_1 <- birth_pattern$data_supracyc[birth_pattern$data_supracyc$cycle%in%1, c("tu", "nb_births")]

## Extract the two first years of the phenology of births
birth_pattern_2 <- birth_pattern$data_supracyc[birth_pattern$data_supracyc$cycle%in%c(1, 2), c("cycle", "tu", "nb_births")]

## Extract the complete phenology of birth
birth_pattern_3 <- birth_pattern$data_supracyc

# ARGUMENTS FOR THE METRICS ########################################################################################################################################

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

# TEST THE METRICS INDIVIDUALLY ON THE SIMULATED PHENOLOGY OF BIRTHS (SEVERAL OUTPUTS PER METRIC) ##################################################################

## When one year of data is available
pheno.mean.med.skew.sd(pattern=birth_pattern_1, graph=opt.graph)

pheno.var.varcor(pattern=birth_pattern_1, interval=opt.interval, corr=opt.corr)

pheno.bgper.per(pattern=birth_pattern_1, graph=opt.graph)

pheno.nbtu(pattern=birth_pattern_1, percent_min=opt.percent_min, consecutive=opt.consecutive, graph=opt.graph)

pheno.pielou(pattern=birth_pattern_1)

pheno.bgthper.interq(pattern=birth_pattern_1, first_quantile=opt.first_quantile, last_quantile=opt.last_quantile, percent=NA, percent_min=NA, percent_min_end=NA,
                     graph=opt.graph)

pheno.bgthper.interq(pattern=birth_pattern_1, first_quantile=NA, last_quantile=NA, percent=opt.percent, percent_min=NA, percent_min_end=NA, graph=opt.graph)

pheno.bgthper.interq(pattern=birth_pattern_1, first_quantile=NA, last_quantile=NA, percent=NA, percent_min=opt.percent_min, percent_min_end=NA, graph=opt.graph)

pheno.bgthper.interq(pattern=birth_pattern_1, first_quantile=NA, last_quantile=NA, percent=NA, percent_min=NA, percent_min_end=opt.percent_min_end, graph=opt.graph)

pheno.centre(pattern=birth_pattern_1, graph=opt.graph)

pheno.mode(pattern=birth_pattern_1, graph=opt.graph)

pheno.propmed.propmode(pattern=birth_pattern_1, reference="median", period=opt.period, graph=opt.graph)

pheno.propmed.propmode(pattern=birth_pattern_1, reference="mode", period=opt.period, graph=opt.graph)

pheno.perhdr(pattern=birth_pattern_1, percent=opt.percent, graph=opt.graph)

pheno.maxprop.minprop(pattern=birth_pattern_1, period=opt.period, object="max", graph=opt.graph)

pheno.maxprop.minprop(pattern=birth_pattern_1, period=opt.period, object="min", graph=opt.graph)

pheno.diffmima(pattern=birth_pattern_1, period=opt.period, graph=opt.graph)

pheno.meanvl.meanvo(pattern=birth_pattern_1, graph=opt.graph)

pheno.skinner(pattern=birth_pattern_1, period=opt.period, percent=opt.percent, graph=opt.graph)

pheno.minper(pattern=birth_pattern_1, percent=opt.percent, graph=opt.graph)

pheno.zerbe(pattern=birth_pattern_1, percent=opt.percent, graph=opt.graph)

pheno.rutberg(pattern=birth_pattern_1, percent=opt.percent, consecutive=opt.consecutive, graph=opt.graph)

pheno.medprob.sdprob(pattern=birth_pattern_1, graph=opt.graph)

pheno.pergau(pattern=birth_pattern_1, graph=opt.graph)

pheno.rayleigh(pattern=birth_pattern_1, graph=opt.graph)

pheno.kolmogau.kolmouni(pattern=birth_pattern_1, ref="uniform", origin=opt.origin, graph=opt.graph)

pheno.kolmogau.kolmouni(pattern=birth_pattern_1, ref="gaussian", origin=opt.origin, graph=opt.graph)

pheno.comppeaksig.peaksig(pattern=birth_pattern_1, nb_cycles=1, CI=opt.CI, resolution=opt.resolution, graph=opt.graph)

## When two years of data are available

pheno.compmean(pattern=birth_pattern_2, CI=opt.CI, graph=opt.graph)

pheno.watson(pattern=birth_pattern_2, graph=opt.graph)

pheno.mood(pattern=birth_pattern_2, graph=opt.graph)

pheno.diffbgper.diffmed.diffpeak.diffperiod(pattern=birth_pattern_2, graph=opt.graph)

pheno.kolmomult(pattern=birth_pattern_2, graph=opt.graph)

## When more than two years of data are available

pheno.meanmult.permean(pattern=birth_pattern_3, graph=opt.graph)

pheno.cmano(pattern=birth_pattern_3, post_hoc=opt.post_hoc, graph=opt.graph)

pheno.diffmean.meanlin.varlin(pattern=birth_pattern_3, graph=opt.graph)

pheno.bart(pattern=birth_pattern_3, graph=opt.graph)

pheno.slpcomp(pattern=birth_pattern_3, transformation=opt.transformation, graph=opt.graph)

pheno.khi2(pattern=birth_pattern_3, period=opt.period, graph=opt.graph)

pheno.diffslin(pattern=birth_pattern_3, graph=opt.graph)

# TEST THE METRICS ON THE SIMULATED PHENOLOGY OF BIRTHS USING THEGENRIC FUNCTION USED IN THE ARTICLE (ONE OUTPUT PER METRIC) #######################################

metrics.test(pattern=birth_pattern_3, graph=opt.graph, interval=opt.interval, percent_min=opt.percent_min, consecutive=opt.consecutive,
             first_quantile=opt.first_quantile, last_quantile=opt.last_quantile, percent=opt.percent, nsim=opt.nsim,confidence=opt.confidence, CI=opt.CI,
             nboots=opt.nboots, transformation=opt.transformation, post_hoc=opt.post_hoc, resolution=opt.resolution, period=opt.period)

metrics.test.best4(pattern=birth_pattern_3, graph=opt.graph)

####################################################################################################################################################################
####################################################################################################################################################################
