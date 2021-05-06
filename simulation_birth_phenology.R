####################################################################################################################################################################
#### FUNCTION TO SIMULATE BIRTH PHENOLOGY ##########################################################################################################################
####################################################################################################################################################################

# DOCUMENTATION ####################################################################################################################################################

### DESCRIPTION
# This function generates birth patterns of different shapes (uniform, random, with over dispersion of zeros, unimodal, multimodal, ...) thanks to several
# parameters described in "Arguments" section bellow.

### ARGUMENTS
## General parameters:
# - births_distrib : type of the distribution; "uniform" = births are drawn from a uniform distribution, "random" = births are random drawn, "negbinomial" = births
# are drawn from a negative binomial distribution, "peaks" = births are drawn from a normal law (default="uniform")
# - format : format of the output when more than one year is generated, "row" = by lines, "col" = by columns (default="col")
# - graph : plot or not the graph associated to the simulations (default=FALSE)
# - graph_format : "ind" = represent the output year by year on separate graphs, "sup" = superimpose all the years on the same graph (default="ind")
# - nb_cycles : number of cycles (e.g. years) among which simulating birth patterns (default=3)
# - nb_draws_per_tu : number of females that could give birth every day (default=10). Warning: it is better to keep this value at default to maintain similar
# number of potential births between the different distributions available.
# - nb_tu_per_cycle : number of time unit (e.g. days) among which distributing the births (default=365)

## When birth_distrib = "peaks":
# - mean_peaks : mean(s) of the normal distribution(s) (defaut=c(73, 146, 219, 292))
# - nb_peaks : number of peaks in one cycle (defaut=4)
# - probs : probabilities indicating in which probability density to draw the event "a female gives birth/a female does not give birth". The last value is
# necessarily 1. It corresponds to the cumulative proportion of births happening during each peak (default=c(0.2, 0.5, 0.7, 1))
# - sd_peaks : standard deviation of the normale distribution(s) (defaut=c(15, 15, 15, 15))

## When birth_distrib = "negbinomial":
# - off_period_end : end of the off period (when the number of births is null) (defaut=0)
# - off_period_start : start of the off period (when the number of births is null) (defaut=0)
# If no off period is wanted, the arguments off_period_start and off_period_end should be set to 0 both.

## Consistency across cycles:
# It is possible to modify the level of consistency between the cycles for "peaks" with one peak, "peaks" with two peaks and "negbinomial" distributions:
# - "peaks" with one peak: variation of the mean and/or sd of the distribution of births
# - "peaks" with two peaks: variation of the probability to give birth in the first or the second peak
# - "negbinomial" with one off period: variation of the start and/or the end of the off period
# For all other distributions ("uniform", "random", "peaks" with more than two peaks or "negbinomial" with more than one off period), it is not possible to change
# the level of consistency between the cycles (defaut=parameters of the distribution stay the same for all cycles, no supplementary parameter to provide).
# - delta_mean : maximum number of time units to shift the mean around mean_peaks (defaut=0)
# - delta_off_period_end : maximum number of time units to shift around off_period_end (defaut=0)
# - delta_off_period_start : maximum number of time units to shift around off_period_start (defaut=0)
# - delta_prob : maximum  number of units to shift the probability to be born in the first peak around the first value of probs. The second value of probs is always
# 1 (defaut=0)
# - delta_sd :  maximum number of time units to shift the standard deviation around sd_peaks (defaut=0)
# - variability : standard deviation of the normal distributions that define the new parameters (defaut=2). Warning: this parameter should not be set too high,
# otherwise the new parameters could be drawn outside of the cycle limits. The mean of the normal distribution is defined by a random drawing of a date in the
# interval [initial parametre - delta ; initial parametre + delta]. The function will draw values until it find one inside the limits of the cycle.

### VALUES
# The function returns a list of two elements: general_info and data_supracyc.
# data_supracyc: contains the phenology of births in the format selected (number of births for each time unit of each cycle).
# general_info: gathers all the settings chosen to simulate the phenology of births.

### REMARKS
# The function cannot deal with birth patterns spanning more than one cycle (if the number of births at the end of the cycle is not 0, it will be automatically
# reset for the following cycles). So the parameters should be selected to have a birth pattern that fits in one cycle.

### EXAMPLES
# birth_pattern <- simulation.pattern(format="row", graph=T, graph_format="ind",
#                                     births_distrib="peaks", nb_tu_per_cycle=365, nb_draws_per_tu=10, nb_cycles=6,
#                                     delta_prob=0, delta_mean=0, delta_sd=0, delta_off_period_start=0, delta_off_period_end=0, variability=0,
#                                     nb_peaks=2, mean_peaks=c(146, 219), sd_peaks=c(15, 20), probs=c(0.5, 1),
#                                     off_period_start=0, off_period_end=0)

# FUNCTIONS ########################################################################################################################################################

# generate a uniform distribution
unif.patt <- function(nb_tu_per_cycle=365, nb_draws_per_tu=10) {
  one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
  one_cycle_data$nb_births <- round(rep(nb_draws_per_tu, nrow(one_cycle_data))+rnorm(n=nrow(one_cycle_data), mean=0, sd=1), 0)
  return(one_cycle_data)
}

# generate a random distribution
rand.patt <- function(nb_tu_per_cycle=365, nb_draws_per_tu=10) {
  one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
  for (tu in 1:nb_tu_per_cycle) {
    one_cycle_data$nb_births[one_cycle_data$tu==tu] <- sample(0:nb_draws_per_tu, 1, replace=T)
  }
  return(one_cycle_data)
}

# generate a distribution with over dispersion of zeros and no/one/serveral empty period (following a negative binomial distribution)
negbin.patt <- function(nb_tu_per_cycle=365, nb_draws_per_tu=10, off_period_start=150, off_period_end=200) {
  # check initial conditions
  try(if(length(off_period_start)!=length(off_period_end)) stop(paste("off_period_end length differs from off_period_start length", sep=""), call.=F))
  for (i in 1:length(off_period_start)) {
    try(if(off_period_start[i] > off_period_end[i]) stop(paste("error: off_period_start nb", i, "greater than off_period_end nb", i, sep=" "), call.=F))
  }
  # initialisation of the parametres
  if (off_period_start[1]>0 & off_period_end[1]>0) { # if off_period_start and off_period_end > 0 -> an off period should be generated
    off_period_start2 <- off_period_start+round(rnorm(n=length(off_period_start), mean=0, sd=1), 0) # generate noise at the start/end of the period
    off_period_end2 <- off_period_end+round(rnorm(n=length(off_period_end), mean=0, sd=1), 0) 
    off_period <- NULL
    for (i in 1: length(off_period_start2)) {
      off_period <- c(off_period, seq(off_period_start2[i], off_period_end2[i], 1))
    }
  } else if (off_period_start[1]==0 & off_period_end[1]==0) { # if off_period_start and off_period_end = 0 -> no off period should be generated
    off_period <- NA
  }
  # draw "a female gives birth/a female does not give birth" event
  one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA)
  for (tu in 1:nb_tu_per_cycle) {
    one_cycle_data$nb_births[one_cycle_data$tu==tu] <- rnbinom(n=1, mu=1, size=1)
  }
  one_cycle_data$nb_births[one_cycle_data$tu %in% off_period] <- 0
  return(one_cycle_data)
}

# generate a multimodal distribution (based on a mixture of normal distributions)
norm.patt <- function(nb_peaks=4, nb_tu_per_cycle=365, nb_draws_per_tu=10, probs=c(0.2, 0.5, 0.7, 1), mean_peaks=c(73, 146, 219, 292), sd_peaks=rep(15, 4)) {
  # check initial conditions
  try(if(nb_peaks!=length(probs)) stop(paste("you asked for ", nb_peaks, " peak(s), probs length must be ", nb_peaks, sep=""), call.=F))
  try(if(nb_peaks!=length(mean_peaks)) stop(paste("you asked for ", nb_peaks, " peak(s), mean_peaks length must be ", nb_peaks, sep=""), call.=F))
  try(if(nb_peaks!=length(sd_peaks)) stop(paste("you asked for ", nb_peaks, " peak(s), sd_peaks length must be ", nb_peaks, sep=""), call.=F))
  # initialisation of the parametres
  den_param <- vector("list", length=nb_peaks)
  for (peak in 1:nb_peaks) {
    den_param[[peak]]$moy <- mean_peaks[peak]
    den_param[[peak]]$sd <- sd_peaks[peak]
  }
  # generate the probability densities
  den_rep <- vector("list", length=nb_peaks)
  for (peak in 1:nb_peaks) {
    den_rep[[peak]] <- dnorm(1:nb_tu_per_cycle, den_param[[peak]]$moy, den_param[[peak]]$sd) 
  }
  # draw "a female gives birth/a female does not give birth" event
  nb_draws_per_tu2 <- nb_draws_per_tu*100 # as density values are low, increase it to get enough birth events compared to the other distributions
  one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
  for(tu in 1:nb_tu_per_cycle) {
    dists <- runif(nb_draws_per_tu2, 0, 1)
    one_day_data <- vector(length=nb_draws_per_tu2)
    for (draw in 1:nb_draws_per_tu2) {
      for (peak in 1:(nb_peaks)) { # select in which distribution to draw the event
        if (dists[draw]<probs[peak]) {
          one_day_data[draw] <- rbinom(1, 1, den_rep[[peak]][tu])
          break;
        }
      }
    }
    one_cycle_data$nb_births[one_cycle_data$tu==tu] <- sum(one_day_data) # sum of the drawings for each female for each time unit
  }
  return(one_cycle_data)
}

# selection of the shape of the distribution
select.patt <- function(births_distrib="uniform", nb_peaks=4, nb_tu_per_cycle=365, nb_draws_per_tu=10, probs=c(0.2, 0.5, 0.7, 1), mean_peaks=c(73, 146, 219, 292),
                        sd_peaks=rep(15, 4), off_period_start=150, off_period_end=200) {
  # check initial conditions
  "%!in%" <- function(x, y) ! ("%in%"(x, y)) # exclusion function
  try(if(births_distrib %!in% c("uniform", "random", "negbinomial", "peaks")) stop("Error in births_distrib", call.=F))
  # select the shape of the distribution
  if (births_distrib=="uniform") {
    sortie <- unif.patt(nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu)
  } else if (births_distrib=="random") {
    sortie <- rand.patt(nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu)
  } else if (births_distrib=="peaks") {
    sortie <- norm.patt(nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, probs=probs, mean_peaks=mean_peaks, sd_peaks=sd_peaks)
  } else if (births_distrib=="negbinomial") {
    sortie <- negbin.patt(nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, off_period_start=off_period_start, off_period_end=off_period_end) 
  } else {
    sortie <- "ERROR"
  }
  return(sortie)
}

# simulate phenology of births
simulation.pattern <- function(format="col", graph=F, graph_format="ind", births_distrib="uniform", nb_tu_per_cycle=365, nb_draws_per_tu=10, nb_cycles=3,
                               delta_prob=0, delta_mean=0, delta_sd=0, delta_off_period_start=0, delta_off_period_end=0, variability=2, nb_peaks=4,
                               probs=c(0.2, 0.5, 0.7, 1), mean_peaks=c(73, 146, 219, 292), sd_peaks=rep(15, 4), off_period_start=150, off_period_end=200) {
  ### Initialisations
  library(zoo)
  library(tidyr)
  # check format
  "%!in%" <- function(x, y) ! ("%in%"(x, y)) # exclusion function
  try(if(format %!in% c("row", "col")) stop("Error in format", call.=F))
  # check consistency
  if (births_distrib %in% c("uniform", "random") | (delta_mean==0 & delta_sd==0 & delta_off_period_start==0 & delta_off_period_end==0 & delta_prob==0)) {
    consistency <- T
  } else {
    consistency <- F
  }
  # initialisation of the output table according to the format chosen
  if (format=="row") {
    data_supracyc <- NULL
  } else if (format=="col") {
    data_supracyc <-  data.frame(tu=seq(1, nb_tu_per_cycle, 1))  
  }
  ### Generate phenology of births for shapes "uniform", "random", "negbinomial" or "peaks" (same parameters for each cycle) 
  if (consistency==T) {
    for (cycle in 1:nb_cycles) {
      one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
      one_cycle_data <- select.patt(births_distrib=births_distrib, nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, probs=probs,
                                    mean_peaks=mean_peaks, sd_peaks=sd_peaks, off_period_start=off_period_start, off_period_end=off_period_end)
      if (format=="row") {
        one_cycle_data$cycle <- cycle
        data_supracyc <- rbind(data_supracyc, one_cycle_data)
      } else if (format=="col") {
        colnames(one_cycle_data)[colnames(one_cycle_data) %in% "nb_births"] <- paste("nb_births_cycle", cycle, sep="") # rename according to the cycle
        data_supracyc <- merge(data_supracyc, one_cycle_data, by="tu", all.x=T) # aggregate to the other cycles
        data_supracyc[ , paste("nb_births_cycle", cycle, sep="")][is.na(one_cycle_data[ , paste("nb_births_cycle", cycle, sep="")])==T] <- 0 # convert NA into 0
      }
    }
    ### Generate phenology of births for shapes "negbinomial" with one off period or "peaks" with one or two peaks (variable parameters for each cycle) 
  } else if (consistency==F & nb_cycles>1 & ((births_distrib=="peaks" & nb_peaks%in%c(1, 2)) | (births_distrib=="negbinomial" & length(off_period_start)==1))) {
    ## set initial conditions
    sd_peaks0 <- sd_peaks
    mean_peaks0 <- mean_peaks
    off_period_start0 <- off_period_start
    off_period_end0 <- off_period_end
    probs0 <- probs
    # draw births pattern for the first cycle
    one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
    one_cycle_data <- select.patt(births_distrib=births_distrib, nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, probs=probs,
                                  mean_peaks=mean_peaks, sd_peaks=sd_peaks, off_period_start=off_period_start, off_period_end=off_period_end) 
    if (format=="row") {
      one_cycle_data$cycle <- 1
      data_supracyc <- rbind(data_supracyc, one_cycle_data)
    } else if (format=="col") {
      colnames(one_cycle_data)[colnames(one_cycle_data) %in% "nb_births"] <- "nb_births_cycle1" # rename according to the cycle
      data_supracyc <- merge(data_supracyc, one_cycle_data, by="tu", all.x=T) # aggregate to the other cycles
      data_supracyc$nb_births_cycle1[is.na(data_supracyc$nb_births_cycle1)==T] <- 0 # convert NA into 0
    }
    # draw birth pattern for the following cycles
    for (cycle in 2:nb_cycles) {
      resol_mean <- delta_mean/10
      if (delta_mean>0) { # draw new mean
        new_ref_mean <- seq(mean_peaks0-delta_mean, mean_peaks0+delta_mean, resol_mean) # references for the new parametres
        try(if(min(new_ref_mean)<0 | max(new_ref_mean)>nb_tu_per_cycle) stop("mean(s) births date out of range", call.=F)) # check conditions
        mean_peaks <- 0 # generate new parametres
        while (mean_peaks<=0 | mean_peaks>nb_tu_per_cycle) { # draw mean_peaks until it falls into the adequate interval
          mean_peaks <- round(rnorm(1, mean=sample(new_ref_mean, 1), sd=variability), 0)
        }
      }
      if (delta_sd>0) { # draw new sd
        resol_sd <- delta_sd/5
        new_ref_sd <- seq(sd_peaks0-delta_sd, sd_peaks0+delta_sd, resol_sd)
        try(if(min(new_ref_sd)<0 | max(new_ref_sd)>nb_tu_per_cycle) stop("standard deviation(s) of births date out of range", call.=F))
        sd_peaks <- 0
        while (sd_peaks<=0 | sd_peaks>nb_tu_per_cycle) { 
          sd_peaks <- round(rnorm(1, mean=sample(new_ref_sd, 1), sd=variability), 0)
        }
      }
      if (delta_off_period_start>0) { # draw new start of the off period
        resol_off_period_start <- delta_off_period_start/10
        new_ref_off_period_start <- seq(off_period_start0-delta_off_period_start, off_period_start0+delta_off_period_start, resol_off_period_start)
        try(if(min(new_ref_off_period_start)<0 | max(new_ref_off_period_start)>nb_tu_per_cycle) stop("off period start out of range", call.=F))
        off_period_start <- 0
        while (off_period_start<=0 | off_period_start>nb_tu_per_cycle) {
          off_period_start <- round(rnorm(1, mean=sample(new_ref_off_period_start, 1), sd=variability), 0)
        }
      }
      if (delta_off_period_end>0) { # draw new end of the off period
        resol_off_period_end <- delta_off_period_end/10
        new_ref_off_period_end <- seq(off_period_end0-delta_off_period_end, off_period_end0+delta_off_period_end, resol_off_period_end)
        try(if(min(new_ref_off_period_end)<0 | max(new_ref_off_period_end)>nb_tu_per_cycle) stop("off period end out of range", call.=F))
        off_period_end <- 0
        while (off_period_end<=0 | off_period_end>nb_tu_per_cycle) {
          off_period_end <- round(rnorm(1, mean=sample(new_ref_off_period_end, 1), sd=variability), 0)
        }
      }
      if (delta_prob>0) { # draw the probability to be born during the first or the second peak
        resol_prob <- delta_prob/5
        new_ref_prob <- seq(probs0[1]-delta_prob, probs0[1]+delta_prob, resol_prob)
        try(if(min(new_ref_prob)<0 | max(new_ref_prob)>1) stop("probability of birth in either peak out of range", call.=F))
        probs[1] <- 0
        while (probs[1]<=0 | probs[1]>1) {
          probs[1] <- round(rnorm(1, mean=sample(new_ref_prob, 1), sd=variability/100), 3)
        }
      }
      # draw births pattern with the new parameters
      one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
      one_cycle_data <- select.patt(births_distrib=births_distrib, nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, probs=probs,
                                    mean_peaks=mean_peaks, sd_peaks=sd_peaks, off_period_start=off_period_start, off_period_end=off_period_end) 
      if (format=="row") {
        if (nb_cycles>1) {
          one_cycle_data$cycle <- cycle
        }
        data_supracyc <- rbind(data_supracyc, one_cycle_data)
      } else if (format=="col") {
        colnames(one_cycle_data)[colnames(one_cycle_data) %in% "nb_births"] <- paste("nb_births_cycle", cycle, sep="") # rename according to the cycle
        data_supracyc <- merge(data_supracyc, one_cycle_data, by="tu", all.x=T) # aggregate to the other cycles
        data_supracyc[ , paste("nb_births_cycle", cycle, sep="")][is.na(one_cycle_data[ , paste("nb_births_cycle", cycle, sep="")])==T] <- 0 # convert NA into 0
      }
    }
  }
  ### Final formatting of the outputs
  if (format=="row") {
    data_supracyc <- data_supracyc[ , c("cycle", "tu", "nb_births")]
  }
  ### Graphical representation
  # fonction detecting the highest value of a data frame
  find.max.df <- function(df, col_ref) {
    check <- NULL
    tab <- as.data.frame(df[ , -which(names(df) %in% as.character(col_ref))])
    for (i in 1:ncol(tab)) {check <- c(check, tab[, i])}
    return(max(check))
  }
  # graph representation
  if (graph==T & format=="col") {
    if (graph_format=="sup") {
      plot(data_supracyc$nb_births_cycle1~data_supracyc$tu, type="l", xlab="Cycle Time Unit", ylab="Number of births", xlim=c(1, nb_tu_per_cycle),
           ylim=c(0, find.max.df(data_supracyc, "tu")+1), main=paste("Number of births through", nb_cycles, "cycles", sep=" "))
      if (nb_cycles>1) {
        for (cycle in 2:nb_cycles) {
          lines(data_supracyc[ , paste("nb_births_cycle", cycle, sep="")], col=cycle)
        }
      }
    } else if (graph_format=="ind") {
      data_plot <- read.zoo(data_supracyc)
      plot(na.approx(data_plot, na.rm=F), main="Number of births through each cycle", xlab="Cycle Time Unit", ylab=paste("Cycle ", seq(1, nb_cycles, 1), sep=""),
           ylim=c(0, find.max.df(data_supracyc, "tu")+1))
    }
  } else if (graph==T & format=="row") {
    if (graph_format=="sup") {
      plot(data_supracyc$nb_births[data_supracyc$cycle==1]~data_supracyc$tu[data_supracyc$cycle==1], type="l", xlab="Cycle Time Unit", ylab="Number of births",
           xlim=c(1, nb_tu_per_cycle), ylim=c(0, find.max.df(data_supracyc, "tu")+1), main=paste("Number of births through", nb_cycles, "cycles", sep=" "))
      if (nb_cycles>1) {
        for (cycle in 2:nb_cycles) {
          lines(data_supracyc$nb_births[data_supracyc$cycle==cycle], col=cycle)
        }
      }
    } 
    else if (graph_format=="ind") { 
      data_supracyc2 <- spread(data_supracyc, key="cycle", value="nb_births")
      data_plot <- read.zoo(data_supracyc2)
      plot(na.approx(data_plot, na.rm=F), main="Number of births through each cycle", xlab="Cycle Time Unit", ylab=paste("Cycle ", seq(1, nb_cycles, 1), sep=""),
           ylim=c(0, find.max.df(data_supracyc2, "tu")+1))
    }
  }
  ### Output
  return(list(general_info=list(format=format, births_distrib=births_distrib, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, nb_cycles=nb_cycles,
                                delta_mean=delta_mean, delta_sd=delta_sd, nb_peaks=nb_peaks, probs=probs, mean_peaks=mean_peaks, sd_peaks=sd_peaks,
                                off_period_start=off_period_start, off_period_end=off_period_end), data_supracyc=data_supracyc))
}

####################################################################################################################################################################
####################################################################################################################################################################
