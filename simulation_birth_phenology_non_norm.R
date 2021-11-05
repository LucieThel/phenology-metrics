####################################################################################################################################################################
#### SIMULATION : PHENOLOGY OF BIRTHS FROM NON-NORMAL DISTRIBUTIONS ################################################################################################
####################################################################################################################################################################

# DOCUMENTATION ####################################################################################################################################################

### DESCRIPTION
# This group of functions generates patterns of births following several kinds of non-normal distributions.

### ARGUMENTS
## General parameters:
# - births_distrib : type of the distribution; "peaks" = births are drawn from a normal law; "random" = births are randomly distributed; "skewnorm" = births are
# drawn from a skewd normal distribution; "cauchy" = births are drawn from a Cauchy distribution (default="peaks")
# - format : format of the output when more than one reproductive cycle is generated, "row" = by lines, "col" = by columns (default="col")
# - graph : plot or not the graph associated to the simulations (default=FALSE)
# - graph_format : "ind" = represents the output cycle by cycle on separate graphs, "sup" = superimposes all the reproductive cycles on the same graph (default=
# "ind")
# - nb_cycles : number of cycles (e.g. years) during which simulating patterns of births (default=3)
# - nb_draws_per_tu : number of females that could give birth every time unit (default=10)
# - nb_tu_per_cycle : number of time unit (e.g. days) among which distributing births (default=365)

## specific parameters:
# - delta_mean1_mean2 : duration between first and second mean date of births when choosing a bimodal distribution defined according to the first peak (in the case
# of a normal distribution) (default=0)
# - mean_peaks : mean(s) of the normal distribution(s) (default=c(73, 146, 219, 292))
# - nb_peaks : number of peaks during one cycle (default=4)
# - probs : probabilities indicating in which probability density to draw the event "a female gives birth/a female does not give birth". The last value is
# necessarily 1. It corresponds to the cumulative proportion of births occurring during each peak (default=c(0.2, 0.5, 0.7, 1))
# - sd_peaks : standard deviation of the normal distribution(s) (default=c(15, 15, 15, 15))
# - skewn : degree of skewness in the case of a skewed distribution (default=0)

## Consistency across reproductive cycles:
# It is possible to modify the level of consistency between the reproductive cycles for some parameters of the distributions:
# - delta_mean : maximum number of time units to shift the mean around mean_peaks (default=0)
# - delta_prob : maximum  number of units to shift the probability to be born in the first peak around the first value of probs. The second value of probs is always
# 1 (default=0)
# - delta_sd :  maximum number of time units to shift the standard deviation around sd_peaks (default=0)
# - variability : standard deviation of the normal distributions that define the new parameters (default=2). Warning: this parameter should not be set too high,
# otherwise the new parameters could be drawn outside of the cycle limits. The mean of the normal distribution is defined by a random drawing of a date in the
# interval [initial parameter - delta ; initial parameter + delta]. The function will draw values until it find one inside the limits of the cycle.

### VALUES
# The function returns a list of two elements: general_info and data_supracyc.
# data_supracyc: contains the phenology of births in the format selected (number of births for each time unit of each reproductive cycle).
# general_info: gathers all the settings chosen to simulate the phenology of births.

### REMARKS
# The function cannot deal with birth patterns spanning more than one cycle (if the number of births at the end of the cycle is not 0, it will be automatically
# reset for the following cycles). So the parameters should be selected to create a pattern of births that fits in one cycle duration.
# The function is not designed to allow variations of the skewness for a skewed normal distribution across several reproductive cycles.

### EXAMPLES
# random distribution of births:
# birth_patt1 <- simulation.pattern(format="row", graph=T, graph_format="ind", births_distrib="random", nb_tu_per_cycle=365, nb_draws_per_tu=10, nb_cycles=6,
#                                   delta_mean=0, delta_sd=0, variability=0, nb_peaks=0, probs=0, mean_peaks=0, sd_peaks=0, skewn=0, delta_mean1_mean2=0)
# skewed normal distribution of births:
# birth_patt2 <- simulation.pattern(format="row", graph=T, graph_format="ind", births_distrib="skewnorm", nb_tu_per_cycle=365, nb_draws_per_tu=10, nb_cycles=6,
#                                   delta_mean=0, delta_sd=0, variability=0, nb_peaks=1, probs=1, mean_peaks=183, sd_peaks=30, skewn=2, delta_mean1_mean2=0)
# bimodal distribution of births:
# birth_patt3 <- simulation.pattern(format="row", graph=T, graph_format="ind", births_distrib="peaks", nb_tu_per_cycle=365, nb_draws_per_tu=10, nb_cycles=6,
#                                   delta_mean=0, delta_sd=0, variability=0, nb_peaks=2, probs=c(0.8, 1), mean_peaks=122, sd_peaks=c(30, 15), skewn=0,
#                                   delta_mean1_mean2=70)
# Cauchy distribution
# birth_patt4 <- simulation.pattern(format="row", graph=T, graph_format="ind", births_distrib="cauchy", nb_tu_per_cycle=365, nb_draws_per_tu=10, nb_cycles=6,
#                                   delta_mean=0, delta_sd=0, variability=0, nb_peaks=1, probs=1, mean_peaks=183, sd_peaks=30, skewn=0, delta_mean1_mean2=0)

# FUNCTIONS ########################################################################################################################################################

# # generate a random distribution
rand.patt <- function(nb_tu_per_cycle=365, nb_draws_per_tu=10) {
  one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
  for (tu in 1:nb_tu_per_cycle) {
    one_cycle_data$nb_births[one_cycle_data$tu==tu] <- sample(0:nb_draws_per_tu, 1, replace=T)
  }
  return(one_cycle_data)
}

# generate a multimodal distribution (based on a mixture of normal distributions)
norm.patt <- function(nb_peaks=4, nb_tu_per_cycle=365, nb_draws_per_tu=10, probs=c(0.2, 0.5, 0.7, 1), mean_peaks=c(73, 146, 219, 292), sd_peaks=rep(15, 4),
                      delta_mean1_mean2=0) {
  # check initial conditions 
  try(if(nb_peaks!=length(probs)) stop(paste("you asked for ", nb_peaks, " peak(s), probs length must be ", nb_peaks, sep=""), call.=F))
  try(if(nb_peaks!=length(mean_peaks) & delta_mean1_mean2>0) stop(paste("you asked for ", nb_peaks, " peak(s), mean_peaks length must be ", nb_peaks, sep=""),
                                                                  call.=F)) # si mean2 pas definie par mean1
  try(if(nb_peaks!=length(sd_peaks)) stop(paste("you asked for ", nb_peaks, " peak(s), sd_peaks length must be ", nb_peaks, sep=""), call.=F))
  # initialisation of the parameters
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

# generate a skewed normal distribution
skewnorm.patt <- function(nb_tu_per_cycle=365, nb_draws_per_tu=10, mean_peaks=183, sd_peaks=15, skewn=2) {
  library(fGarch)
  # generate the probability densities
  den <- dsnorm(x=1:nb_tu_per_cycle, mean=mean_peaks, sd=sd_peaks, xi=skewn, log=F)
  # draw "a female gives birth/a female does not give birth" event
  nb_draws_per_tu2 <- nb_draws_per_tu*100 # as density values are low, increase it to get enough birth events compared to the other distributions
  one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
  for(tu in 1:nb_tu_per_cycle) {
    one_day_data <- vector(length=nb_draws_per_tu2)
    for (draw in 1:nb_draws_per_tu2) {
      one_day_data[draw] <- rbinom(1, 1, den[tu])
    }
    one_cycle_data$nb_births[one_cycle_data$tu==tu] <- sum(one_day_data) # sum of the drawings for each female for each time unit
  }
  return(one_cycle_data)
}

# generate a Cauchy distribution
cauchy.patt <- function(nb_tu_per_cycle=365, nb_draws_per_tu=10, mean_peaks=183, sd_peaks=15) {
  # generate the probability densities
  den <- dcauchy(x=1:nb_tu_per_cycle, location=mean_peaks, scale=sd_peaks)
  # draw "a female gives birth/a female does not give birth" event
  nb_draws_per_tu2 <- nb_draws_per_tu*100 # as density values are low, increase it to get enough birth events compared to the other distributions
  one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
  for(tu in 1:nb_tu_per_cycle) {
    one_day_data <- vector(length=nb_draws_per_tu2)
    for (draw in 1:nb_draws_per_tu2) {
      one_day_data[draw] <- rbinom(1, 1, den[tu])
    }
    one_cycle_data$nb_births[one_cycle_data$tu==tu] <- sum(one_day_data) # sum of the drawings for each female for each time unit
  }
  return(one_cycle_data)
}

# selection of the shape of the distribution (here, the only choice is "peaks")
select.patt <- function(births_distrib="peaks", nb_peaks=4, nb_tu_per_cycle=365, nb_draws_per_tu=10, probs=c(0.2, 0.5, 0.7, 1), mean_peaks=c(73, 146, 219, 292),
                        sd_peaks=rep(15, 4), skewn=0) {
  # check initial conditions
  "%!in%" <- function(x, y) ! ("%in%"(x, y)) # exclusion function
  try(if(births_distrib %!in% c("random", "peaks", "skewnorm", "cauchy")) stop("Error in births_distrib", call.=F))
  # select the shape of the distribution
  if (births_distrib=="random") {
    sortie <- rand.patt(nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu)
  } else if (births_distrib=="peaks") {
    sortie <- norm.patt(nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, probs=probs, mean_peaks=mean_peaks, sd_peaks=sd_peaks)
  } else if (births_distrib=="skewnorm") {
    sortie <- skewnorm.patt(nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, mean_peaks=mean_peaks, sd_peaks=sd_peaks, skewn=skewn)
  } else if (births_distrib=="cauchy") {
    sortie <- cauchy.patt(nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, mean_peaks=mean_peaks, sd_peaks=sd_peaks)
  } else {
    sortie <- "ERROR"
  }
  return(sortie) 
}

# simulate phenology of births
simulation.pattern <- function(format="col", graph=F, graph_format="ind", births_distrib="peaks", nb_tu_per_cycle=365, nb_draws_per_tu=10, nb_cycles=3,
                               delta_mean=0, delta_sd=0, variability=2, nb_peaks=4, probs=c(0.2, 0.5, 0.7, 1), mean_peaks=c(73, 146, 219, 292),
                               sd_peaks=rep(15, 4), skewn=0, delta_mean1_mean2=0) {
  ### Initialisations
  library(zoo)
  library(tidyr)
  # check format
  "%!in%" <- function(x, y) ! ("%in%"(x, y)) # exclusion function
  try(if(format %!in% c("row", "col")) stop("Error in format", call.=F))
  # check consistency
  if (births_distrib %in% c("random") | (delta_mean==0 & delta_sd==0)) {
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
  ### Generate phenology of births with the same parameters for each reproductive cycle
  if (consistency==T) {
    for (cycle in 1:nb_cycles) {
      if (delta_mean1_mean2>0) { mean_peaks <- c(mean_peaks, mean_peaks+delta_mean1_mean2) }
      one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA)
      one_cycle_data <- select.patt(births_distrib=births_distrib, nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, probs=probs,
                                    mean_peaks=mean_peaks, sd_peaks=sd_peaks, skewn=skewn)
      if (format=="row") {
        one_cycle_data$cycle <- cycle
        data_supracyc <- rbind(data_supracyc, one_cycle_data)
      } else if (format=="col") {
        colnames(one_cycle_data)[colnames(one_cycle_data) %in% "nb_births"] <- paste("nb_births_cycle", cycle, sep="") # rename according to the cycle
        data_supracyc <- merge(data_supracyc, one_cycle_data, by="tu", all.x=T) # aggregate to the other cycles
        data_supracyc[ , paste("nb_births_cycle", cycle, sep="")][is.na(one_cycle_data[ , paste("nb_births_cycle", cycle, sep="")])==T] <- 0 # convert NA into 0
      }
    }
    ### Generate phenology of births with one or two peaks and variable parameters for each cycle
  } else if (consistency==F & nb_cycles>1 & ((births_distrib=="peaks" & nb_peaks%in%c(1, 2)) | (births_distrib=="skewnorm") | (births_distrib=="cauchy"))) {
    if (nb_peaks==2 & (delta_mean>0 | delta_sd>0)) { # in the case of a bimodal distribution with variations on the first peak only
      sd_peaks0 <- sd_peaks[1]
      sd_peaks2 <- sd_peaks[2]
      mean_peaks0 <- mean_peaks
      # draw births pattern for the first cycle
      one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
      one_cycle_data <- select.patt(births_distrib=births_distrib, nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, probs=probs,
                                    mean_peaks=c(mean_peaks, mean_peaks+delta_mean1_mean2), sd_peaks=sd_peaks, skewn=skewn)
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
          new_ref_mean <- seq(mean_peaks0-delta_mean, mean_peaks0+delta_mean, resol_mean) # references for the new parameters
          try(if(min(new_ref_mean)<0 | max(new_ref_mean)>nb_tu_per_cycle) stop("mean(s) births date out of range", call.=F)) # check conditions
          mean_peaks <- 0 # generate new parameters
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
          sd_peaks <- c(sd_peaks, sd_peaks2)
        }
        # draw births pattern with the new parameters
        one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
        one_cycle_data <- select.patt(births_distrib=births_distrib, nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu,
                                      probs=probs, mean_peaks=c(mean_peaks, mean_peaks+delta_mean1_mean2), sd_peaks=sd_peaks, skewn=skewn)
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
    } else { # in the other cases
      ## set initial conditions
      sd_peaks0 <- sd_peaks
      mean_peaks0 <- mean_peaks
      probs0 <- probs
      # draw births pattern for the first cycle
      one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
      one_cycle_data <- select.patt(births_distrib=births_distrib, nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, probs=probs,
                                    mean_peaks=mean_peaks, sd_peaks=sd_peaks, skewn=skewn)
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
          new_ref_mean <- seq(mean_peaks0-delta_mean, mean_peaks0+delta_mean, resol_mean) # references for the new parameters
          try(if(min(new_ref_mean)<0 | max(new_ref_mean)>nb_tu_per_cycle) stop("mean(s) births date out of range", call.=F)) # check conditions
          mean_peaks <- 0 # generate new parameters
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
        # draw births pattern with the new parameters
        one_cycle_data <- data.frame(tu=seq(1, nb_tu_per_cycle, 1), nb_births=NA) 
        one_cycle_data <- select.patt(births_distrib=births_distrib, nb_peaks=nb_peaks, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu,
                                      probs=probs, mean_peaks=mean_peaks, sd_peaks=sd_peaks, skewn=skewn)
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
    } else if (graph_format=="ind") { 
      data_supracyc2 <- spread(data_supracyc, key="cycle", value="nb_births")
      data_plot <- read.zoo(data_supracyc2)
      plot(na.approx(data_plot, na.rm=F), main="Number of births through each cycle", xlab="Cycle Time Unit", ylab=paste("Cycle ", seq(1, nb_cycles, 1), sep=""),
           ylim=c(0, find.max.df(data_supracyc2, "tu")+1))
    }
  }
  ### Output
  return(list(general_info=list(format=format, births_distrib=births_distrib, nb_tu_per_cycle=nb_tu_per_cycle, nb_draws_per_tu=nb_draws_per_tu, nb_cycles=nb_cycles,
                                delta_mean=delta_mean, delta_sd=delta_sd, delta_mean1_mean2=delta_mean1_mean2, nb_peaks=nb_peaks, probs=probs,
                                mean_peaks=mean_peaks, sd_peaks=sd_peaks, skewn=skewn), 
              data_supracyc=data_supracyc))
}

####################################################################################################################################################################
####################################################################################################################################################################
