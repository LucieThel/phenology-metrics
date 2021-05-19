####################################################################################################################################################################
#### TIMING - SYNCHRONY - RHYTHMICITY - REGULARITY METRICS EXTRACTED FROM THE LITERATURE ASSOCIATED TO BIRTH PHENOLOGY IN LARGE HERBIVORES #########################
####################################################################################################################################################################

# DOCUMENTATION ####################################################################################################################################################

### DESCRIPTION
# Code of the different metrics used to evaluate timing, synchrony, rhythmicity and regularity in the phenology of births in large herbivores.

### ARGUMENTS - VALUES - DETAILS
# For each function, a "graph" argument allows to choose if a plot should be returned or not.
# The functions pheno.calabrese(), pheno.caughley(), pheno.bowyer() and pheno.riedmanmeng() led to code adjustments and could properly work only in the context of
# our analysis.
# The function pheno.zerbe() calls Perl, so necessitates Perl (run with perl 5, version 26, subversion 1 (v5.26.1)) to be installed on the computer running the
# function. An intermediary file is created, so a working directory should be specified too.
# See details in each function.

### WORKING DIRECTORY
working_dir <- "/home/lthel/Bureau/Rscripts metrics"

# FUNCTIONS ########################################################################################################################################################

# MEAN, MEAN WEIGHTED, MEDIAN, RANDOMNESS INDEX, SKEWNESS 1, SKEWNESS 2, STANDARD DEVIATION
pheno.findlaylambin <- function(pattern, graph=F) { 
  # DOCUMENTATION
  ## Description: basic description of the birth pattern (mean, weigthed mean, median, standard deviation, two different skewness, randomness index).
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: mean, mean weighted, median, randomness index, standard deviation, two different skewness.
  ## Details: none.
  
  # INITIALISATIONS
  library(e1071)
  colnames(pattern) <- c("unit", "nb_births")
  data_int <- NULL
  
  # FUNCTION
  mean_birth_date <- sum(pattern$unit*pattern$nb_births, na.rm=T)/sum(pattern$nb_births, na.rm=T)
  var_weighted <- sum(pattern$nb_births*((pattern$unit-mean_birth_date)^2), na.rm=T)/sum(pattern$nb_births, na.rm=T)
  r_index <- var_weighted/mean_birth_date
  stand_dev <- sqrt(var_weighted)
  pattern$cum <- cumsum(pattern$nb_births)
  med_birth <- which(pattern$cum>=sum(pattern$nb_births)/2)[1]
  mode_birth_date <- which(pattern$nb_births==max(pattern$nb_births)) # may be more than one. Take the first if more than one
  skew1 <- (mean_birth_date-mode_birth_date[1])/stand_dev
  skew2 <- skewness(rep(pattern$unit, pattern$nb_births))
  mean_weighted <- mean_birth_date*(1/var_weighted)
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    abline(v=c(mean_birth_date, med_birth), col=c("brown", "blue"))
    text(mean_birth_date, max(pattern$nb_births)-1, "Mean birth date", col="brown")
    text(med_birth, max(pattern$nb_births)-1, "Median birth date", col="blue")
  }
  
  # OUTPUT
  return(list(mean=mean_birth_date, mean_weighted=mean_weighted, standard_deviation=stand_dev, median=med_birth, skewness_mean=skew1,
              skewness_var=skew2, randomness_index=r_index))
}

# VARIANCE, VARIANCE CORR
pheno.johnson <- function(pattern, interval="check", corr=T, graph=F) { 
  # DOCUMENTATION
  ## Description: weighted variance corrected by the Sheppard method (generalized to irregularly sampled data or not) or not.
  ## Arguments: 
  # check = specify if the function should select the adapted correction automatically.
  # corr = specify if the correction should be applied or not.
  # interval = specify which correction to apply: "regular" for a regular sampling and a Sheppard correction, "irregular" for a irregular sampling and a Sheppard
  # correction generalised to irregular data.
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: weighted variance or weighted variance corrected by Sheppard method (generalized to irregularly sampled data or not).
  ## Details: no graphical representation available.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  n <- sum(pattern$nb_births)
  h <- pattern$unit[2]-pattern$unit[1] 
  pattern$hi[1] <- 0
  for (i in 2:nrow(pattern)) {
    pattern$hi[i] <- pattern$unit[i]-pattern$unit[i-1]
  }
  if (interval=="check") {
    if (length(unique(pattern$hi))>2) { interval <- "irregular" } else { interval <- "regular" }
  }
  pattern <- rbind(c(0, 0), pattern)
  
  # FUNCTION
  if (corr==T) {
    if (interval=="irregular") {
      pattern$unitbis[1] <- 0
      pattern$int_mean[1] <- 0
      for (i in 2:nrow(pattern)) {
        pattern$unitbis[i] <- pattern$unit[i-1]+((pattern$unit[i]-pattern$unit[i-1])/2) 
      }
      mean_birth_datebis <- sum(pattern$unitbis*pattern$nb_births, na.rm=T)/n
      var_weighted <- sum(pattern$nb_births*((pattern$unitbis-mean_birth_datebis)^2), na.rm=T)/n
      pattern$nihi <- pattern$nb_births*pattern$hi
      variance <- var_weighted-sum(pattern$nihi, na.rm=T)/(12*n)
    } else if (interval=="regular") {
      mean_birth_date <- sum(pattern$unit*pattern$nb_births, na.rm=T)/n
      var_weighted <- sum(pattern$nb_births*((pattern$unit-mean_birth_date)^2), na.rm=T)/n
      variance <- var_weighted-h/12
    }
  } else if (corr==F) {
    mean_birth_date <- sum(pattern$unit*pattern$nb_births, na.rm=T)/sum(pattern$nb_births, na.rm=T)
    variance <- sum(pattern$nb_births*((pattern$unit-mean_birth_date)^2), na.rm=T)/sum(pattern$nb_births, na.rm=T)
  }
  
  # PLOT
  if (graph==T) {
    # no representation
  }
  
  # OUTPUT
  return(variance) # according to the correction, return variance or variance_corr
}

# BEGINNING PERIOD, FREQUENCY, PERIOD
pheno.bunnellyomtov <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: first and last date of birth, birth period length, birth frequency.
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: birth frequency, birth period length, first birth date, last birth date.
  ## Details: birth frequency is calculated as 1/(end-start).
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  
  # FUNCTION
  first_birth <- which(pattern$nb_births>0)[1] 
  last_birth <- which(pattern$nb_births>0)[length(which(pattern$nb_births>0))] 
  period <- last_birth-first_birth+1
  freq <- 1/period
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    abline(v=c(first_birth, last_birth), col="purple")
    text(first_birth, max(pattern$nb_births)-1, "Start birth period", col="purple")
    text(last_birth, max(pattern$nb_births)-1, "End birth period", col="purple")
  }
  
  # OUTPUT
  return(list(period=period, beginning_period=first_birth, end_births_period=last_birth, frequency=freq))
}

# NB TU MINIMAL BIRTHS
pheno.moe1 <- function(pattern, percent_min=1, consecutive=F, graph=F) {
  # DOCUMENTATION
  ## Description: number and position of time units with at least a certain percent of births.
  ## Arguments: 
  # consecutive = specify if time units should be counted only if all the following time units gather at least a certain percentage or even if they are mixed with
  # time units under the certain percentage of births.
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # percent_min = specify the percentage threshold.
  ## Values: number of time units with at least a certain percent of births, position of the time units with at least a certain percent of births.
  ## Details: none.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  perc <- percent_min/100*sum(pattern$nb_births)
  count <- 0
  all_count <- NULL
  when <- NULL
  
  # FUNCTION  
  if (consecutive==F) {
    for (i in 1:nrow(pattern)) {
      if (pattern$nb_births[i]>=perc) {
        count <- count+1
        when <- c(when, i) 
      }
      all_count <- count
    }
  } else if (consecutive==T) {
    for (i in 1:nrow(pattern)) {
      if (pattern$nb_births[i]>=perc) {
        count <- count+1
        when <- c(when, i) 
      } else if (pattern$nb_births[i]<perc) {
        all_count <- c(all_count, count)
        count <- 0
      }
    }
    all_count <- c(all_count, count) # to count the last time unit if it is not null, but the counter was reintialised just before
    all_count <- all_count[all_count>0] 
  }
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    points(rep(0, length(when))~when, col="red", pch=16)
  }
  
  # OUTPUT
  return(list(nb_tu_minimal_births=all_count, which_temporal_units=when))
}

# PIELOU
pheno.sinclair <- function(pattern, graph=F) { 
  # DOCUMENTATION
  ## Description: evenness index of the distribution of births.
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: evenness index.
  ## Details: if the evenness index is closed to 0, the births are highly synchronous. If it is closed to 1, they are not (they can be randomly, uniformly, ...
  # distributed). No graphical representation available.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  
  # FUNCTION
  pattern$prop_births <- pattern$nb_births/sum(pattern$nb_births, na.rm=T)
  pattern$interm <- pattern$prop_births*log(pattern$prop_births)
  pielou <- abs(sum(pattern$interm, na.rm=T)/log(nrow(pattern)))
  
  # PLOT
  if (graph==T) {
    # no representation
  }
  
  # OUTPUT
  return(pielou)
}

# BEGINNING PERIOD THRESHOLD, INTERQUANTILES PERIOD Q, PROP BIRTHS INTERQUANTILE PERIOD Q
pheno.gaillardmajluf <- function(pattern, first_quantile=NA, last_quantile=NA, percent=NA, percent_min=NA, percent_min_end=NA, graph=F) { 
  # DOCUMENTATION
  ## Description: duration of the birth period between the two selected quantiles, number and proportion of births taking place during this period, start and end of
  # the birth period when a certain percentage of births have been exceeded, duration of the period gathering a certain percentage of births based on the detection
  # of the associated quantiles.
  ## Arguments: 
  # first_quantile = specifiy the lower limit of the birth period to spot. Give the value in percent.
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # last_quantile = specifiy the upper limit of the birth period to spot. Give the value in percent.
  # percent = specifiy the percentage threshold to estimate quantiles determining the birth period.
  # percent_min = specifiy the minimum percentage of births to exceed to detect the start of the birth period.
  # percent_min_end = specifiy the minimum percentage of births to exceed to detect the end of the birth period.
  ## Values: duration of the birth period between the two selected quantiles, duration of the period gathering a certain percentage of births, end of the birth
  # period, number of births taking place inside the birth period, proportion of births taking place inside the birth period, start of the birth period.
  ## Details: none.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  pattern$cum <- cumsum(pattern$nb_births)
  Qa_birth <- NA
  Qb_birth <- NA
  nb_births <- NA
  prop_birth <- NA
  diff <- NA
  Qmin_birth <- NA
  Qmax_birth <- NA
  
  # FUNCTION
  if (is.na(first_quantile)==F | is.na(last_quantile)==F | is.na(percent)==F) {
    if (is.na(first_quantile)==F & is.na(last_quantile)==F & is.na(percent)==T) { # if ask for quantiles
      Qa_birth <- which(pattern$cum>=(first_quantile/100)*sum(pattern$nb_births, na.rm=T))[1]
      Qb_birth <- which(pattern$cum>=(last_quantile/100)*sum(pattern$nb_births, na.rm=T))[1]
    } else if (is.na(percent)==F & is.na(first_quantile)==T & is.na(last_quantile)==T) { # if ask for a proportion of births
      lim_inf <- (100-percent)/2 
      lim_inf_data <- lim_inf/100*sum(pattern$nb_births)
      lim_sup <- 100-((100-percent)/2) 
      lim_sup_data <- lim_sup/100*sum(pattern$nb_births)
      i <- 1
      while (pattern$cum[i]<lim_inf_data) { i <- i+1 }
      Qa_birth <- pattern$unit[i]
      i <- 1
      while (pattern$cum[i]<lim_sup_data) { i <- i+1 }
      Qb_birth <- pattern$unit[i]
    } else { # if error in the arguments
      try(if(is.na(percent)==F & is.na(first_quantile)==F & is.na(last_quantile)==F) stop("Try with only quantiles or only percent", call.=F)) 
    }
    diff <- Qb_birth-Qa_birth+1
    nb_births <- sum(pattern$nb_births[Qa_birth:Qb_birth], na.rm=T)
    prop_birth <- (nb_births/sum(pattern$nb_births, na.rm=T))*100
  }
  # start of the birth period
  if (is.na(percent_min)==F) {
    Qmin_birth <- which(pattern$cum>=(percent_min/100)*sum(pattern$nb_births, na.rm=T))[1]
  }
  # end of the birth period
  if (is.na(percent_min_end)==F) {
    Qmax_birth <- which(pattern$cum>=(sum(pattern$nb_births, na.rm=T)-(percent_min_end/100)*sum(pattern$nb_births, na.rm=T)))[1]
  }
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    if (is.na(first_quantile)==F | is.na(last_quantile)==F | is.na(percent)==F) {
      abline(v=c(Qa_birth, Qb_birth), col="purple")
      text(Qa_birth, max(pattern$nb_births)-1, "First quantile", col="purple")
      text(Qb_birth, max(pattern$nb_births)-1, "Last quantile", col="purple")
    } else if (is.na(percent_min)==F) {
      abline(v=c(Qmin_birth), col="purple")
      text(Qmin_birth, max(pattern$nb_births)-1, "First births", col="purple")
    } else if (is.na(percent_min_end)==F) {
      abline(v=c(Qmax_birth), col="purple")
      text(Qmax_birth, max(pattern$nb_births)-1, "Last births", col="purple")
    }
  }
  
  # OUTPUT
  return(list(interquantiles_period=diff, first_quantile=Qa_birth, last_quantile=Qb_birth, nb_births=nb_births, 
              prop_births_interquantiles_period=prop_birth, beginning_period_threshold=Qmin_birth, end_period_threshold=Qmax_birth))
}

# CENTRE
pheno.sigouin <- function(pattern, graph=F) { 
  # DOCUMENTATION
  ## Description: centre of the birth distribution (between first and last birth).
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: centre of the birth distribution.
  ## Details: none.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  
  # FUNCTION
  first_birth <- which(pattern$nb_births>0)[1] 
  last_birth <- which(pattern$nb_births>0)[length(which(pattern$nb_births>0))] 
  centre <- (first_birth+last_birth)/2 
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    abline(v=c(first_birth, centre, last_birth), col=c("purple", "pink", "purple"))
    text(centre, max(pattern$nb_births)-1, "Birth period centre", col="pink")
  }
  
  # OUTPUT
  return(centre)
}

# MODE, PROP BIRTHS MODE
pheno.jemison <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: time unit gathering the highest number of births and proportion of births during this time unit.
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: time unit gathering the highest number of births, proportion of births during the time unit gathering the highest number of births.
  ## Details: works better with unimodal distributions, otherwise could find more than one mode. Returns a warning when more than one peak is detected.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  
  # FUNCTION
  max_birth <- which(pattern$nb_births==max(pattern$nb_births, na.rm=T))
  if (length(max_birth)>1) { warning("More than one peak found") }
  nb_births <- unique(pattern$nb_births[max_birth])
  prop <- nb_births/sum(pattern$nb_births, na.rm=T)
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    abline(v=max_birth, col="green")
  }
  
  # OUTPUT
  return(list(mode=max_birth, prop_births_mode=prop))
}

# PROP BIRTHS AROUND MEDIAN, PROP BIRTHS AROUND MODE
pheno.adamsjemison <- function(pattern, reference="mode", period=3, graph=F) { 
  # DOCUMENTATION
  ## Description: proportion of births in a given period, centred around the median or the mode birth date.
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # period = specify the length of the period to consider.
  # reference = specify if the period should be centred around the median or the mode.
  ## Values: proportion of births around the median, proportion of birth around the mode.
  ## Details: if "period" is odd, the function looks for the period [reference - floor(period/2) ; reference + floor(period/2)]. If "period" is even, the function
  # looks for the period [reference - period/2 ; reference + period/2]. If there is more than one mode, the function takes the first one as a reference.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  tot_birth <- sum(pattern$nb_births)
  ref <- floor(period/2)
  
  # FUNCTION
  pattern$cum <- cumsum(pattern$nb_births)
  med_birth <- which(pattern$cum>=sum(pattern$nb_births)/2)[1]
  max_birth <- which(pattern$nb_births==max(pattern$nb_births, na.rm=T))[1] # if more than one mode, take the first one
  if (length(max(pattern$nb_births, na.rm=T))>1) { warning("More than one peak found") }
  if (reference=="mode") { central <- max_birth } else if (reference=="median") { central <- med_birth }
  borne_inf <- central-ref # initially, consider that the index does not go beyond the interval [1 ; nrow(pattern)]
  borne_sup <- central+ref
  if (central-ref<=0) { # if the index goes beyond 1
    borne_inf <- 1
  }
  if (central+ref>nrow(pattern)) { # if the index goes beyond nrow(pattern)
    borne_sup <- nrow(pattern)
  }
  birth_int <- sum(pattern$nb_births[borne_inf:borne_sup])
  prop_birth <- birth_int/tot_birth*100
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    abline(v=c(borne_inf, borne_sup), col="purple")
    text(borne_inf, max(pattern$nb_births)-1, "Lower boundary", col="purple")
    text(borne_sup, max(pattern$nb_births)-1, "Upper boundary", col="purple")
  }
  
  # OUTPUT
  return(list(prop_births=prop_birth, nb_births=birth_int, central_temporal_unit=med_birth, low_bound=borne_inf, upp_bound=borne_sup))
}

# PERIOD HDR
pheno.calabrese <- function(pattern, percent=80, graph=F) {
  # DOCUMENTATION
  ## Description: time period when a certain percentage of births happened using high density regions.
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # percent = specify the threshold percentage of births.
  ## Values: length of the time period(s) when a certain percentage of births happened, time period(s) when a certain percentage of births happened.
  ## Details: works well with unimodal distributions (otherwise, can find more than one interval and even individual time units leading to problems when
  # interpreting the output).
  
  # INITIALISATIONS
  library(hdrcde)
  colnames(pattern) <- c("unit", "nb_births")
  warning("Works well with unimodal like distributions only!") 
  
  # FUNCTION
  data_int <- rep(pattern$unit, pattern$nb_births)
  ints <- hdr(data_int, prob=percent, nn=nsim)
  periods <- ints$hdr # return the lower and upper successive boundaries for all the intervals gathering together the proportion of births selected
  # A unique time unit could contain births to take into account in the calculation: in this context the function does not return two boundaries but only the
  # reference of the time unit
  idx <- 1:length(periods)
  if ((length(idx)/2)==floor(length(idx)/2)) { # if even number of boundaries: delimited periods can be found (warning: except if even number of individual time
    odd <- idx%%2==1                           #  units, then the following calculation is wrong)
    evn <- idx%%2==0
    period_length_hdr <- sum(ints$hdr[evn]-ints$hdr[odd]) # sum of all the periods
    x1 <- periods[which(idx%%2==1)] # for graphical representation
    x2 <- periods[which(idx%%2==0)]
  } else { 
    warning(paste("one period containing a part of the", percent, "% of births is a complete temporal unit (meaning uneven number of boundraies) :",
                  "period_length cannot be calculated", sep=" ")) 
    period_length_hdr <- NA
    x1 <- NA
    x2 <- NA
  }
  
  # PLOT
  if (graph==T) {
    if (unique(is.na(x1))==F) { # plot available only if boundaries are identified, otherwise impossible to know which boundary to link together
      plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
      rect(x1, rep(0, length(x1)), x2, rep(max(pattern$nb_births+1), length(x1)), border=T, col=rgb(1, 0, 0, 0.2))
    }
  }
  
  # OUTPUT
  return(list(period_hdr=period_length_hdr, which_periods=periods))
}

# MAX PROP BIRTHS PERIOD GIVEN, MIN PROP BIRTHS PERIOD GIVEN
pheno.owensmith1 <- function(pattern, period=2, object="max", graph=F) { 
  # DOCUMENTATION
  ## Description: proportion of births occuring during the period gathering the most and the least births.
  ## Arguments: 
  # object = specify if the minimum or the maximum birth period should be returned.
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # period = specify the period length to analyse.
  ## Values: proportion of births in the period gathering the least births, proportion of births in the period gathering the most births.
  ## Details: several periods can meet the asumption and overlap.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  tot_births <- sum(pattern$nb_births)
  
  # FUNCTION
  sum_int <- data.frame(indice=seq(1, nrow(pattern)-period+1, 1), tot=NA)
  for (i in 1:(nrow(pattern)-period+1)) {
    sum_int$tot[i] <- sum(pattern$nb_births[i:(i+period-1)], na.rm=T)
  }
  if (object=="max") {
    sum_fin <- max(sum_int$tot, na.rm=T)
  } else if (object=="min") {
    sum_fin <- min(sum_int$tot, na.rm=T)
  }
  prop_b <- sum_fin/tot_births*100
  period_object_birth <- data.frame(beginning=rep(NA, length(which(sum_int$tot==sum_fin))), end=NA)
  for (j in 1:nrow(period_object_birth)) {
    period_object_birth$beginning[j] <- which(sum_int$tot==sum_fin)[j]
    period_object_birth$end[j] <- which(sum_int$tot==sum_fin)[j]+period-1
  }
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    rect(period_object_birth$beginning, rep(0, length(period_object_birth)), period_object_birth$end, 
         rep(max(pattern$nb_births+1), length(period_object_birth)), border=T, col=rgb(1, 0, 0, 0.2))
  }
  
  # OUTPUT
  return(list(prop_births=prop_b, nb_births=sum_fin, period=period_object_birth))
}

# DIFF MIN MAX PROP BIRTHS
pheno.owensmith2 <- function(pattern, period=3, graph=F) { 
  # DOCUMENTATION
  ## Description: difference between the proportion of births occurring during the period gathering the most births and the least births.
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # period = specify the period length to analyse (same for min and max).
  ## Values: difference between the proportion of births occurring during the period gathering the most births and the least births.
  ## Details: no graphical representation available.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  
  # FUNCTION
  owen_max <- pheno.owensmith1(pattern=pattern, period=period, object="max")
  owen_min <- pheno.owensmith1(pattern=pattern, period=period, object="min")
  differ <- owen_max$prop_births-owen_min$prop_births
  
  # PLOT
  if (graph==T) {
    # no representation
  }
  
  # OUTPUT
  return(diff_min_max_prop_births=differ)
}

# MEAN VECTOR LENGTH, MEAN VECTOR ORIENTATION
pheno.campos <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: mean vector orientation and mean vector length of the births distribution, once converted into cicular data (radians).
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births)
  ## Values: mean birth date, mean vector length of the birth distribution, mean vector orientation of the birth distribution.
  ## Details: mean vector orientation corresponds to the mean birth date, mean vector length corresponds to the spread of the birth distribution (closed to 0: the
  # distribution is uniform, closed to 1: the distribution is unimodal).
  
  # INITIALISATIONS
  library(circular)
  library(CircStats)
  colnames(pattern) <- c("unit", "nb_births")
  ref <- nrow(pattern)
  pattern$deg <- seq(0, 359.999, 360/ref)
  pattern$rad <- (pattern$deg*pi)/180
  data_rad <- rep(pattern$rad, pattern$nb_births)
  data_rad <- as.circular(data_rad, units="radians", info=F)
  
  # FUNCTION
  mean_vect_orientation <- mean.circular(data_rad) # in radians
  if (mean_vect_orientation<0) {mean_vect_orientation <- mean_vect_orientation+2*pi} # conversion to get values >0 only
  mean_birth_date <- as.numeric(as.character((((mean_vect_orientation*180/pi)*nrow(pattern))/360)+1)) # estimate mean birth date from mean vector orientation
  mean_vect_length <- est.rho(data_rad)
  
  # PLOT
  if (graph==T) {
    plot.circular(data_rad, stack=T, bins=360, shrink=3) 
    lines.circular(rep(mean_vect_orientation, 2), c(0, mean_vect_length), col="brown", lwd=3)
  }
  
  # OUTPUT
  return(list(mean_vector_orientation=mean_vect_orientation, mean_birth_date=mean_birth_date, mean_vect_length=mean_vect_length))
}

# SKINNER
pheno.skinner <- function(pattern, period=4, percent=75, graph=F) { 
  # DOCUMENTATION
  ## Description: does a given proportion of births occur during a given period or not?
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # percent = specify the threshold percentage of birth.
  # period = specify the period length to analyse.
  ## Values: births synchronous or not.
  ## Details: none.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  nb_births_lim <- percent/100*sum(pattern$nb_births) 
  summ <- data.frame(begin=rep(NA, length(1:(nrow(pattern)-period+1))), end=NA, nb_births=NA)
  
  # FUNCTION
  for (i in 1:(nrow(pattern)-period+1)) {
    test <- sum(pattern$nb_births[i:(i+period-1)], na.rm=T)
    if (test<nb_births_lim) {
      duration <- NA
    } else {
      duration <- paste(i, i+period-1, sep="-") 
    }
    summ$nb_births[i] <- test
    summ$begin[i] <- i
    summ$end[i] <- i+period-1
  }
  summ <- na.omit(summ)
  summ <- summ[order(summ$nb_births, decreasing=T), ]
  summ <- summ[summ$nb_births>nb_births_lim, ] # keep only the periods satifying the criteria
  if (nrow(summ)>=1) { ccl <- TRUE } else { ccl <- FALSE }
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    rect(summ$begin, rep(0, nrow(summ)), summ$end, rep(max(pattern$nb_births+1), nrow(summ)), border=T, col=rgb(1, 0, 0, 0.2))  
  }
  
  # OUTPUT
  return(list(skinner=ccl, summ=summ))
}

# MIN PERIOD PROP BIRTHS GIVEN
pheno.meng <- function(pattern, percent=75, graph=F) { 
  # DOCUMENTATION
  ## Description: shortest period(s) gathering at least a certain percentage of births, if existing (otherwise, the complete birth period is returned).
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # percent = specify the threshold percentage of births.
  ## Values: all the periods gathering at least the threshold percentage of births, the shortest period gathering the highest percentage of births above the
  # threshold.
  ## Details: none.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  tot_births <- sum(pattern$nb_births)
  ref <- percent/100*tot_births
  sum_int <- data.frame(begin=seq(1, nrow(pattern), 1), end=NA, nb_tot_births=NA)
  
  # FUNCTION
  for (i in 1:nrow(pattern)) {
    births_sum <- pattern$nb_births[i]
    j <- 0
    while (j+1<=nrow(pattern)-i & births_sum<ref) {
      j <- j+1
      births_sum <- births_sum+pattern$nb_births[i+j]
    }
    sum_int$nb_tot_births[i] <- births_sum
    sum_int$end[i] <- i+j
  }
  sum_int$length <- sum_int$end-sum_int$begin+1
  sum_int$prop_births <- (sum_int$nb_tot_births/sum(pattern$nb_births, na.rm=T))*100
  sum_tri <- sum_int[which(sum_int$nb_tot_births>ref), ] 
  sum_tri <- sum_tri[order(sum_tri$length), ]
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    rect(sum_tri$begin, rep(0, nrow(sum_tri)), sum_tri$end, rep(max(pattern$nb_births+1), nrow(sum_tri)), border=T, col=rgb(1, 0, 0, 0.2))  
  }
  
  # OUTPUT
  return(list(min_period_prop_births_given=sum_tri[1, ], summ=sum_tri))
}

# ZERBE PERL
pheno.zerbe <- function(pattern, percent=80, graph=F) { 
  # DOCUMENTATION
  ## Description: period length gathering a certain percentage of births around the mode of the distribution. Cf. Zerbe P., Clauss M., Codron D., Bingaman Lackey
  # L., Rensch E., Streich J. W., Hatt J.-M. & Müller D. W. (2012). Reproductive seasonality in captive wild ruminants: implications for biogeographical adaptation,
  # photoperiodic control, and life history. Biological Reviews, 87, 965-990. DOI: 10.1111/j.1469-185x.2012.00238.x
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # percent = specify the threshold percentage of births.
  ## Values: all the periods gathering a certain percentage of births around the mode, the shortest period gathering a certain percentage of births around the mode.
  ## Details: rounds to the smallest integer when a decimal value is returned. See Zerbe et al. 2012 for a description of the function. Require a working directory
  # to work.
  
  # INITIALISATIONS
  setwd(working_dir)
  library(stringr)
  colnames(pattern) <- c("unit", "nb_births")
  c <- ncol(pattern)-1
  n <- nrow(pattern)
  write.table(pattern, file="pattern.txt", sep=";", row.names=F, col.names=F)
  
  # FUNCTION
  cmd <- paste("perl ", "peakfinder.pl", "--t", percent, "--c", c, "--n", n, "pattern.txt", sep=" ")
  object <- system(cmd, intern=T) # call perl
  # extraction of the output of perl
  summ <- data.frame(begin=rep(NA, length(object)), end=NA, nb_births=NA, prop_births=NA)
  summ$begin <- as.numeric(str_extract(object, pattern="(?<=block )\\d+"))+1 # +1 because initialisation at 0
  summ$end <- as.numeric(str_extract(object, pattern="(?<=to )\\d+"))+1
  summ$nb_births <- as.numeric(str_extract(object, pattern="\\d+(?= births)"))
  summ$prop_births <- as.numeric(str_extract(object, pattern="(?<=prop of )\\d+\\.\\d+"))
  summ$length_period <- summ$end-summ$begin+1
  summ <- na.omit(summ)
  test <- unique(str_extract(object, pattern="No window found"))
  if (length(test)>1) {
    summ <- data.frame(begin=NA, end=NA, nb_births=NA, prop_births=NA, length_period=NA)
    warning("No window found")
  }
  
  # PLOT
  if (graph==T) {
    par(mfrow=c(1, 1))
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l")
    rect(summ$begin, rep(0, nrow(summ)), summ$end, rep(max(pattern$nb_births+1), nrow(summ)), border=T, col=rgb(1, 0, 0, 0.2))  
  }
  
  # OUTPUT
  return(list(zerbe_perl=summ[1, ], summ=summ))
}

# RUTBERG
pheno.rutberg <- function(pattern, percent=80, consecutive=F, graph=F) {
  # DOCUMENTATION
  ## Description: period length gathering a certain percentage of births since first birth.
  ## Arguments: 
  # consecutive = specify if time units should be counted until the first birth followed by time units gathering at least one birth only (no time unit without
  # births in between) or not.
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: first birth date, last birth date, number of births in the period, percentage of births in the period, period length gathering a certain percentage of
  # births since first birth.
  ## Details: when consecutive is set to TRUE, the function could not detect any period gathering enough births.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  ref <- percent/100*sum(pattern$nb_births)
  deb <- which(pattern$nb_births>0)[1] # start count when first birth is detected
  fin <- NA
  count <- 0
  
  # FUNCTION
  for (i in deb:nrow(pattern)) {
    count <- count+pattern$nb_births[i]
    if (consecutive==TRUE) {
      if (pattern$nb_births[i]==0) { # if a time unit without births is encountered, back to 0
        count <- 0
        deb <- i+1
      }
    } else if (consecutive==FALSE) {
      count <- count # if consecutive = F and a tile unit without births is encountered, do nothing
    }
    if (count>=ref) { # when the given percentage is met, stop
      fin <- i 
      break;
    }
  }
  if(is.na(fin)==T) { warning("No period with 80% of births found with arg consecutive=TRUE, try consecutive=FALSE") }
  deg_synch <- fin-deb+1
  prop <- (count/sum(pattern$nb_births, na.rm=T))*100
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l") 
    abline(v=c(deb, fin), col="purple")
  }
  
  # OUTPUT
  return(list(rutberg=deg_synch, beginning_birth_period=deb, end_birth_period=fin, nb_births=count, prop_births=prop))
}

# MEDIAN PROBIT, STANDARD DEVIATION PROBIT
pheno.caughley <- function(pattern, graph=F) { 
  # DOCUMENTATION
  ## Description: median birth date and standard error of the birth distribution via probit analysis.
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: median birth date, standard deviation of the birth distribution, standard error of the median birth date.
  ## Details: works with unimodal distributions.
  
  # INITIALISATIONS
  library(ecotoxicology)
  colnames(pattern) <- c("unit", "nb_births")
  pattern$prop_births <- (pattern$nb_births/sum(pattern$nb_births))*100
  pattern$cum <- cumsum(pattern$prop_births)
  
  # FUNCTION
  pattern$probit <- PercentageToProbit(pattern$cum)
  for (i in 1:nrow(pattern)) {
    if (pattern$probit[i] %in% c("Inf", "-Inf")) {
      pattern$probit[i] <- NA
    }
  }
  pattern$probit_weight <- Probitw(pattern$probit)
  mod <- lm(data=pattern, probit~unit)
  recap <- mod
  a <- recap$coefficients[1]
  b <- recap$coefficients[2]
  m <- (5-a)/b # median
  sd <- 1/b # standard deviation
  
  # PLOT
  if (graph==T) {
    plot(pattern$probit~pattern$unit, col="black", pch=16, xlab="unit", ylab="probit of proportion of births")
    lines(mod$fitted.values~mod$model$unit, col="red", pch=16, cex=1.5, xlab="unit", ylab="cumulated data")
    abline(v=m, col="green")
  }
  
  # OUTPUT
  return(list(median_probit=m, standard_deviation_probit=sd))
}

# PERIOD GAUSSIAN
pheno.paoli <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: duration of the birth period using the formula 4*standard error.
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: duration of the birth period.
  ## Details: works with unimodal distributions.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  data_int <- rep(pattern$unit, pattern$nb_births)
  
  # FUNCTION
  peak <- mean(data_int, na.rm=T)
  sd <- sd(data_int, na.rm=T)
  period_length <- 2*2*sd 
  period_gaussian <- period_length
  
  # PLOT
  if (graph==T) {
    hist(data_int, xlim=c(0, nrow(pattern)))
    hist(rnorm(length(data_int), mean=peak, sd=sd), add=T, col=rgb(1, 0, 0, 0.2))
    abline(v=peak, lwd=5, col="green")
    lines(c((peak-period_length/2), (peak+period_length/2)), c(0, 0), lwd=5, col="cyan")
  }
  
  # OUTPUT
  return(period_gaussian)
}

# RAYLEIGH
pheno.dibitetti <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: is the distribution unimodal? rayleigh test (H0: random distribution, H1: unimodal distribution).
  ## Arguments: 
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  ## Values: is the distribution unimodal or not, mean resultant length, p-value.
  ## Details: none.
  
  # INITIALISATIONS
  library(circular)
  library(CircStats)
  colnames(pattern) <- c("unit", "nb_births")
  ref <- nrow(pattern) 
  pattern$deg <- seq(0, 359.999, 360/ref)
  pattern$rad <- (pattern$deg*pi)/180
  data_rad <- rep(pattern$rad, pattern$nb_births)
  data_rad_origin <- as.circular(data_rad, units="radians", info=F)
  
  # FUNCTION
  rayleigh_test <- r.test(data_rad)
  if (rayleigh_test$p.value>0.05) {
    signif_seasonal <- F 
  } else {
    signif_seasonal <- T  
  }
  
  # PLOT
  if (graph==T) {
    plot.circular(data_rad_origin, stack=T, bins=360, shrink=3) 
  }
  
  # OUTPUT
  return(list(rayleigh=signif_seasonal, p_value=rayleigh_test$p.value, mean_resultant_length=rayleigh_test$r.bar))  
}

# KOLMOGOROV SMIRNOV GAUSSIAN, KOLMOGOROV SMIRNOV UNIFORM
pheno.riedmanmeng <- function(pattern, ref="uniform", origin="corr", graph=F) { 
  # DOCUMENTATION
  ## Description: is the birth distribution similar to a normal or a uniform distribution? Kolmogorov Smirnov test (H0: the distributions are similar).
  ## Arguments: 
  # origin = specify if the origin of the dataset should be the first birth date or the first time unit of the cycle.
  # pattern = a data frame with two columns (col1=time unit, col2=number of births).
  # ref = specify the distribution of reference to compare to the birth distribution.
  ## Values: is the birth distribution similar to the distribution of reference, distribution used as a reference, p value, statistic of the test.
  ## Details: data should not contain any ex aequo. Using the parameter "origin" allows to get better results when comparing birth distriubtions with high
  # resolution. No graphical representation available.
  
  # INITIALISATIONS
  colnames(pattern) <- c("unit", "nb_births")
  if (origin=="corr") {
    first_birth <- which(pattern$nb_births>0)[1] 
    last_birth <- which(pattern$nb_births>0)[length(which(pattern$nb_births>0))] 
    pattern <- pattern[first_birth:last_birth, ]
  }
  
  # FUNCTION
  if (ref=="uniform") {
    kolmo <- ks.test(pattern$nb_births, "punif", min=0, max=max(pattern$nb_births, na.rm=T), alternative="two.sided")
    val_ref <- rep(1/nrow(pattern), nrow(pattern))*sum(pattern$nb_births, na.rm=T)
  } else if (ref=="gaussian") {
    kolmo <- ks.test(pattern$nb_births, "pnorm", mean=mean(pattern$nb_births, na.rm=T), sd=sd(pattern$nb_births, na.rm=T), alternative="two.sided")
    val_ref <- rnorm(1:length(pattern$nb_births), mean=mean(pattern$nb_births, na.rm=T), sd=sd(pattern$nb_births, na.rm=T))
  } 
  if (kolmo$p.value>0.05) {
    signif_distrib <- T 
  } else {
    signif_distrib <- F  
  }
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births~pattern$unit, xlab="Time", ylab="Nb births", type="l") 
    if (ref=="uniform") {
      lines(val_ref~pattern$unit, col="red")
    } else if (ref=="gaussian") {
      # no representation
    }
  }
  
  # OUTPUT
  return(list(distrib_tested=ref, kolmogorov_smirnov=signif_distrib, p_value=kolmo$p.value, statistic_d=kolmo$statistic))
}

# BEGINNING PERIOD MEAN MULTI, MEAN MEAN MULTI, MEDIAN MEAN MULTI, PERIOD MEAN MULTI
pheno.millar <- function(pattern, graph=F) { 
  # DOCUMENTATION
  ## Description: mean date of the start of the birth periods, end of the birth periods, duration of the birth periods, mean birth dates, median birth dates, and
  # their standard deviation.
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: mean end time unit, mean mean birth time unit, mean median birth date, mean period length,  mean start time unit, standard deviation of the end time
  # units, standard deviation of the mean birth dates, standard deviation of the median birth dates, standard deviation of the period length, standard deviation of
  # the start time units.
  ## Details: none.
  
  # INITIALISATIONS
  library(zoo)
  library(tidyr)
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle)
  first_birth <- NULL
  last_birth <- NULL
  mean_birth <- NULL
  median_birth <- NULL
  period <- NULL
  j <- 0
  
  # FUNCTION
  for (i in all_cycles) {
    j <- j+1
    first_birth[j] <- which(pattern$nb_births[pattern$cycle==i]>0)[1] 
    last_birth[j] <- which(pattern$nb_births[pattern$cycle==i]>0)[length(which(pattern$nb_births[pattern$cycle==i]>0))] 
    mean_birth[j] <- pheno.findlaylambin(pattern[pattern$cycle==i, c("unit", "nb_births")])$mean
    median_birth[j] <- pheno.findlaylambin(pattern[pattern$cycle==i, c("unit", "nb_births")])$median
    period[j] <- last_birth[j]-first_birth[j]+1
  }
  mean_first_birth <- mean(first_birth, na.rm=T)
  sd_first_birth <- sd(first_birth, na.rm=T)
  mean_last_birth <- mean(last_birth, na.rm=T)
  sd_last_birth <- sd(last_birth, na.rm=T)
  mean_mean_birth <- mean(mean_birth, na.rm=T)
  sd_mean_birth <- sd(mean_birth, na.rm=T)
  mean_median_birth <- mean(median_birth, na.rm=T)
  sd_median_birth <- sd(median_birth, na.rm=T)
  mean_birth_period_length <- mean(period, na.rm=T)
  sd_birth_period_length <- sd(period, na.rm=T)
  
  # PLOT
  if (graph==T) {
    patternbis <- as.data.frame(spread(data=pattern, key="cycle", value="nb_births"))
    data_plot <- read.zoo(patternbis)
    pnl <- function(x, y, ...) {
      panel.number <- parent.frame()$panel.number
      abline(v=c(mean_first_birth, mean_last_birth, mean_mean_birth, mean_median_birth), col=c("purple", "purple", "brown", "blue"))
      rect(mean_first_birth-sd_first_birth, 0, mean_first_birth+sd_first_birth, max(pattern$nb_births), col=rgb(0.6, 0, 1, 0.2), border=F)
      rect(mean_last_birth-sd_last_birth, 0, mean_last_birth+sd_last_birth, max(pattern$nb_births), col=rgb(0.6, 0, 1, 0.2), border=F)
      rect(mean_mean_birth-sd_mean_birth, 0, mean_mean_birth+sd_mean_birth, max(pattern$nb_births), col=rgb(0.8, 0.2, 0, 0.2), border=F)
      rect(mean_median_birth-sd_median_birth, 0, mean_median_birth+sd_median_birth, max(pattern$nb_births), col=rgb(0, 0, 1, 0.2), border=F)
      lines(x, y)
    }
    plot(data_plot, panel=pnl)
  }
  
  # OUTPUT
  return(list(beginning_birth_season=list(mean_multi=mean_first_birth, standard_deviation_multi=sd_first_birth), 
              end_birth_season=list(mean_multi=mean_last_birth, standard_deviation_multi=sd_last_birth), 
              period=list(period_mean_multi=mean_birth_period_length, period_standard_deviation_multi=sd_birth_period_length),
              mean=list(mean_multi=mean_mean_birth, standard_deviation_multi=sd_mean_birth),
              median=list(mean_multi=mean_median_birth, standard_deviation_multi=sd_median_birth)))
}

# COMP MEAN CI
pheno.whiting <- function(pattern, CI=95, graph=F) { 
  # DOCUMENTATION
  ## Description: are the two mean birth dates similar? by checking if the confidence intervals overlap.
  ## Arguments: 
  # CI = specifiy the percentage of the confidence intervals.
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: are the two confidence interval overlaping or not, the two confidence intervals, the two mean birth dates to compare.
  ## Details: none.
  
  # INITIALISATIONS
  library(zoo)
  library(tidyr)
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle) 
  if(length(all_cycles)>2) { warning(paste("More than two cycles detected, the function works with two cycles at once", 
                                           "(the function is applied on the two first cycles only).", sep=" "))}
  demi_alpha <- -(CI/100-1)/2
  ref <- 1-demi_alpha
  pattern.list <- vector("list", 2)
  for (i in 1:2) {
    pattern.list[[i]] <- pattern[pattern$cycle==all_cycles[i], c("unit", "nb_births")]
  }
  names(pattern.list) <- all_cycles[1:2]
  sub_list <- vector("list", 2)
  
  # FUNCTION
  for (i in 1:2) { # find start, end, mean and median birth date and the birth period for both distributions
    data_int <- data.frame(mean=NA, sd=NA, low=NA, up=NA) 
    data_int$mean <- pheno.findlaylambin(pattern.list[[i]])$mean
    data_int$sd <- pheno.findlaylambin(pattern.list[[i]])$standard_deviation
    data_int$up <- data_int$mean+(qnorm(ref)*(data_int$sd/sqrt(sum(pattern.list[[i]]$nb_births))))
    data_int$low <- data_int$mean-(qnorm(ref)*(data_int$sd/sqrt(sum(pattern.list[[i]]$nb_births))))
    if (data_int$low<0) {
      data_int$low <- 0
    }
    if (data_int$up>nrow(pattern.list[[i]])) {
      data_int$up <- nrow(pattern.list[[i]])
    }
    sub_list[[i]] <- data_int
  } 
  if (sub_list[[2]]$low<=sub_list[[1]]$low & sub_list[[2]]$up>=sub_list[[1]]$up) { # are confidence intervals overlaping or not
    seasonal <- TRUE
  } else if (sub_list[[2]]$low>=sub_list[[1]]$low & sub_list[[2]]$up<=sub_list[[1]]$up) {
    seasonal <- TRUE
  } else if (sub_list[[2]]$low<=sub_list[[1]]$low & sub_list[[2]]$up<=sub_list[[1]]$up & sub_list[[2]]$up>=sub_list[[1]]$low) {
    seasonal <- TRUE
  } else if (sub_list[[2]]$low>=sub_list[[1]]$low & sub_list[[2]]$low<=sub_list[[1]]$up & sub_list[[2]]$up>=sub_list[[1]]$up) {
    seasonal <- TRUE
  } else {
    seasonal <- FALSE
  }
  
  # PLOT
  if (graph==T) {
    pattern <- pattern[ , c("cycle", "unit", "nb_births")]
    patternbis <- as.data.frame(spread(data=pattern, key="cycle", value="nb_births"))
    data_plot <- read.zoo(patternbis)[ , 1:2]
    pnl <- function(x, y, ...) {
      panel.number <- parent.frame()$panel.number
      for (i in 1:2) {
        if (panel.number==i) { # for the first birth distribution
          abline(v=sub_list[[i]]$mean, col="brown")
          rect(sub_list[[i]]$low, 0, sub_list[[i]]$up, max(pattern$nb_births), col=rgb(0.8, 0.2, 0, 0.2), border=F)
        }
      }
      lines(x, y)
    }
    plot(data_plot, panel=pnl)
  }
  
  # OUTPUT
  return(list(comp_mean_ci=seasonal, values_per_cycle=sub_list))
}

# WATSON WILLIAMS
pheno.pare <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: are the mean birth dates similar between the two birth distributions. Watson-Williams test (H0 : the means are similar).
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: are the mean similar or not, statistics of the test, p-value of the test.
  ## Details: the two distributions must be drawn from a von Mises distribution, the concentration parameters have the same value in both distributions and are
  # large enough (>1) (checked in the test).
  
  # INITIALISATIONS
  library(circular)
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle)
  pattern$deg <- NA
  pattern$rad <- NA
  data_rad <- vector("list", 2)
  if(length(all_cycles)>2) { warning(paste("More than two cycles detected, the function works with two cycles at once", 
                                           "(the function is applied on the two first cycles only).", sep=" "))}
  
  # FUNCTION
  for (i in 1:2) {
    ref <- nrow(pattern[pattern$cycle %in% all_cycles[i], ])
    pattern$deg[pattern$cycle %in% all_cycles[i]] <- seq(0, 359.999, 360/ref)
    pattern$rad[pattern$cycle %in% all_cycles[i]] <- (pattern$deg[pattern$cycle %in% all_cycles[i]]*pi)/180
    data_rad[[i]] <- rep(pattern$rad[pattern$cycle %in% all_cycles[i]], pattern$nb_births[pattern$cycle %in% all_cycles[i]])
    data_rad[[i]] <- as.circular(data_rad[[i]], units="radians", info=F)
  }
  watson_test <- watson.williams.test(data_rad)
  if (watson_test$p.value<0.05) {
    signif <- F 
  } else {
    signif <- T
  }
  mean1 <- mean.circular(data_rad[[1]])
  mean2 <- mean.circular(data_rad[[2]])
  
  # PLOT
  if (graph==T) {
    par(mfrow=c(1, 2))
    plot.circular(data_rad[[1]], stack=T, bins=360, shrink=3) 
    points.circular(mean1, col="brown", lwd=3)
    plot.circular(data_rad[[2]], stack=T, bins=360, shrink=3) 
    points.circular(mean2, col="brown", lwd=3)
  }
  
  # OUTPUT
  return(list(watson_williams=signif, p_value=watson_test$p.value, statistic_f=watson_test$statistic))  
}

# COMP PEAK SIGMOID CI, PEAK SIGMOID
pheno.moe2 <- function(pattern, nb_cycles=2, CI=95, resolution=0.1, graph=F) {
  # DOCUMENTATION
  ## Description: when one reproduction cycle is given: peak birth date using a logistic regression on the birth distribution. When two reproductive cycles are
  # given: evaluate the level of consistency of the peak birth date.
  ## Arguments: 
  # CI = specify the confidence interval.
  # nb_cycles = specify the number of cycles (one or two).
  # pattern = a data frame with two columns (col1=time unit, col2=number of births) or three columns (col1=first time unit, col2=second time unit, col3=number of
  # births).
  # resolution = specify the size of the steps between each point for the prediction.
  ## Values: consitency of the peak birth date (if two cycles provided), peak birth date.
  ## Details: works for unimodal distributions. Based on the equation A/(1+c*exp(-k*t)) with A the asymptote, c the integration constant, k the growth rate.
  
  # INITIALISATIONS
  library(drc)
  seasonal <- NA
  if (nb_cycles==1) {
    colnames(pattern) <- c("unit", "nb_births")
    pattern$cum <- cumsum(pattern$nb_births)
  } else if (nb_cycles==2) {
    list_int <- vector("list", 2)
    colnames(pattern) <- c("cycle", "unit", "nb_births")
    all_cycles <- unique(pattern$cycle)
    for (i in 1:2) {
      pattern$cum[pattern$cycle==all_cycles[i]] <- cumsum(pattern$nb_births[pattern$cycle==all_cycles[i]])
    }
  } else if (nb_cycles<1 | nb_cycles>2) {
    warning(paste("More than two cycles given, the function works with one or two cycles at once", 
                  "(the function is applied on the two first cycles only).", sep=" "))
  }
  
  # FUNCTION  
  if (nb_cycles==1) { # if only one cycle provided, find birth peak date
    mod <- drm(cum~unit, data=pattern, fct=L.3(), type="continuous")
    summod <- summary(mod) # b=lower asymptote, d=upper asymptote, e=inflection point
    pattern$fitthreshold <- predict(mod, type="response") # predictions
    pattern$ci_lo <- predict(mod, interval="confidence")[ , 2] # lower confidence interval
    pattern$ci_up <- predict(mod, interval="confidence")[ , 3] # upper confidence interval
    peak_sigmoid <- summod$coefficients[3]
  } else if (nb_cycles==2) { # if two cycles provided, find both birth peak dates and compare them
    for (i in 1:2) {
      mod <- drm(cum~unit, data=pattern[pattern$cycle==all_cycles[i], ], fct=L.3(), type="continuous")
      list_int[[i]]$mod <- mod
      summod <- summary(mod) # b=lower asymptote, d=upper asymptote, e=inflection point
      list_int[[i]]$peak_sigmoid <- summod$coefficients[3]
      list_int[[i]]$pattern <- data.frame(unit=seq(1, max(unique(pattern$unit[pattern$cycle==all_cycles[i]])), resolution)) # prepare data frame for predictions
      list_int[[i]]$pattern$fitthreshold <- predict(mod, newdata=list_int[[i]]$pattern, type="response") # predictions
      list_int[[i]]$pattern$ci_lo <- predict(mod, newdata=list_int[[i]]$pattern, interval="confidence")[ , 2] # lower confidence interval
      list_int[[i]]$pattern$ci_up <- predict(mod, newdata=list_int[[i]]$pattern, interval="confidence")[ , 3] # upper confidence interval
      list_int[[i]]$low <- confint(mod, level=CI/100)[3, 1] # confidence intervals on the birth peak dates
      list_int[[i]]$up <- confint(mod, level=CI/100)[3, 2]
    }
    peak_sigmoid <- c(list_int[[1]]$peak_sigmoid, list_int[[2]]$peak_sigmoid) # birth peak dates
    if (list_int[[2]]$low<=list_int[[1]]$low & list_int[[2]]$up>=list_int[[1]]$up) { # does confidence interval overlap?
      seasonal <- TRUE
    } else if (list_int[[2]]$low>=list_int[[1]]$low & list_int[[2]]$up<=list_int[[1]]$up) {
      seasonal <- TRUE
    } else if (list_int[[2]]$low<=list_int[[1]]$low & list_int[[2]]$up<=list_int[[1]]$up & list_int[[2]]$up>=list_int[[1]]$low) {
      seasonal <- TRUE
    } else if (list_int[[2]]$low>=list_int[[1]]$low & list_int[[2]]$low<=list_int[[1]]$up & list_int[[2]]$up>=list_int[[1]]$up) {
      seasonal <- TRUE
    } else {
      seasonal <- FALSE
    }
  }
  
  # PLOT
  if (graph==T) {
    if (nb_cycles==1) {
      par(mfrow=c(1, 1))
      plot(pattern$cum~pattern$unit, pch=16, cex=1.5, xlab="unit", ylab="cumulated data")
      points(pattern$fitthreshold~pattern$unit, col="red", type="l", lwd=2)
      abline(v=peak_sigmoid, col="green")
    } else {
      par(mfrow=c(2, 1))
      for (i in 1:2) {
        plot(pattern$cum[pattern$cycle==all_cycles[i]]~pattern$unit[pattern$cycle==all_cycles[i]], pch=16, cex=1.5, xlab="unit", 
             ylab="cumulated data")
        points(list_int[[i]]$pattern$fitthreshold~list_int[[i]]$pattern$unit, col="red", type="l", lwd=2)
        abline(v=list_int[[i]]$peak_sigmoid, col="green")
        abline(v=c(list_int[[i]]$low, list_int[[i]]$up), col="orange")
      }
    }
  }
  
  # OUTPUT
  return(list(peak_sigmoid=peak_sigmoid, comp_peak_sigmoid_ci=seasonal))
}

# MOOD
pheno.berger <- function(pattern, graph=F) { 
  # DOCUMENTATION
  ## Description: are the median birth date similar between the two cycles? Mood test (H0: the medians are identical).
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: are the two medians equals, p value, statistic of the test.
  ## Details: none.
  
  # INITIALISATIONS
  library(coin)
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle) 
  if(length(all_cycles)>2) { warning(paste("More than two cycles detected, the function works with two cycles at once", 
                                           "(the function is applied on the two first cycles only).", sep=" "))}
  sub_list <- vector("list", 2)
  
  # FUNCTION
  sub_df1 <- rep(pattern$unit[pattern$cycle==1], pattern$nb_births[pattern$cycle==1])
  sub_df2 <- rep(pattern$unit[pattern$cycle==2], pattern$nb_births[pattern$cycle==2])
  sub_df <- data.frame(cycle=c(rep(1, length(sub_df1)), rep(2, length(sub_df2))), unit=c(sub_df1, sub_df2))
  sub_df$cycle <- as.factor(sub_df$cycle)
  mood <- median_test(formula=unit~cycle, data=sub_df)
  pvalue <- pvalue(mood)
  if (pvalue<0.05) {
    signif <- F
  } else {
    signif <- T
  }
  
  # PLOT
  if (graph==T) {
    par(mfrow=c(2, 1))
    plot(pattern$nb_births[pattern$cycle==all_cycles[1]]~pattern$unit[pattern$cycle==all_cycles[1]], xlab="Time", ylab="Nb births", type="l")
    abline(v=pheno.findlaylambin(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$median, col="blue")
    text(pheno.findlaylambin(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$median, 
         max(pattern$nb_births[pattern$cycle==all_cycles[1]])-1, "Median birth date", col="blue")
    plot(pattern$nb_births[pattern$cycle==all_cycles[2]]~pattern$unit[pattern$cycle==all_cycles[2]], xlab="Time", ylab="Nb births", type="l")
    abline(v=pheno.findlaylambin(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$median, col="blue")    
    text(pheno.findlaylambin(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$median, 
         max(pattern$nb_births[pattern$cycle==all_cycles[2]])-1, "Median birth date", col="blue")
  }
  
  return(list(mood=signif, p_value=pvalue, statistic_z=mood@statistic@teststatistic))  
}

# COMP MEAN ANOVA
pheno.linnell <- function(pattern, post_hoc=T, graph=F) {
  # DOCUMENTATION
  ## Description: are the mean birth date similar between the two cycles? One way anova test.
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  # post_hoc = specify if a post hoc test should be processed.
  ## Values: are the two means equal, output of the post hoc test, p value, statistic of the test.
  ## Details: post hoc test is a mutlicomparison of Tuckey. Conditions of application: normal distributions of births, homoscedasticity, random and independant
  # samples.
  
  # INITIALISATIONS
  library(agricolae)
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  patternbis <- data.frame(cycle=rep(pattern$cycle, pattern$nb_births), unit=rep(pattern$unit, pattern$nb_births))
  comp_moy2 <- NA
  
  # FUNCTION
  results <- aov(unit~cycle, patternbis)
  summ <- summary(results)
  pvalue <- summ[[1]]$`Pr(>F)`[1]
  if (post_hoc==T) { # should be applied only if the anova is significant
    if (pvalue<0.05) {
      comp_moy <- HSD.test(results, "cycle", "unit", console=F)
      comp_moy2 <- comp_moy$groups
    } else if (pvalue>=0.05) {
      warning("No significant difference according to the ANOVA, the post hoc test is useless.")
    }
  } 
  if (pvalue>0.05) {
    signif_seasonality <- T
  } else {
    signif_seasonality <- F 
  }
  
  # PLOT
  if (graph==T) {
    boxplot(patternbis$unit~patternbis$cycle, xlab="Cycle", ylab="Mean birth unit")
  }
  
  # OUTPUT
  return(list(comp_mean_anova=signif_seasonality, p_value=pvalue, statistic_f=summ[[1]]$`F value`[1], post_hoc=comp_moy2))
}

# DIFF MEAN LINEAR, MEAN LINEAR RANDOM, SEASONALITY LINEAR RANDOM
pheno.loe <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: mean birth date and inter-cycle variance of the birth distributions based on the predictions of a linear model with random effects. The function
  # also estimates the evolution of mean birth date through time using a linear regression.
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: coefficient of determination of the relationship between mean birth date and cycle, confidence interval of the mean birth date, correlation between
  # mean birth date and the cycle, mean birth date, slope coefficient of the relationship between mean birth date and the cycle.
  ## Details: if there is no trend in the relationship between meand birth date and the cycle, the descriptors of this relationship are not meaningful.
  
  # INITIALISATION
  library(lme4)
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle)
  data_int <- vector("list", length(all_cycles))
  
  for (i in 1:length(all_cycles)) {
    data_int[[i]] <- data.frame(cycle=rep(all_cycles[i], sum(pattern$nb_births[pattern$cycle==all_cycles[i]])), birth_date=NA)
    data_int[[i]]$birth_date <- rep(pattern$unit[pattern$cycle==all_cycles[i]], pattern$nb_births[pattern$cycle==all_cycles[i]])
  }
  patternbis <- do.call(rbind, data_int)
  
  # FUNCTION 
  # mean birth date and variance of all the cycles
  mod <- lmer(birth_date~1+(1|cycle), data=patternbis) 
  summ <- summary(mod)
  mean_birth_date <- summ$coefficients[1]
  conf_int <- summ$coefficients[2] 
  seaso <- summ$varcor$cycle[1] # between cycle variance
  # relationship between annual mean birth date and time
  patternbis$pred <- predict(mod, type="response")
  mod2 <- lm(pred~cycle, data=patternbis)
  rsquared <- summary(mod2)$r.squared
  patternbis$pred2 <- predict(mod2, type="response") 
  evolution <- mod2$coefficients[2] # slope coefficient informs about the evolution of mean birth date across time
  
  # PLOT
  if (graph==T) {  
    par(mfrow=c(1, 2))
    plot(pattern$nb_births~pattern$unit, type="l", xlab="TU", ylab="Nb births")
    abline(v=mean_birth_date, col="brown")
    
    plot(patternbis$pred~patternbis$cycle, xlab="Cycle", ylab="Mean nb of births")
    lines(patternbis$pred2~patternbis$cycle, col="red")
  }
  
  # OUTPUT
  return(list(mean_linear_random=mean_birth_date, CI_mean_linear_random=conf_int, seasonality_linear_random=seaso, diff_mean_linear=evolution,
              rsquared_diff_mean_linear=rsquared))
}

# BARTLETT
pheno.hass <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: are the birth distribution variances similar between the cycles? Bartlett test (H0: the variances are the same).
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: are the variance the same, p value, statistic of the test.
  ## Details: birth distribution should follow a normal distribution.
  
  # INITIALISATIONS
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle)
  patternbis <- data.frame(cycle=rep(pattern$cycle, pattern$nb_births), unit=rep(pattern$unit, pattern$nb_births))
  
  # FUNCTION
  bartlett <- bartlett.test(patternbis$unit~patternbis$cycle)
  pvalue <- bartlett$p.value
  if (pvalue>0.05) {
    signif <- T
  } else {
    signif <- F 
  }
  
  # PLOT
  if (graph==T) {
    all_cycles <- unique(pattern$cycle)
    plot(pattern$nb_births[pattern$cycle==all_cycles[1]]~pattern$unit[pattern$cycle==all_cycles[1]], xlab="unit", ylab="Nb births", type="l", 
         col="red")
    for (i in 2:length(all_cycles)) {
      lines(pattern$nb_births[pattern$cycle==all_cycles[i]]~pattern$unit[pattern$cycle==all_cycles[i]], col=(i+1))
    }
    legend("topright", legend=all_cycles, col=c(2, seq(3, length(all_cycles)+1)), lty=1)
  }
  
  # OUTPUT
  return(list(bartlett=signif, p_value=pvalue, k_squared=bartlett$statistic))
}

# SLOPE COMPARISON
pheno.bowyer <- function(pattern, transformation="log", graph=F) {
  # DOCUMENTATION
  ## Description: are the slope coefficient similar between the cycles? using the predictions of a linear model and t test (H0: the medians are identical).
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: are the slopes similar, couples of slope coefficients that are different.
  ## Details: works with unimodal distributions. Warning: several preliminary transformations are done and may not be adapted to any dataset.
  
  # INITIALISATIONS
  library(lsmeans) 
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle)
  unequal <- NULL
  pattern$unit2 <- NA
  for (i in all_cycles) {
    # add a second temporal unit initializing the original one at first birth date and ending at last birth date for each cycle
    ref <- data.frame(original=seq(pattern$unit[pattern$nb_births>0 & pattern$cycle==i][1], max(pattern$unit[pattern$cycle==i], na.rm=T), 1))
    ref$new <- seq(1, nrow(ref), 1)
    for (j in 1:nrow(ref)) {
      pattern$unit2[pattern$unit==ref$original[j] & pattern$cycle==i] <- ref$new[j]
    }
    # calculate the proportion of births and the cumulated proportion of births for each cycle
    pattern$prop_births[pattern$cycle==i] <- (pattern$nb_births[pattern$cycle==i]/sum(pattern$nb_births[pattern$cycle==i], na.rm=T))*100
    pattern$cum[pattern$cycle==i] <- cumsum(pattern$prop_births[pattern$cycle==i]) 
  }
  if (transformation=="log") { # log transformation of the cumulated proportions
    pattern$transfo <- log(pattern$cum)
  } else if (transformation=="none") {
    pattern$transfo <- pattern$cum
  }
  # removal of the dates before first and after last birth date, and of the dates without any birth
  if (length(pattern$cycle[-(which(is.na(pattern$unit2)==T))])>0) {  
    patternbis <- pattern[-(which(is.na(pattern$unit2)==T)), ]
  } else {
    patternbis <- pattern
  }
  patternbis <- patternbis[-(which(patternbis$prop_births==0)), ]
  patternbis$cycle <- as.factor(patternbis$cycle)
  patternbis$unit3 <- log(patternbis$unit2) # log transformation of the dates
  
  # FUNCTION
  model <- lm(transfo~unit3*cycle, data=patternbis)
  # comparison of the slopes of each cycle
  slopes <- lstrends(model, "cycle", var="unit3")
  comp <- as.data.frame(pairs(slopes))
  for (i in 1:nrow(comp)) {
    if (comp$p.value[i]<0.05) {
      unequal <- c(unequal, as.character(comp$contrast[i]))
    }
  }
  if (is.null(unequal)==T) {
    repeat_synchrony <- T
  } else {
    repeat_synchrony <- F
  }
  
  # PLOT
  if (graph==T) {
    patternbis$pred <- predict(model, type="response")
    plot(patternbis$transfo[patternbis$cycle==all_cycles[1]]~patternbis$unit3[patternbis$cycle==all_cycles[1]], pch=16, xlab="log(unit)", 
         ylab="Cumulated percentage of births", xlim=c(min(patternbis$unit3, na.rm=T), max(patternbis$unit3, na.rm=T)), 
         ylim=c(min(patternbis$transfo, na.rm=T), max(patternbis$transfo, na.rm=T)))
    lines(patternbis$pred[patternbis$cycle==all_cycles[1]]~patternbis$unit3[patternbis$cycle==all_cycles[1]], lwd=2)
    for (i in 2:length(all_cycles)) {
      points(patternbis$transfo[patternbis$cycle==all_cycles[i]]~patternbis$unit3[patternbis$cycle==all_cycles[i]], pch=16, col=i+1)
      lines(patternbis$pred[patternbis$cycle==all_cycles[i]]~patternbis$unit3[patternbis$cycle==all_cycles[i]], col=i+1, lwd=2)
    }
    legend("topright", legend=all_cycles, col=c(1, seq(3, length(all_cycles)+1)), lty=1)
  }
  
  # OUTPUT
  return(list(slope_comparison=repeat_synchrony, cycles_different=unequal))
}

# DIFF BEGINNING PERIOD, DIFF FREQUENCY, DIFF MEAN, DIFF MEDIAN, DIFF PEAK, DIFF PERIOD
pheno.diff <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: absolute period length between 2 descriptors of birth distribution (first birth date, peak birth date, mean birth date, median birth date, period
  # length, frequency of the birth distribution)
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: difference between first birth dates, difference between frequency of the birth distribution, difference between mean birth dates, difference between
  # median birth dates, difference between peak birth dates, difference between period length of the birth distribution.
  ## Details: none.
  
  # INITIALISATIONS
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle) 
  if(length(all_cycles)>2) { warning(paste("More than two cycles detected, the function works with two cycles at once", 
                                           "(the function is applied on the two first cycles only).", sep=" "))}
  tab_int <- data.frame(diff_begin=NA, diff_peak=NA, diff_period=NA, diff_freq=NA, diff_mean=NA, diff_median=NA)
  
  # FUNCTION
  tab_int$diff_begin <- abs(pheno.bunnellyomtov(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$beginning_period-
                              pheno.bunnellyomtov(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$beginning_period)
  tab_int$diff_peak <- abs(pheno.jemison(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$mode[1]-
                             pheno.jemison(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$mode[1])
  tab_int$diff_period <- abs(pheno.bunnellyomtov(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$period-
                               pheno.bunnellyomtov(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$period)
  tab_int$diff_freq <- abs(pheno.bunnellyomtov(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$frequency-
                             pheno.bunnellyomtov(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$frequency)
  tab_int$diff_mean <- abs(pheno.findlaylambin(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$mean-
                             pheno.findlaylambin(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$mean)
  tab_int$diff_median <- abs(pheno.findlaylambin(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$median-
                               pheno.findlaylambin(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$median)
  
  # PLOT
  if (graph==T) {
    par(mfrow=c(2, 1))
    plot(pattern$nb_births[pattern$cycle==all_cycles[1]]~pattern$unit[pattern$cycle==all_cycles[1]], xlab="Time", ylab="Nb births", type="l")
    abline(v=c(pheno.bunnellyomtov(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$beginning_period, 
               pheno.jemison(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$mode[1],
               pheno.findlaylambin(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$mean,
               pheno.findlaylambin(pattern[pattern$cycle==all_cycles[1], c("unit", "nb_births")])$median), col=c("purple", "green", "brown", 
                                                                                                                 "blue"))
    plot(pattern$nb_births[pattern$cycle==all_cycles[2]]~pattern$unit[pattern$cycle==all_cycles[2]], xlab="Time", ylab="Nb births", type="l")
    abline(v=c(pheno.bunnellyomtov(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$beginning_period, 
               pheno.jemison(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$mode[1],
               pheno.findlaylambin(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$mean,
               pheno.findlaylambin(pattern[pattern$cycle==all_cycles[2], c("unit", "nb_births")])$median), col=c("purple", "green", "brown", 
                                                                                                                 "blue"))
  }
  
  # OUTPUT
  return(list(diff_beginning_period=tab_int$diff_begin, diff_peak=tab_int$diff_peak, diff_period=tab_int$diff_period,
              diff_frequency=tab_int$diff_freq, diff_mean=tab_int$diff_mean, diff_median=tab_int$diff_median))
}

# KHI2 PROP BIRTHS MEDIAN
pheno.adams <- function(pattern, period=3, graph=F) {
  # DOCUMENTATION
  ## Description: are the proportions of births around median birth date similar between the cycles? Khi 2 test of independance (H0: the proportion of births around
  # median birth date is independant from the cycle).
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  # period = specify the period to time to consider around median birth date to caluclate the proportion of births.
  ## Values: are the proportions of births similar or not, p value, statistic of the test.
  ## Details: tests if the distribution of births proportions according to the cycle follows a uniform distribution, but method not detailed in the reference
  # article.
  
  # INITIALISATIONS
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle)
  mat.med <- matrix(data=NA, nrow=2, ncol=length(all_cycles))
  
  # FUNCTION
  for (i in 1:length(all_cycles)) {
    # proportion of births around the median birth date
    mat.med[1, i] <- pheno.adamsjemison(pattern[pattern$cycle==all_cycles[i], c("unit", "nb_births")], period=period, graph=F)$nb_births
    # proportion of births out from the period around median birth date
    mat.med[2, i] <- sum(pattern$nb_births[pattern$cycle==all_cycles[i]], na.rm=T)-mat.med[1, i]
  }
  khi2 <- chisq.test(mat.med, simulate.p.value=T) # simulate.p.value=T in case less than 5 births per date
  pvalue <- khi2$p.value
  if (is.na(pvalue)==F) { # when all births occur in the period selected around median birth date
    if (pvalue<=0.05) {
      recurr_synchro <- F
    } else if (pvalue>0.05) {
      recurr_synchro <- T
    }
  } else {
    recurr_synchro <- NA
  }
  
  # PLOT
  if (graph==T) {
    par(mfrow=c(1, 1))
    plot(mat.med[1, ]~all_cycles, pch=16, col="black", xlab="Cycle", ylab="Number of births around the median birth date", 
         ylim=c(0, max(mat.med[1, ])+1)) 
  }
  
  # OUTPUT
  return(list(khi2_prop_births_median=recurr_synchro, p_value=pvalue, x_squared=khi2$statistic)) 
}

# KOLMOGOROV SMIRNOV MULTI
pheno.schaik <- function(pattern, graph=F) { 
  # DOCUMENTATION
  ## Description: are the birth distributions similar between the two cycles? Kolmogorov Smirnov test (H0: the distributions are identical).
  ## Arguments: 
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: are the births distributions similar, p value, statistic of the test.
  ## Details: the data should not contain any ex aequo.
  
  # INITIALISATIONS
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle) 
  if(length(all_cycles)>2) { warning(paste("More than two cycles detected, the function works with two cycles at once", 
                                           "(the function is applied on the two first cycles only).", sep=" "))}
  if (length(which(duplicated(pattern$nb_births)==T))>0) { warning("Ex aequos, test can be inaccurate") }
  sub_list <- vector("list", 2)
  
  # FUNCTION
  kolmo_test <- ks.test(pattern$nb_births[pattern$cycle %in% all_cycles[1]], pattern$nb_births[pattern$cycle %in% all_cycles[2]], 
                        alternative="two.sided")
  pvalue <- kolmo_test$p.value
  if (pvalue<0.05) {
    signif <- F
  } else {
    signif <- T
  }
  
  # PLOT
  if (graph==T) {
    plot(pattern$nb_births[pattern$cycle %in% all_cycles[1]]~pattern$unit[pattern$cycle %in% all_cycles[1]], type="l", xlab="unit", 
         ylab="nb births")
    lines(pattern$nb_births[pattern$cycle %in% all_cycles[2]]~pattern$unit[pattern$cycle %in% all_cycles[2]], col="grey")
    legend("topright", legend=c("cycle 1", "cycle 2"), col=c("black", "grey"), lty=1)
  }
  
  # OUTPUT
  return(list(kolmogorov_smirnov_multi=signif, p_value=pvalue, statistic_d=kolmo_test$statistic))
}

# DIFF SYNCHRONY LINEAR
pheno.paoli2 <- function(pattern, graph=F) {
  # DOCUMENTATION
  ## Description: is there a trend in the evolution of birth period duration across time? using a linear regression.
  ## Arguments:
  # pattern = a data frame with three columns (col1=first time unit, col2=second time unit, col3=number of births).
  ## Values: coefficient of determination of the relationship between birth period length and time, slope coefficient.
  ## Details: works with unimodal distribution. The function uses pheno.paoli(). The coefficient of determination allows to check if there is a trend or not, and if
  # the slope coefficient can be interpreted or not (if there is no trend, the slope coefficient is meaningless).
  
  # INITIALISATIONS
  colnames(pattern) <- c("cycle", "unit", "nb_births")
  all_cycles <- unique(pattern$cycle)
  
  # FUNCTION
  patternbis <- data.frame(cycle=all_cycles, synchrony=NA)
  for (i in 1:length(all_cycles)) {
    patternbis$synchrony[i] <- pheno.paoli(pattern[pattern$cycle==all_cycles[i], c("unit", "nb_births")])
  }
  mod <- lm(synchrony~cycle, data=patternbis)
  patternbis$pred <- predict(mod, type="response") 
  diff_synchrony_linear <- mod$coefficients[2]
  rsquared <- summary(mod)$r.squared
  
  # PLOT
  if (graph==T) {
    plot(patternbis$synchrony~patternbis$cycle, xlab="Cycle", ylab="Synchrony of births")
    lines(patternbis$pred~patternbis$cycle, col="red")
  }
  
  # OUTPUT
  return(list(diff_synchrony_linear=diff_synchrony_linear, rsquared_diff_mean_linear=rsquared))
}

####################################################################################################################################################################
####################################################################################################################################################################
