####################################################################################################################################################################
#### BEST TIMING - SYNCHRONY - RHYTHMICITY - REGULARITY METRIC IDENTIFIED VIA THE EIGHT CRITERIA (FROM BIRTH PHENOLOGY IN LARGE HERBIVORES LITERATURE) #############
####################################################################################################################################################################

# DOCUMENTATION ####################################################################################################################################################

### DESCRIPTION
# Code of the four metrics identified as the "best" metrics to evaluate timing, synchrony, rhythmicity and regularity of the phenology of births in the large
# herbivores literature.
# One can find a generic function at the end of the code, which returns one output for each of the four metrics for a given phenology of births at the same time.
# see details in the "phenology_metrics" script.

### ARGUMENTS - VALUES - DETAILS
# For each function, a "graph" argument allows to choose if a plot should be returned or not.
# See details in each function.

# FUNCTIONS ########################################################################################################################################################

# MEAN VECTOR LENGTH (short name: meanvl), MEAN VECTOR ORIENTATION (short name: meanvo)
pheno.meanvl.meanvo.best <- function(pattern, graph=F) {
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
  return(list(mean_vector_orientation=mean_vect_orientation, #### meanvo ####
              mean_birth_date=mean_birth_date, 
              mean_vect_length=mean_vect_length))            #### meanvl ####
}

# MOOD (short name: mood)
pheno.mood.best <- function(pattern, graph=F) { 
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
    pattern_cum1 <- cumsum(pattern$nb_births[pattern$cycle==all_cycles[1]])
    pattern_cum2 <- cumsum(pattern$nb_births[pattern$cycle==all_cycles[2]])
    med_birth1 <- which(pattern_cum1>=sum(pattern$nb_births[pattern$cycle==all_cycles[1]])/2)[1]
    med_birth2 <- which(pattern_cum2>=sum(pattern$nb_births[pattern$cycle==all_cycles[2]])/2)[1]
    par(mfrow=c(2, 1))
    plot(pattern$nb_births[pattern$cycle==all_cycles[1]]~pattern$unit[pattern$cycle==all_cycles[1]], xlab="Time", ylab="Nb births", type="l")
    abline(v=med_birth1, col="blue")
    text(med_birth1, max(pattern$nb_births[pattern$cycle==all_cycles[1]])-1, "Median birth date", col="blue")
    plot(pattern$nb_births[pattern$cycle==all_cycles[2]]~pattern$unit[pattern$cycle==all_cycles[2]], xlab="Time", ylab="Nb births", type="l")
    abline(v=med_birth2, col="blue")    
    text(med_birth2, max(pattern$nb_births[pattern$cycle==all_cycles[2]])-1, "Median birth date", col="blue")
    par(mfrow=c(1, 1))
  }
  
  return(list(mood=signif, 
              p_value=pvalue, 
              statistic_z=mood@statistic@teststatistic)) #### mood ####
}

# KOLMOGOROV SMIRNOV MULTI-YEAR (short name: kolmomult)
pheno.kolmomult.best <- function(pattern, graph=F) { 
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
  return(list(kolmogorov_smirnov_multi=signif, 
              p_value=pvalue, 
              statistic_d=kolmo_test$statistic)) #### kolmomult ####
}

################################################################

# Generic function returning one output for each of the four "best" metrics, as done in the associated article
metrics.test.best4 <- function(pattern, graph=F) {
  
  ### Reshape the phenology of births to adjust to each kind of metrics (apply to one year, two years or more than two years)
  colnames(pattern) <- c("cycle", "unite", "nb_births") 
  pattern$cycle <- as.numeric(as.factor(pattern$cycle)) # this line initiate the counting of the years to 1
  # metrics which apply to one year : pattern.one
  pattern.one <- pattern[pattern$cycle==1, c("unite", "nb_births")] # draw first year
  # metrics which apply to two years : pattern.two
  tab_int1 <- pattern[pattern$cycle==1, ] # draw first two years
  tab_int2 <- pattern[pattern$cycle==2, ]
  tab_int <- rbind(tab_int1, tab_int2)
  tab_int$cycle <- rep(c(1, 2), each=(length(tab_int$unite)/2))
  pattern.two <- tab_int[ , c("cycle", "unite", "nb_births")]
  # metrics which apply to more than two years : pattern
  
  ### Data.frame of the outputs
  fin_sum <- data.frame(pattern=NA)
  
  ### Calculate the value of each metric for the phenology of births tested
  # One year
  meanvl_meanvo <- pheno.meanvl.meanvo.best(pattern.one, graph=graph)
  # Two years
  mood <- pheno.mood.best(pattern.two, graph=graph)
  kolmomult <- pheno.kolmomult.best(pattern.two, graph=graph) 

  ### Fill in the data.frame of the outputs
  fin_sum$pattern <- "BirthPhenology"
  
  fin_sum$kolmomult <- kolmomult$statistic_d
  fin_sum$meanvl <- meanvl_meanvo$mean_vect_length
  fin_sum$meanvo <- meanvl_meanvo$mean_vector_orientation
  fin_sum$mood <- mood$statistic_z
  
  ### Convert TRUE/FALSE to 1/0
  for (i in 1:ncol(fin_sum)) {
    fin_sum[ , i][fin_sum[ , i]==T] <- 1
  }
  
  ### Outputs
  return(fin_sum) # a data.frame with one row, and one column for each metric
}

####################################################################################################################################################################
####################################################################################################################################################################
