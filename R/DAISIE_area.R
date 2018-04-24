# Function to describe changes in area through time. Adapted from
# Valente et al 2014 ProcB
# 
# Convert Apars[2] into % of time to be easily applicable. Add function to convert
# time % to Topt

island_area <- function(t, time, Apars, shape){

  # Constant
  if (shape == "constant"){
    return(Apars[1])
  }	
  if(shape == "quadratic"){
    if(Apars[2] > time){
      stop("Apars[2] > time: Peak position cannot be higher than total time.")
    }
    aTmax <- time # total time
    Amax <- Apars[1] # maximum area
    Topt <- Apars[2] * time # peak position in %
    peak <- Apars[3] # peakiness - we specify a value of 1 but this is flexible.
    proptime<- t/Tmax	
    f <- Topt/Tmax / (1 - Topt/Tmax)
    a <- f*peak/(1+f)
    b <- peak/(1+f) 
    At <- Amax * proptime^a * (1 - proptime)^ b / ((a / (a + b))^a * (b / (a + b))^b)
    return(At)}
  
  #Linear decline
  if(shape == "linear"){
    Tmax <- time
    proptime<- t/Tmax
    b <- Apars[1] # intercept (peak area)
    m <- -(b / Tmax) # slope
    At <- m * t + b
    return(At)
  }
}


# Function to describe changes in extinction rate through time. From
# Valente et al 2014 ProcB
getExtRate <- function(t, time, Apars, Epars, shape, extcutoff){
  X <- log(Epars[1] / Epars[2]) / log(0.1)
  extrate <- Epars[1]/((island_area(t, time, Apars, shape) / Apars[1])^X)
  extrate[which(extrate > extcutoff)] <- extcutoff
  extrate[which(extrate > extcutoff)] <- extcutoff
  return(extrate)
}


