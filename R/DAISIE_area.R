# Function to describe changes in area through time. Adapted from
# Valente et al 2014 ProcB

island_area <- function(t, totaltime, Apars, shape){

  Tmax <- totaltime # total time
  Amax <- Apars[1] # maximum area
  Topt <- Apars[2] # peak position in %
  peak <- Apars[3] # peakiness - we specify a value of 1 but this is flexible.
  proptime<- t/Tmax	
  # Constant
  if (is.null(shape)){
    return(Apars[1])
  }	
  if(shape == "quadratic"){

    f <- Topt / (1 - Topt)
    a <- f * peak/ ( 1 + f)
    b <- peak / (1 + f) 
    At <- Amax * proptime^a * (1 - proptime)^ b/ ((a / (a + b))^a * (b / (a + b))^b)
    return(At)}
  
  #Linear decline
  if(shape == "linear"){
    b <- Amax # intercept (peak area)
    m <- -(b / Topt) # slope
    At <- m * t + b
    return(At)
  }
}


# Function to describe changes in extinction rate through time. From
# Valente et al 2014 ProcB
getExtRate <- function(t, totaltime, Apars, Epars, shape, extcutoff){
  X <- log(Epars[1] / Epars[2]) / log(0.1)
  extrate <- Epars[1]/((island_area(t, totaltime, Apars, shape) / Apars[1])^X)
  extrate[which(extrate > extcutoff)] <- extcutoff
  extrate[which(extrate > extcutoff)] <- extcutoff
  return(extrate)
}
