# Function to describe changes in area through time. Adapted from
# Valente et al 2014 ProcB
island_area<-function(t, time, Apars, shape = NULL){
  if(is.null(shape)){return(Apars[1])}	
  if(shape == "quadratic"){
    if(Apars[2] > time){
      stop("Apars[2] > time: Peak position cannot be higher than total time.")
    }
    Tmax <- time # total time
    Amax <- Apars[1] # maximum area
    Topt<-Apars[2] # peak position
    peak<-Apars[3] # peakiness - we specify a value of 1 but this is flexible.
    proptime<- t/Tmax	
    f <-Topt/Tmax / (1 - Topt/Tmax)
    a <-f*peak/(1+f)
    b <-peak/(1+f) 
    At<- Amax * proptime^a * (1 - proptime)^ b / ((a / (a + b))^a * (b / (a + b))^b)
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

##### DOUBLE CHECK APARS INDICES

# Function to describe changes extinction rate through time. Adapted from
# Valente et al 2014 ProcB
getExtRate <- function(t, time, Apars, mu, shape, extcutoff, mu_version){
  Epars <- c()
  Epars[1] <- max(0.01, mu - mu/2) # Epars contains mu_min and mu_max for quadratic and logistic
  Epars[2] <- mu + mu/2
  
  if(mu_version == "procb"){
    X <- log(Epars[1] / Epars[2]) / log(0.1)
    extrate <- Epars[1] / ((island_area(t, time, Apars, shape) / time)^X) 
  }else if(mu_version == "muMax_power"){
    X <- mu
    extrate <- -((island_area(t, time, Apars, shape) / time))^(X) + 1 + Epars[1] # Epars1 is basal extinction
  }else if(mu_version == "logistic"){
    extrate <- Epars[1] + 1 / (island_area(t, time, Apars, shape) + 10) # logistic function
  }else{
    stop("Please insert valid mu function version.")
  }
  # print(paste0("area/areaMax at t=", t, ": ", (island_area(t, Apars, shape) / Apars[2])))
  # print(paste0("mu = ",extrate))
  extrate[which(extrate > extcutoff)] <- extcutoff
  # print(island_area(t, Apars, shape))
  return(extrate)
}


