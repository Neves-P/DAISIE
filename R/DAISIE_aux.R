island_area<-function(t,Apars,shape = NULL){
  if(is.null(shape)){return(Apars[2])}	
  if(shape == "quadratic"){
    if(Apars[3] > Apars[1]){
      stop("Apars[3] > Apars[1]: Peak position cannot be higher than total time.")
    }
    Tmax <- Apars[1] # total time
    Amax <- Apars[2] # maximum area
    Topt<-Apars[3] # peak position
    peak<-Apars[4] # peakiness - we specify a value of 1 but this is flexible.
    proptime<- t/Tmax	
    f <-Topt/Tmax / (1 - Topt/Tmax)
    a <-f*peak/(1+f)
    b <-peak/(1+f) 
    At<- Amax * proptime^a * (1 - proptime)^ b / ((a / (a + b))^a * (b / (a + b))^b)
    return(At)}
  
  #Linear decline
  if(shape == "linear"){
    Tmax <- Apars[1]
    proptime<- t/Tmax
    b <- Apars[2] # intercept (peak area)
    m <- -(b / Tmax) # slope
    At <- m * t + b
    return(At)
  }
}



#Function to describe relationship between area and extinction rate
getExtRate <- function(t,Apars,mu,shape,extcutoff, mu_version){
  Epars[1] <- min(0.01, mu - mu/2) # Epars contains mu_min and mu_max for quadratic and logistic
  Epars[2] <- mu + mu/2
  
  if(mu_version == "procb"){
    X <- log(Epars[1] / Epars[2]) / log(0.1)
    extrate <- Epars[1] / ((island_area(t, Apars, shape) / Apars[2])^X) 
  }else if(mu_version == "muMax_power"){
    X <- mu
    extrate <- -((island_area(t, Apars, shape) / Apars[2]))^(X) + 1 + Epars[1] # Epars1 is basal extinction
  }else if(mu_version == "logistic"){
    extrate <- Epars[1] + 1 / (island_area(t, Apars, shape) + 10) # logistic function
  }else{
    stop("Please insert valid mu function version.")
  }
  # print(paste0("area/areaMax at t=", t, ": ", (island_area(t, Apars, shape) / Apars[2])))
  # print(paste0("mu = ",extrate))
  extrate[which(extrate > extcutoff)] <- extcutoff
  # print(island_area(t, Apars, shape))
  return(extrate)
}

t_vec <- 0:99

# Test functions
res_area <- c()
for(i in 1:length(t_vec)){res_area[i] <- island_area(t_vec[[i]], Apars = c(100, 100, 25, 1), shape = "linear")}
res_ext <- c()
for(i in 1:length(t_vec)){res_ext[i] <- getExtRate(t_vec[[i]], Apars = c(100, 100, 25, 1), mu, shape = "linear", 1000, mu_version ="logistic")}
plot(res_area, main = paste0("Area. Apars = ", Apars))
plot(res_ext, main = "Ext 2018. Apars = 10, 100, 5, 90")

