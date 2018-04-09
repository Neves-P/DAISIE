island_area<-function(t,Apars,shape){
  if(shape==0){return(Apars[2])}	
  if(shape==1){
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
  if(shape == 2){
    Tmax <- Apars[1]
    proptime<- t/Tmax
    b <- Apars[2] # intercept (peak area)
    m <- -(b / Tmax) # slope
    At <- m * t + b
    return(At)
  }
}

extcutoff<-max(1000,1000*(ana0+clado0+immig0))

#Function to describe relationship between area and extinction rate
getExtRate <- function(t,Apars,Epars,shape,extcutoff){
  # X <- log(Epars[1] / Epars[2]) / log(0.1)
  X <- Epars[1]/Epars[2]
  extrate <- Epars[1] / ((island_area(t, Apars, shape) / Apars[2])^X) # / inverts function shape
  extrate[which(extrate > extcutoff)] <- extcutoff
  return(extrate)
}
res_area <- c()
for(i in 1:length(vec)){res_area[i] <- island_area(vec[[i]], Apars = c(10, 100, 1, 90), shape = 1)}
res_ext <- c()
for(i in 1:length(vec)){res_ext[i] <- getExtRate(vec[[i]], Apars = c(10, 100, 1, 90), Epars, shape = 1, 1000)}
plot(res_area)
plot(res_ext)
