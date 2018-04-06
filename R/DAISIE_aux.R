island_area<-function(t,Apars,shape){
  if(shape==0){return(Apars[2])}	
  if(shape==1){
    Tmax <- Apars[1] # total time
    Amax <- Apars[2] # maximum area
    Topt<-Apars[3] #peak position
    peak<-Apars[4] #peakiness - we specify a value of 1 but this is flexible.
    proptime<- t/Tmax	
    f <-Topt/Tmax / (1 - Topt/Tmax)
    a <-f*peak/(1+f)
    b <-peak/(1+f) 
    At<- Amax * proptime^a * (1 - proptime)^b / ((a/(a+b))^a * (b/(a+b))^b)
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

#Function to describe relationship between area and extinction rate
ext<-function(t,Apars,Epars,shape,extcutoff){
  X<-log(Epars[1]/Epars[2])/log(0.1)
  extrate<-Epars[1]/((island_area(t,Apars,shape)/Apars[2])^X)
  extrate[which(extrate>extcutoff)]<-extcutoff
  return(extrate)}

