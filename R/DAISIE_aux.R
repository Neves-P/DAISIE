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
 # USAR EXPRESSÃO ANTIGA E DEPOIS MUDAR
extcutoff<-max(1000,1000*(ana0+clado0+immig0))

#Function to describe relationship between area and extinction rate
getExtRate <- function(t,Apars,Epars,shape,extcutoff){
  # X <- log(Epars[1] / Epars[2]) / log(0.1)
  X <- (Epars[1]/Epars[2])
  # X <- 0.5
  extrate <- (island_area(t, Apars, shape) / Apars[2])^-X
  print(paste0("area/areaMax at t=", t, ": ", (island_area(t, Apars, shape) / Apars[2])))
  print(paste0("mu = ",extrate))
  extrate[which(extrate > extcutoff)] <- extcutoff
  return(extrate)
}
res_area <- c()
for(i in 1:length(t_vec)){res_area[i] <- island_area(t_vec[[i]], Apars = c(100, 100, 25, 1), shape = 1)}
res_ext <- c()
for(i in 1:length(t_vec)){res_ext[i] <- getExtRate(t_vec[[i]], Apars = c(100, 100, 25, 1), Epars, shape = 1, 1000)}
plot(res_area, main = "Area. Apars = 10, 100, 5, 90")
plot(res_ext, main = "Ext 2018. Apars = 10, 100, 5, 90")


