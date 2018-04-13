# getExtRate <- function(t,Apars,Epars,shape,extcutoff){
#   # X <- log(Epars[1] / Epars[2]) / log(0.1)
#   # X <- (Epars[1]/Epars[2])
#   # X <- 0.5
#   X <- Epars[1]
#   extrate <- (island_area(t, Apars, shape) / Apars[2])^-X
#   print(paste0("area/areaMax at t=", t, ": ", (island_area(t, Apars, shape) / Apars[2])))
#   print(paste0("mu = ",extrate))
#   extrate[which(extrate > extcutoff)] <- extcutoff
#   return(extrate)
# }
mu_loop <- c()
mu_loop1 <- c()
for(i in 1:100){
  mu_loop[i] <- getExtRate(i, Apars, Epars, shape = "quadratic", extcutoff, "logistic")
  mu_loop1[i] <- getExtRate1(i, Apars, Epars, shape = "quadratic", extcutoff)
}

area_loop <- c()
for(i in 1:100){
  area_loop[i] <- island_area(i, Apars, shape = "quadratic")
}
plot(mu_loop)
plot(mu_loop1)
plot(-area_loop)
Apars
