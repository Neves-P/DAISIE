
mu_loop <- c()
for(i in 1:100){
  mu_loop[i] <- getExtRate(i, Apars, mu, shape = "quadratic", extcutoff, "logistic")
}

area_loop <- c()
for(i in 1:100){
  area_loop[i] <- island_area(i, Apars, shape = "quadratic")
}
plot(mu_loop)
plot(mu_loop1)
plot(-area_loop)
Apars
