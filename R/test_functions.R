test_ext <- function(Apars, mu, shape, extcutoff, mu_shape) {
  mu_loop <- c()
  for(i in 1:100){
    mu_loop[i] <- getExtRate(i, Apars, mu, shape, extcutoff, mu_version)
  }
}

test_area <- function(Apars, shape) {
  area_loop <- c()
  for(i in 1:100){
    area_loop[i] <- island_area(i, Apars, shape)
  }
}
