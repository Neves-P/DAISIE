
library(DAISIEontogeny)
pars <- c(0.5, 0.2, 40, 0.1, 0.2)
Apars <- c(1000, 5, 1)

parabola_island <- DAISIE_sim(10, 1000, pars, 1, island_ontogeny = "quadratic", Apars = Apars, 
           mu_version = "procb")

new_DAISIE_no_ontogeny <- DAISIE_sim(10, 1000, pars, 1, island_ontogeny = NULL, Apars = NULL, 
                              mu_version = NULL)

old_DAISIE_no_ontogeny <- DAISIE::DAISIE_sim(10, 1000, pars, 1)

parabola_island[[1]][[1]]$stt_all

new_DAISIE_no_ontogeny[[1]][[1]]$stt_all

old_DAISIE_no_ontogeny[[1]][[1]]$stt_all
