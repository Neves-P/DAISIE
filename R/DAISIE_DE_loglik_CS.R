DAISIE_DE_loglik <- function(pars1,
                             brts,
                             missnumspec,
                             stac,
                             methode,
                             reltolint = 1e-15,
                             abstolint = 1e-15,
                             sampling = 'rho') {
  if( stac == 0) {
    loglikelihood <- DAISIE_DE_logp0(island_age = brts[1],
                                     pars1 = pars1,
                                     reltolint = abstolint,
                                     abstolint = reltolint,
                                     methode = methode)
  } else if (stac == 1 || stac == 4 || stac == 8) {
    loglikelihood <- DAISIE_DE_logpNE(brts = brts,
                                      pars1 = pars1,
                                      stac = stac,
                                      methode = methode,
                                      reltolint = reltolint,
                                      abstolint = abstolint)

  } else if (stac == 2 && length(brts) == 2 || stac == 3 && length(brts) == 2 || stac == 5 && length(brts) == 2 || stac == 9) {
    if(sampling == 'rho' || missnumspec == 0)
      loglikelihood <- DAISIE_DE_logpES(brts = brts,
                                        missnumspec = missnumspec,
                                        stac = stac,
                                        pars1 = pars1,
                                        methode = methode,
                                        reltolint = reltolint,
                                        abstolint = abstolint) + fac
    else
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_function = DAISIE_DE_logpES,
                                   brts = brts,
                                   missnumspec = missnumspec,
                                   stac = stac,
                                   pars1 = pars1,
                                   methode = methode,
                                   reltolint = reltolint,
                                   abstolint = abstolint)
  } else if (stac == 2 && length(brts) > 2 || stac == 3 && length(brts) > 2 || stac == 6) {
    if(sampling == 'rho' || missnumspec == 0)
      loglikelihood <- DAISIE_DE_logpEC(brts = brts,
                                        missnumspec = missnumspec,
                                        stac = stac,
                                        pars1 = pars1,
                                        methode = methode,
                                        reltolint = reltolint,
                                        abstolint = abstolint) + fac
    else
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_function = DAISIE_DE_logpEC,
                                   brts = brts,
                                   missnumspec = missnumspec,
                                   stac = stac,
                                   pars1 = pars1,
                                   methode = methode,
                                   reltolint = reltolint,
                                   abstolint = abstolint)
  } else {
    stop("Unknown stac value: ", stac)
  }
  return(loglikelihood)
}

DAISIE_DE_loglik_CS <- function( pars1,
                                 pars2,
                                 datalist,
                                 methode = "odeint::runge_kutta_cash_karp54",
                                 abstolint = 1e-15,
                                 reltolint = 1e-15,
                                 equal_extinction = TRUE,
                                 sampling = 'n') {

  # Apply equal extinction condition AFTER initializing pars1
  if (equal_extinction) {
    pars1[3] <- pars1[2]
  }
  cond <- pars2[3]
  island_age <- datalist[[1]]$island_age
  parameter <- pars1
  if (length(parameter) == 5) {
    logp0 <- DAISIE_DE_logp0(island_age = island_age,
                             pars1 = pars1,
                             reltolint = reltolint, #was 1E-12
                             abstolint = abstolint, #was 1E-12
                             methode = methode)
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0))
    for (i in 2:length(datalist)) {
      datalist[[i]]$type1or2 <- 1
    }
  } else {
    numimm <- datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 <- length(which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2))
    numimm_type1 <- length(datalist) - 1 - numimm_type2
    if (!is.na(parameter[11])) {
      if (parameter[11] < numimm_type2 / numimm | parameter[11] > (1 - numimm_type1 / numimm)) {
        return(-Inf)
      }
      datalist[[1]]$not_present_type2 <- max(0, round(parameter[11] * numimm) - numimm_type2)
      datalist[[1]]$not_present_type1 <- numimm - (length(datalist) - 1) - datalist[[1]]$not_present_type2
    }
    logp0_type1 <- DAISIE_DE_logp0(island_age = island_age,
                                   pars1 = pars1,
                                   methode = methode,
                                   reltolint = reltolint,
                                   abstolint = abstolint)
    logp0_type2 <- DAISIE_DE_logp0(island_age = island_age,
                                   pars1 = pars1,
                                   methode = methode,
                                   reltolint = reltolint,
                                   abstolint = abstolint)
    loglik <- datalist[[1]]$not_present_type1 * logp0_type1 + datalist[[1]]$not_present_type2 * logp0_type2
    logcond <- (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 +
                                         (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
  }

  loglik <- loglik - logcond
  vec_loglikelihood <- c()

  for (i in 2:length(datalist)) {
    stac <- datalist[[i]]$stac
    brts <- datalist[[i]]$branching_times
    missnumspec <- datalist[[i]]$missing_species
    S <- length(brts) - 1
    fac <- S * (log(S) - log(S + missnumspec))
    loglikelihood <- DAISIE_DE_loglik(pars1 = pars1,
                                      brts = brts,
                                      missnumspec = missnumspec,
                                      stac = stac,
                                      methode = methode,
                                      reltolint = reltolint,
                                      abstolint = abstolint,
                                      sampling = sampling)

    vec_loglikelihood <- c(vec_loglikelihood, loglikelihood)

    print_parameters_and_loglik(
      pars = c(stac, parameter),
      loglik = loglikelihood,
      verbose = pars2[4],
      parnames = c("lambda^c", "mu", "K", "gamma", "lambda^a", "prob_init_pres"),
      type = 'clade_loglik'
    )
  }
  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)
}
