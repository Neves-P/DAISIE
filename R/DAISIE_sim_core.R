DAISIE_sim_core <- function(time, 
                            mainland_n, 
                            pars, 
                            Apars = NULL, 
                            Epars = NULL, 
                            island_ontongey = NULL)
{
  
  # This initial chunk is the same as in DAISIE
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  
  if (pars[4] == 0)
  {
    stop("Rate of colonization is zero. Island cannot be colonized")
  }
  
  timeval <- 0
  
  mainland_spec <- seq(1, mainland_n)
  maxspecID <- mainland_n
  
  island_spec <- c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) = c("Time", "nI", "nA", "nC")
  stt_table[1,] = c(time, 0, 0, 0)
 
  
  # while(timeval < time)
  # {  	
  #   ext_rate = mu * length(island_spec[,1])
  #   ana_rate = laa * length(which(island_spec[,4] == "I"))
  #   clado_rate = max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
  #   immig_rate = max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
  #   
  #   totalrate = ext_rate + clado_rate + ana_rate + immig_rate
  #   dt = rexp(1,totalrate)
  #   
  #   timeval  =  timeval  + dt
  #   
  #   possible_event = sample(1:4,1,replace=FALSE,c(immig_rate,ext_rate,ana_rate,clado_rate))
  #   
  
  # Gillespie algorithm 
  while (timeval < time){
    N <- length(island_spec[,1])
    
    # No ontogeny rates - (Style edited for consistency)
    if (is.null(island_ontogeny)){
    ext_rate <-mu * N
    ana_rate <- laa * N
    clado_rate <- max(N * (lac * (1 - N / K)), 0, na.rm = T)
    immig_rate <- max(c(mainland_n * gam * (1 - N / K), 0), na.rm = T)
    totalrate <- ext_rate + clado_rate + ana_rate + immig_rate
    dt <- rexp(1, totalrate)
    
    } else {
      
      ext_rate <- getExtRate(t = timeval, time = time, Apars = Apars, 
                             Epars = Epars, shape = island_ontogeny, 
                             extcutoff = extcutoff) * N
      ana_rate <- laa * length(which(island_spec[,4] == "I"))
      clado_rate <- max(c(N * (lac * (1 - N) / K), 0), na.rm = T)
      immig_rate = max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
      if (K == Inf){
        
      }
    }
  }
}