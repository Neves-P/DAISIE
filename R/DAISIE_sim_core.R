DAISIE_sim_core <- function(totaltime,
                            mainland_n,
                            pars,
                            Apars = NULL,
                            Epars = NULL,
                            island_ontogeny = NULL)
{
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  
  extcutoff <- max(1000, 1000 * (laa + lac + gam))
  ext_multiplier <- 0.5
  
  if(pars[4] == 0) 
  {
    stop('Rate of colonisation is zero. Island cannot be colonised.')
  }  
  
  timeval <- 0
  
  mainland_spec <- seq(1,mainland_n,1)
  maxspecID <- mainland_n
  
  island_spec = c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time","nI","nA","nC")
  stt_table[1,] <- c(totaltime,0,0,0)
  
  #### Algorithm with no area changes ####
  if (is.null(island_ontogeny)) {
  while(timeval < totaltime)
  {  	
  	ext_rate <- mu * length(island_spec[,1])
  	ana_rate <- laa * length(which(island_spec[,4] == "I"))
  	clado_rate <- max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
  	immig_rate <- max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
  
  	totalrate <- ext_rate + clado_rate + ana_rate + immig_rate
  	dt <- rexp(1,totalrate)
  	
  	timeval <- timeval + dt
  	
  	possible_event <- sample(1:4,1,replace=FALSE,c(immig_rate,ext_rate,ana_rate,clado_rate))
  
    ##############
    if(timeval <= totaltime)
  	{ 
  	  new_state <- DAISIE_sim_update_state(timeval, possible_event,maxspecID,mainland_spec,island_spec)
  	  island_spec <- new_state$island_spec
  	  maxspecID <- new_state$maxspecID
  	}
    stt_table <- rbind(stt_table,
      c(totaltime - timeval,
        length(which(island_spec[,4] == "I")),
        length(which(island_spec[,4] == "A")),
        length(which(island_spec[,4] == "C"))
        )
      )
  }
  
  stt_table[nrow(stt_table),1] <- 0
  
  ######## Algorithm if area changes ########
  } else if (!is.null(island_ontogeny) && (xor(island_ontogeny != "quadratic", 
                                               island_ontogeny != "linear") ||
                                           xor(island_ontogeny != "linear",
                                               island_ontogeny !=  "constant"))) {
    

    # Determine totaltime where A is max and thor (horizon totaltime to change rates)
    time_area_max <- Apars[2] * totaltime # totaltime where area is max
    thor <- time_area_max
    
    # Determine rates
    ext_rate <- getExtRate(t = timeval, totaltime = totaltime, Apars = Apars, 
                           Epars = Epars, shape = island_ontogeny, 
                           extcutoff = extcutoff) * length(island_spec[, 1])
    ext_rate_max <- getExtRate(t = timeval, totaltime = totaltime,
                               Apars = Apars, 
                               Epars = Epars, shape = island_ontogeny, 
                               extcutoff = extcutoff) * length(island_spec[, 1])
    ana_rate = laa * length(which(island_spec[,4] == "I"))
    clado_rate = max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
    immig_rate = max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
    
    # Pick timeval
    
    totalrate <- ext_rate_max + ana_rate + clado_rate + immig_rate
    dt <- rexp(1, totalrate)
    timeval <- timeval + dt
    
    # Checks if timeval is larger than thor from the start and jumps simulation
    # to thor if that's the case from the start 
    if (timeval > thor) {
      timeval <- thor
    }
    
    while(timeval < totaltime) {
      if (timeval < thor) {
        
        # Determine event
        event <- DDD::sample2(1:4, 1, prob = c(immig_rate, ext_rate_max,
                                               ana_rate, clado_rate),
                              replace = FALSE)
        
        if (event == 2) {
          if (ext_rate_max == 0) {
            possible_event <- 0
          } else {
            possible_event <- 2 * (runif(1) < ext_rate / ext_rate_max) # event is indeed extinction
          }
        } else {
          possible_event <- event
        }
        
        # Run event
        new_state <- DAISIE_sim_update_state(timeval, possible_event, maxspecID, mainland_spec, island_spec)
        island_spec <- new_state$island_spec
        maxspecID <- new_state$maxspecID
        
        # Update state
        stt_table <- rbind(stt_table, c(totaltime - timeval,
                                        length(which(island_spec[, 4] == "I")),
                                        length(which(island_spec[, 4] == "A")),
                                        length(which(island_spec[, 4] == "C"))))
        
        # Recalculate rates
        ext_rate_max <- getExtRate(t = timeval, totaltime = totaltime,
                                   Apars = Apars, 
                                   Epars = Epars, shape = island_ontogeny,
                                   extcutoff = extcutoff) * length(island_spec[, 1])
        ext_rate <- getExtRate(t = timeval, totaltime = totaltime, Apars = Apars,
                               Epars = Epars, shape = island_ontogeny,
                               extcutoff = extcutoff) * length(island_spec[, 1])
        ana_rate = laa * length(which(island_spec[,4] == "I"))
        clado_rate = max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
        immig_rate = max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
        
        totalrate <- ext_rate_max + ana_rate + immig_rate + clado_rate
        dt <- rexp(1, totalrate)
        timeval <- timeval + dt
        
        ##### After thor is reached ####
        
      } else {
        
        
        # Update timeval
        
        timeval <- thor
        
        # Recalculate thor
        thor <- timeval + ext_multiplier * (totaltime - timeval)
        
        # Recalculate rates
        ext_rate <- getExtRate(t = timeval, totaltime = totaltime, Apars = Apars, 
                               Epars = Epars, shape = island_ontogeny, 
                               extcutoff = extcutoff) * length(island_spec[,1])
        
        ext_rate_max <- getExtRate(t = thor, totaltime = totaltime,
                                   Apars = Apars, 
                                   Epars = Epars, shape = island_ontogeny, 
                                   extcutoff = extcutoff) * length(island_spec[,1])
        
        ana_rate = laa * length(which(island_spec[,4] == "I"))
        clado_rate = max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
        immig_rate = max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
        
        # Determine timeval
        totalrate <- ext_rate_max + ana_rate + immig_rate + clado_rate
        dt <- rexp(1, totalrate)
        timeval <- timeval + dt
        
        
        # Update stt table
        stt_table = rbind(stt_table,
                          c(totaltime - timeval,
                            length(which(island_spec[, 4] == "I")),
                            length(which(island_spec[, 4] == "A")),
                            length(which(island_spec[, 4] == "C"))))
      }
    }
    
    stt_table[nrow(stt_table),1] <- 0
    
  } else {
    
    stop("Please select valid island area function.")
  } 
  
  
  ############# 
  ### if there are no species on the island branching_times = island_age, stac = 0, missing_species = 0 
  if(length(island_spec[,1]) == 0)
  {
    island <- list(stt_table = stt_table, branching_times = totaltime, stac = 0, missing_species = 0)
  } else
  {
    cnames <- c("Species","Mainland Ancestor","Colonisation totaltime (BP)",
      "Species type","branch_code","branching totaltime (BP)","Anagenetic_origin")
    colnames(island_spec) <- cnames

    ### set ages as counting backwards from present
    island_spec[,"branching totaltime (BP)"] <- totaltime - as.numeric(island_spec[,"branching totaltime (BP)"])
    island_spec[,"Colonisation totaltime (BP)"] <- totaltime - as.numeric(island_spec[,"Colonisation totaltime (BP)"])
      
    if(mainland_n == 1)
    {
      island <- DAISIE_ONEcolonist(totaltime,island_spec,stt_table)
    }
      
    if(mainland_n > 1)
    {  
      ### number of colonists present
      colonists_present <- sort(as.numeric(unique(island_spec[,'Mainland Ancestor'])))
      number_colonists_present <- length(colonists_present) 
        
      island_clades_info <- list()  
      for(i in 1:number_colonists_present)
      {
        subset_island <- island_spec[which(island_spec[,'Mainland Ancestor']==colonists_present[i]),] 
        if(class(subset_island) != 'matrix')
        {
          subset_island <- rbind(subset_island[1:7])
          colnames(subset_island) <- cnames
        }
        island_clades_info[[i]] <- DAISIE_ONEcolonist(totaltime,island_spec=subset_island,stt_table=NULL)
        island_clades_info[[i]]$stt_table <- NULL
      }
      island <- list(stt_table = stt_table, taxon_list = island_clades_info)
    }
  }
  return(island) 
}

DAISIE_sim_update_state <- function(timeval, possible_event,maxspecID,mainland_spec,island_spec)
{  
  ##########################################
  #IMMIGRATION
  if(possible_event == 1)
  {  	
    colonist = DDD::sample2(mainland_spec,1)
    
    if(length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
    } else
    {
      isitthere = c()
    }
    
    if(length(isitthere) == 0)
    {
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA))
    }
    
    if(length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA)
    }
  }
  
  ##########################################
  #EXTINCTION
  if(possible_event == 2)
  { 	
    extinct = DDD::sample2(1:length(island_spec[,1]),1)
    #this chooses the row of species data to remove
    
    typeofspecies = island_spec[extinct,4]
    
    if(typeofspecies == "I")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove immigrant
    
    if(typeofspecies == "A")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove anagenetic
    
    if(typeofspecies == "C")
    {
      #remove cladogenetic
      
      #first find species with same ancestor AND arrival totaltime
      sisters = intersect(which(island_spec[,2] == island_spec[extinct,2]),which(island_spec[,3] == island_spec[extinct,3]))
      survivors = sisters[which(sisters != extinct)]
      
      if(length(sisters) == 2)
      {
        #survivors status becomes anagenetic	
        island_spec[survivors,4] = "A"
        island_spec[survivors,c(5,6)] = c(NA,NA)
        island_spec[survivors,7] = "Clado_extinct"
        island_spec = island_spec[-extinct,]
      }
      
      if(length(sisters) >= 3)
      {		
        numberofsplits = nchar(island_spec[extinct,5])
        
        mostrecentspl = substring(island_spec[extinct,5],numberofsplits)
        
        if(mostrecentspl=="B")
        { 
          sistermostrecentspl = "A"
        }
        if(mostrecentspl=="A")
        {
          sistermostrecentspl = "B"
        }
        
        motiftofind = paste(substring(island_spec[extinct,5],1,numberofsplits-1),sistermostrecentspl,sep = "")
        
        possiblesister = survivors[which(substring(island_spec[survivors,5],1,numberofsplits) == motiftofind)]
        
        #different rules depending on whether a B or A is removed. B going extinct is simpler because it only 
        #carries a record of the most recent speciation			
        if(mostrecentspl == "A")
        {								
          #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
          tochange = possiblesister[which(island_spec[possiblesister,6] == max(as.numeric(island_spec[possiblesister,6])))]
          island_spec[tochange,6] = island_spec[extinct,6]	
        }
        
        #remove the offending A/B from these species
        island_spec[possiblesister,5] = paste(substring(island_spec[possiblesister,5],1,numberofsplits - 1),
                                              substring(island_spec[possiblesister,5],numberofsplits + 1,
                                                        nchar(island_spec[possiblesister,5])),sep = "")	
        island_spec = island_spec[-extinct,]
      }
    }
    island_spec = rbind(island_spec)	
  }
  
  ##########################################
  #ANAGENESIS
  if(possible_event == 3)
  {    
    immi_specs = which(island_spec[,4] == "I")
    
    #we only allow immigrants to undergo anagenesis
    if(length(immi_specs) == 1)
    {
      anagenesis = immi_specs
    }
    if(length(immi_specs) > 1)
    {
      anagenesis = DDD::sample2(immi_specs,1)
    }
    
    maxspecID = maxspecID + 1
    island_spec[anagenesis,4] = "A"
    island_spec[anagenesis,1] = maxspecID
    island_spec[anagenesis,7] = "Immig_parent"
  }
  
  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive 
  if(possible_event == 4)
  { 		
    tosplit = DDD::sample2(1:length(island_spec[,1]),1)
    
    #if the species that speciates is cladogenetic
    if(island_spec[tosplit,4] == "C")
    {
      #for daughter A
      
      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      oldstatus = island_spec[tosplit,5]
      island_spec[tosplit,5] = paste(oldstatus,"A",sep = "")
      #island_spec[tosplit,6] = timeval
      island_spec[tosplit,7] = NA
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],
                                        "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      
      maxspecID = maxspecID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      
      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      island_spec[tosplit,5] = "A"
      island_spec[tosplit,6] = island_spec[tosplit,3]
      island_spec[tosplit,7] = NA
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA))
      
      maxspecID = maxspecID + 2
    }
  }
  return(list(island_spec = island_spec,maxspecID = maxspecID))
}

DAISIE_ONEcolonist <- function(totaltime,island_spec,stt_table)
{
  
  ### number of independent colonisations
  uniquecolonisation <- as.numeric(unique(island_spec[,"Colonisation totaltime (BP)"]))
  number_colonisations <- length(uniquecolonisation) 
  
  ### if there is only one independent colonisation - anagenetic and cladogenetic
  #species are classed as stac=2; immigrant classed as stac=4: 
  if(number_colonisations == 1)
  {
    if(island_spec[1,"Species type"] == "I")
    {
      descendants <- list(stt_table = stt_table, 
                          branching_times = c(totaltime,as.numeric(island_spec[1,"Colonisation totaltime (BP)"])),
                          stac = 4,
                          missing_species = 0)
    }
    if(island_spec[1,"Species type"] == "A")
    {
      descendants <- list(stt_table = stt_table,
                         branching_times = c(totaltime,as.numeric(island_spec[1,"Colonisation totaltime (BP)"])),
                         stac = 2,
                         missing_species = 0)
    } 
    if(island_spec[1,"Species type"] == "C")
    {
      descendants <- list(stt_table = stt_table,
                         branching_times = c(totaltime,rev(sort(as.numeric(island_spec[,"branching totaltime (BP)"])))),
                         stac = 2,
                         missing_species = 0)
    }
  }
  
  ### if there are two or more independent colonisations, all species are classed as stac=3 and put within same list item: 
  else if(number_colonisations > 1)
  {
    descendants <- list(stt_table = stt_table,
                       branching_times = NA,stac = 2,missing_species = 0,
                       other_clades_same_ancestor = list())
 
    ### create table with information on other clades with same ancestor, but this information is not used in DAISIE_ML
    oldest <- which(as.numeric(island_spec[,"Colonisation totaltime (BP)"]) == max(as.numeric(island_spec[,"Colonisation totaltime (BP)"])))
    
    oldest_table <- island_spec[oldest,]
    if(class(oldest_table) == 'character')
    { 
      oldest_table <- t(as.matrix(oldest_table))
    }
    if(oldest_table[1,'Species type'] == 'A')
    {
       descendants$branching_times <- c(totaltime, as.numeric(oldest_table[1,"Colonisation totaltime (BP)"]))
    } else if(oldest_table[1,'Species type'] == 'C')
    {
      descendants$branching_times <- c(totaltime, rev(sort(as.numeric(oldest_table[,'branching totaltime (BP)']))))
    }
    
    youngest_table = island_spec[-oldest,]
    if(class(youngest_table) == 'character')
    {
      youngest_table <- t(as.matrix(youngest_table))
    }
    
    uniquecol <- as.numeric(unique(youngest_table[,"Colonisation totaltime (BP)"]))
    
    descendants$missing_species <- length(which(youngest_table[,"Species type"]!='I'))
    
    for(colonisation in 1:length(uniquecol))
    {
      descendants$other_clades_same_ancestor[[colonisation]] <- list(brts_miss = NA,species_type = NA)	
      
      samecolonisation <- which(as.numeric(youngest_table[,"Colonisation totaltime (BP)"]) == uniquecol[colonisation])
      
      if(youngest_table[samecolonisation[1],"Species type"] == "I")
      {
        descendants$stac <- 3
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss <- as.numeric(youngest_table[samecolonisation,"Colonisation totaltime (BP)"])
        descendants$other_clades_same_ancestor[[colonisation]]$species_type <- "I"
      } else if(youngest_table[samecolonisation[1],"Species type"] == "A")
      {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss <- as.numeric(youngest_table[samecolonisation,"Colonisation totaltime (BP)"])
        descendants$other_clades_same_ancestor[[colonisation]]$species_type <- "A"
      } else if (youngest_table[samecolonisation[1],"Species type"] == "C")
      {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss <- rev(sort(as.numeric(youngest_table[samecolonisation,"branching totaltime (BP)"])))
        descendants$other_clades_same_ancestor[[colonisation]]$species_type <- "C"
      }
    }
  }
  return(descendants)  
}