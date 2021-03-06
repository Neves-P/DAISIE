\name{DAISIE_sim}
\alias{DAISIE_sim}
\title{Simulate islands with given parameters.}
\description{ This function simulates islands with given cladogenesis, extinction,
Kprime, immigration and anagenesis parameters. If a single parameter 
set is provided (5 parameters) it simulates islands where all species 
have the same macro-evolutionary process. If two paramater sets (10 parameters) are
provided, it simulates islands where two different macro-evolutionary processes 
operate, one applying to type 1 species and other to type 2 species.}
\usage{
DAISIE_sim(
   time,
   M,
   pars,
   replicates,
   divdepmodel = 'CS',
   prop_type2_pool = NA,
   replicates_apply_type2 = TRUE,
   sample_freq = 25,
   plot_sims = TRUE,
   island_ontogeny = NULL,
   Apars = NULL
   ...)
}

\arguments{
  \item{time}{Length of the simulation in time units. For example, 
if an island is know to be 4 million years old, setting time = 4 will 
simulate entire life span of the island; setting time = 2 will stop 
the simulation at the mid-life of the island.}
\item{M}{The size of the mainland pool, i.e the number of species that can
 potentially colonize the island }
\item{pars}{ Contains the model parameters: \cr \cr
   \code{pars[1]} corresponds to lambda^c (cladogenesis rate) \cr
   \code{pars[2]} corresponds to mu (extinction rate) \cr
   \code{pars[3]} corresponds to K (clade-level carrying capacity). Set K=Inf 
for non-diversity dependence.\cr
   \code{pars[4]} corresponds to gamma (immigration rate) \cr
   \code{pars[5]} corresponds to lambda^a (anagenesis rate) \cr   
   \code{pars[6]} corresponds to lambda^c (cladogenesis rate) for type 2 species \cr
   \code{pars[7]} corresponds to mu (extinction rate) for type 2 species\cr
   \code{pars[8]} corresponds to K (clade-level carrying capacity) for type 2 species.
 Set K=Inf for non-diversity dependence.\cr
   \code{pars[9]} corresponds to gamma (immigration rate) for type 2 species\cr
   \code{pars[10]} corresponds to lambda^a (anagenesis rate) for type 2 species\cr   
   The elements 6:10 are optional and are required only when type 2 species are included.
   }  
\item{replicates}{Number of island replicates to be simulated.}
\item{divdepmodel}{Option divdepmodel='CS' runs model with clade-specific carrying capacity, where diversity-dependence operates only within single clades, i.e. only among species originating from the same mainland colonist. Option divdepmodel= 'IW' runs model with island-wide carrying capacity, where diversity-dependence operates within and among clades.}
\item{prop_type2_pool}{Fraction of mainland species that belongs to the 
second subset of species (type 2). Applies only when two types of species are simulated (length(pars)=10).}
 \item{replicates_apply_type2}{ Applies only when two types of species are being simulated. 
Default replicates_apply_type2=TRUE
 runs simulations until the number of islands where a type 2 species has colonised is equal to the specified number of replicates. This is recommended if prop_type2_pool is small or if the rate of immigration of type two species (pars[9]) is low, meaning that more 
replicates are needed to achieve an adequate sample size of islands with type 2 species. Setting replicates_apply_type2=FALSE, simulates islands up 
to specified number of replicates regardless of whether type 2 species have colonised or not.}
\item{sample_freq}{Specifies the number of units time should be divided by for 
plotting purposes. Larger values will lead to 
plots with higher definition, but will also run slower.}
\item{plot_sims}{Default=TRUE plots species-through-time (STT) plots. It detects how many types
 of species are present. If only one type of species is present, STT is plotted for all species. If two types are present, 
three plots are produced: STT for all, STT for type 1 and STT for type 2. }
\item{island_ontogeny}{The type of possible island ontogeny functions. Setting quadratic will model a function of time as a negative quadratic function which can be shaped by Apars. Setting linear will model a linear function with negative slope.}
\item{Apars}{Vector of parameters controlling island ontogeny area function. }
  \code{Apars[1]}{Maximum area of the island.}
  \code{Apars[2]}{Peak position of quadratic function. If in middle of 
  time, then a negative quadratic function is produced. If lower than
  total time then a hump shaped function resembling volcanic area is 
  produced}
  \code{Apars[3]}{Peakiness. How bell shaped the function is.}
\item{mu_version}{Selects function to relate extinction rate to area ontogeny. procb uses function implemented in 2014 procb paper. logistic implements function similar to logistic function. muMax_version uses a simple power law of mu rate + mu0 * 0.5}
\item{...}{ Any arguments to pass on to plotting functions.}
}
\details{Returns R list object that contains the simulated islands.}
\value{
Each simulated dataset is an element of the list, which can be called using [[x]]. For example if the object is called island_replicates, the 1st replicate can be called using island_replicates[[1]]
Each of the island replicates is a list in itself. The first (e.g. island_replicates[[x]][[1]]) element of that list has the following components: \cr
  \code{$island_age} - the island age \cr
  Then, depending on whether a distinction between types is made, we have:\cr
  \code{$not_present} - the number of mainland lineages that are not present on the island \cr  
  or:\cr
  \code{$not_present_type1} - the number of mainland lineages of type 1 that are not present on the island \cr  
  \code{$not_present_type2} - the number of mainland lineages of type 2 that are not present on the island \cr 
  \code{$stt_all} - STT table for all species on the island (nI - number of non-endemic species; nA - number of anagenetic species, nC - number of cladogenetic species, present - number of independent colonisations present )\cr
  \code{$stt_stt_type1} - STT table for type 1 species on the island - only if 2 types of species were simulated (nI - number of non-endemic species; nA - number of anagenetic species, nC - number of cladogenetic species, present - number of independent colonisations present )\cr   
  \code{$stt_stt_type2} - STT table for type 2 species on the island - only if 2 types of species were simulated (nI - number of non-endemic species; nA - number of anagenetic species, nC - number of cladogenetic species, present - number of independent colonisations present )\cr   
  \code{$brts_table} - Only for simulations under 'IW'. Table containing information on order of events in the data, for use in maximum likelihood optimization.)\cr   

  The subsequent elements of the list each contain information on a single colonist lineage on the island and has 4 components:\cr
  \code{$branching_times} - island age and stem age of the population/species in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species these should be island age and branching times of the radiation including the stem age of the radiation.\cr                          
  \code{$stac} - the status of the colonist \cr
  * Non_endemic_MaxAge: 1 \cr
  * Endemic: 2 \cr
  * Endemic&Non_Endemic: 3 \cr
  * Non_endemic: 4 \cr
  \code{$missing_species} - number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
  \code{$type_1or2} - whether the colonist belongs to type 1 or type 2 \cr
 
}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\author{
Luis Valente and Albert Phillimore
}
\seealso{
\code{\link{DAISIE_format_CS}}
\code{\link{DAISIE_plot_sims}}
}
\examples{ 
cat("
## Simulate 40 islands for 4 million years, where all species have equal 
## rates, and plot the species-through-time plot. Pool size 1000.

pars_equal = c(2.550687345,2.683454548,Inf,0.00933207,1.010073119)
island_replicates_equal = DAISIE_sim(
   time = 4,
   M = 1000,
   pars = pars_equal,
   replicates = 40
   )

## Simulate 15 islands for 4 million years with two types of species (type1 
## and type 2), and plot the species-through-time plot. Pool size 1000. Fraction
## of type 2 species in source pool is 0.163. Function will simulate until number of islands
## where type 2 species has colonised is equal to number specified in replicates.

pars_type1 = c(0.195442017,0.087959583,Inf,0.002247364,0.873605049)
pars_type2 = c(3755.202241,8.909285094,14.99999923,0.002247364,0.873605049)
island_replicates_2types = DAISIE_sim(
   time = 4,
   M = 1000,
   pars = c(pars_type1,pars_type2),
   replicates = 15,
   prop_type2_pool = 0.163
   )
")
}
\keyword{ models }