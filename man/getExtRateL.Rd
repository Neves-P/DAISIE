\name{getExtRate}
\alias{getExtRate}
\title{Function to obtain extention rate from current area.
}
\description{This function evaluates island area at current time and calculates 
the extinction rate according to one of three possible functions of area.
}
\usage{
getExtRate(
  t,
  time,
  Apars,
  mu,
  shape,
  extcutoff,
  mu_version)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{Current time}
  \item{time}{Total time of simulation}
\item{Apars}{Vector of parameters controlling island ontogeny area function. }
  \code{Apars[1]}{Maximum area of the island.}
  \code{Apars[2]}{Peak position of quadratic function. If in middle of 
  time, then a negative quadratic function is produced. If lower than
  total time then a hump shaped function resembling volcanic area is 
  produced}
  \code{Apars[3]}{Peakiness. How bell shaped the function is.}
\item{mu}{Basal extinction rate (mu0)}
\item{shape}{Type of function. Can be quadratic, linear or constant (if NULL))}
\item{extcutoff}{Maximum possible value of mu to prevent simulation taking to long}
\item{mu_version}{Selects function to relate extinction rate to area ontogeny. procb uses function implemented in 2014 procb paper. logistic implements function similar to logistic function. muMax_version uses a simple power law of mu rate + mu0 * 0.5}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric representing current area of island at time t.
}
\references{Valente, L.M., R.S. Etienne and A.B. Phillimore (2014). The Effects of Island Ontogeny on Species Diversity and Phylogeny. Proceedings of the Royal Society B: Biological Sciences 281 (1784): 20133227–20133227. 
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.
%% ~put references to the literature/web site here ~
}
\author{Valente, Etienne, Santos Neves, Phillimore 
}
\note{Code was adapted from code available in supplementary information of Valente et al. 2014
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DAISIE_sim}}
\code{\link{island_area}}
}
\examples{

Apars <-  c(100, 5, 1)

res <- c()
for(i in 1:100){
  res[i] <- getExtRate(t = i, time = 100, Apars, mu = 0.2, shape = "quadratic",
             extcutoff = 1100, mu_version = "muMax_power")
}

plot(res)

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
