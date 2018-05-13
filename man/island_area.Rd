\name{island_area}
\alias{island_area}
\title{Function to calculate current island area
}
\description{This function evaluates calculates current area at time "t". Island can be constant, vary linearly or 
vary according to a beta function.
}
\usage{
island_area(
  t,
  time,
  Apars,
  shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{Current time}
  \item{time}{Total time of simulation}
\item{Apars}{Vector of parameters controlling island ontogeny area function.}
  \code{Apars[1]}{ Maximum area of the island.\cr}
  \code{Apars[2]}{ Peak position of quadratic function. If in middle of 
  time, then a negative quadratic function is produced. If lower than
  total time then a hump shaped function resembling volcanic area is 
  produced\cr}
  \code{Apars[3]}{ Peakiness. How bell shaped the function is.\cr}
\item{shape}{Type of function. Can be quadratic, linear or constant (if NULL))}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric representing current area of island at time t.
}
\references{Valente, L.M., R.S. Etienne and A.B. Phillimore (2014). The Effects of Island Ontogeny on Species Diversity and Phylogeny. Proceedings of the Royal Society B: Biological Sciences 281 (1784): 20133227-20133227. 
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
# Calculates area at time 2 for an island of total age 10
Apars <-  c(100, 5, 1)

current_area <- island_area(t = 2, time = 10, Apars, "quadratic")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
