\name{datagen}
\alias{datagen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation of random binary genotype 
}
\description{
Generation of random binary genotype satisfying both provided ploidy levels and allele frequency. It may present an interest to cheack the validity of an observed phenomenon, doing some bootstrap for signifiance etc… 
}
\usage{
datagen(frec, ploidy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frec}{
a vector giving the allele frequency of all allele, every single frequency have to be between 0 and 1}
  \item{ploidy}{
a vector giving the ploidy of every single sample, should be a non-null natural number 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Benjamin Leduc
}
\note{
Whenever you make something out of this function you have to precise that there are generated data. 
You shall not use this for getting fake experimental data (it constitues a violation of this package lisence ).
Whenever you use that function for bad purposes I will find you, get you and kick your ass! 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (frec, ploidy) 
{
    fakedata <- matrix(data = NA, nrow = length(ploidy), ncol = length(frec), 
        byrow = FALSE, dimnames = NULL)
    plolev <- unique(ploidy)
    lingen <- function(xx) {
        allfrec <- matrix(data = NA, nrow = length(ploidy), ncol = 1, 
            byrow = FALSE, dimnames = NULL)
        for (i in plolev) allfrec[ploidy == i, ] <- sample(c(0, 
            1), length(ploidy[ploidy == i]), replace = T, prob = c((1 - 
            xx)^i, 1 - ((1 - xx)^i)))
        allfrec
    }
    for (j in 1:length(frec)) fakedata[, j] <- lingen(frec[j])
    fakedata
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
