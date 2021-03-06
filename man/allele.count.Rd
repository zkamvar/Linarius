\name{allele.count}
\Rdversion{0.0.1}
\alias{allele.count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to count theorical allele amount with Dominant markers and mixed ploidy levels 
}
\description{
This function aims at counting allele presence and absence in a ballance population, this considering the ploidy level. 
}
\usage{
allele.count(xx, ploidy = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{
a binary datafram with allele in column and individuals in row. 
}
  \item{ploidy}{
a vector telling the ploidy level of each individual
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
%%  ~~further notes~~
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
function (xx, ploidy = 2) 
{
    plolev <- unique(ploidy)
    alco <- NULL
    names <- NULL
    for (i in plolev) alco <- cbind(alco, apply(as.matrix(xx[ploidy == 
        i, ]) == 0, 2, sum), apply(as.matrix(xx[ploidy == i, 
        ]) == 1, 2, sum))
    for (j in plolev) names <- c(names, paste(j, "x-absence", 
        sep = "", collapse = ""), paste(j, "x-presence", sep = "", 
        collapse = ""))
    colnames(alco) <- names
    alco
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
