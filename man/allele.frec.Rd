\name{allele.frec}
\Rdversion{0.0.1}
\alias{allele.frec}
\title{Calculation of Allele frequency with Dominant markers and mixed ploidy levels 
%%  ~~function to do ... ~~
}
\description{
This function aims at elaluating allele frequency in a ballance population, this considering the ploidy level. 
Then It calculated an weighted avarage 
}
\usage{
allele.frec(xx, ploidy = 2)
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
There is a little biais on small sized sample. 
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
    glob <- NULL
    for (i in plolev) alco <- cbind(alco, (1 - (apply(as.matrix(xx[ploidy == 
        i, ]) == 0, 2, sum)/(apply(as.matrix(xx[ploidy == i, 
        ]) == 0, 2, sum) + apply(as.matrix(xx[ploidy == i, ]) == 
        1, 2, sum)))^(1/i)) * 100)
    for (j in plolev) names <- c(names, paste(j, "x-frequence(\%)", 
        sep = "", collapse = ""))
    colnames(alco) <- names
    for (k in plolev) glob <- c(glob, nrow(xx[ploidy == k, ]))
    alco <- cbind(alco, apply(alco * glob/sum(glob), 1, sum))
    colnames(alco) <- c(names, "Frequency")
    alco
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
