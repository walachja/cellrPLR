\name{cellrPLR_biom}
\alias{cellrPLR_biom}
\alias{print.cellrPLR_biom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{cell-rPLR biomarker identification
%%  ~~function to do ... ~~
}
\description{Algorithm for cellwise outlier diagnostics
using robust pairwise log-ratios - biomarker identification. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cellrPLR_biom(data, type = "biweight", g1, g2, mainGroup = "max",
              biomarker, permutation = FALSE, B = 1000, p.alpha = 0.95)
## S3 method for class 'biomarker'
print(x, ...)              

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{dataset, either 'matrix' or 'data.frame'
%%     ~~Describe \code{data} here~~
}
  \item{type}{type of weighting function, possible values are \code{'biweight'}, \code{'huber'}, \code{'hampel'}
%%     ~~Describe \code{type} here~~
}
  \item{g1}{ vector with positions of samples from group 1
%%     ~~Describe \code{g1} here~~
}
  \item{g2}{ vector with positions of samples from group 2
%%     ~~Describe \code{g2} here~~
}
  \item{mainGroup}{ integer or character: group which is chosen as based. Possible values are: \code{1},\code{'1'} - group 1, \code{2}, \code{'2'} - group2, \code{'all'} - all samples, \code{'max'} - the bigger group
%%     ~~Describe \code{mainGroup} here~~
}
  \item{biomarker}{ variable names of known biomarkers (not necessary)
%%     ~~Describe \code{biomarker} here~~
}
  \item{permutation}{ logical value if permutation tests should be performed
%%     ~~Describe \code{permutation} here~~
}
  \item{B}{ number of permutation iterations, not used if \code{permutation = FALSE}
%%     ~~Describe \code{B} here~~
}
  \item{p.alpha}{cut-off  p-value for permutation tests, not used if \code{permutation = FALSE}
%%     ~~Describe \code{p.alpha} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of \code{Difference}, \code{Biomarker_results} and \code{Permutation_tests}.
\item{Difference}{Difference value is returned. The bigger the value, more important variable is for discrimination. In range <0;2>.}
\item{Biomarker_results}{Data frame with: \code{Ordered_position}, \code{Biomarker_names} and \code{Variable_position}.}
\item{Ordered_position}{Variable importance ordered position.}
\item{Biomarker_names}{If 'true' biomarkers were given as agrument \code{biomarker}, their names are returned.}
\item{Variable_position}{If 'true' biomarkers were given as agrument \code{biomarker}, their ordered positions are returned.}
\item{Permutation_tests}{Data frame with: \code{p_value} ,\code{TRUE_FALSE},\code{biomarkers},\code{biomarkers_names}.}
\item{p_value}{If \code{permutation = TRUE} p_value for each variable is returned.}
\item{TRUE_FALSE}{If \code{permutation = TRUE} infromation if variable is (TRUE) or is not (FALSE) identified as biomarker.}
\item{biomarkers_names}{If \code{permutation = TRUE} names of identified biomarkers.}
\item{statistics}{If \code{permutation = TRUE} and 'true' biomarkers were given: True positive (TP), True Positive Ratio (TPR), False Positive (FP) and False Positive Ratio (FPR) is given.}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{'Cellwise outlier detection and biomarker identification in metabolomics based on pairwise log-ratios', Walach J., Filzmoser P., Kouril S., submitted
%% ~put references to the literature/web site here ~
}
\author{Jan Walach <walach.jan@gmail.com>
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(452)
data <- gendata1_c()$X
colnames(data) <- LETTERS[1:9]
cellrPLR_biom(data, type = "biweight", g1 = 1:20 , g2 = 21:40, mainGroup = "all", biomarker=c('B','C'), permutation = TRUE, B = 100, p.alpha = 0.9)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~biomarker }% use one of  RShowDoc("KEYWORDS")
