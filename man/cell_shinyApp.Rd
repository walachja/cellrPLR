\name{cell_shinyApp}
\alias{cell_shinyApp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Shiny app for cellwise outlier diagnostics
%%  ~~function to do ... ~~
}
\description{Shiny app:Plot function of algorithm for cellwise outlier diagnostics
using robust pairwise log-ratios.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cell_shinyApp(data, g1, g2, mainGroup = "max", mainGroup2 = "all", simul = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{dataset, either 'matrix' or 'data.frame'
%%     ~~Describe \code{data} here~~
}
  \item{g1}{ vector with positions of samples from group 1
%%     ~~Describe \code{g1} here~~
}
  \item{g2}{ vector with positions of samples from group 2
%%     ~~Describe \code{g2} here~~
}
  \item{mainGroup}{ first plotly: integer or character: group which is chosen as based. Possible values are: \code{1},\code{'1'} - group 1, \code{2}, \code{'2'} - group2, \code{'all'} - all samples, \code{'max'} - the bigger group
%%     ~~Describe \code{mainGroup} here~~
}
  \item{mainGroup2}{ second plotly: integer or character: group which is chosen as based. Possible values are: \code{1},\code{'1'} - group 1, \code{2}, \code{'2'} - group2, \code{'all'} - all samples, \code{'max'} - the bigger group
%%     ~~Describe \code{mainGroup2} here~~
}
  \item{simul}{ compute cellwise outlier diagnostics simultaneously, might be hard on memory (RAM)
%%     ~~Describe \code{simul} here~~
}
}
\details{ For an easier exploration and understanding of the data, a Shiny app is included. Shiny , an
open source R package, is a web application and serves as an interactive tool for visualization. The Shiny app allows to
interactively apply different weighting and aggregation functions for the cell-rPLR algorithm, and supports zooming
into regions of the data matrix to see more details. Moreover, the variables can be interactively ordered based on their
importance for the group discrimination.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Shiny app is launched.
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
set.seed(741)
data <- gendata1_c()$X
colnames(data) <- LETTERS[1:9]
cell_shinyApp(data, g1 = 1:20 , g2 = 21:40,simul = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
