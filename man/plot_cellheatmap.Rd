\name{plot_cellheatmap}
\alias{plot_cellheatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{cell-rPLR outlier identification
%%  ~~function to do ... ~~
}
\description{Plot function of algorithm for cellwise outlier diagnostics
using robust pairwise log-ratios.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_cellheatmap(data, type = "biweight", g1, g2, mainGroup = "max", 
             plotly = FALSE, grid=FALSE,title=NULL)
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
  \item{g2}{vector with positions of samples from group 2
%%     ~~Describe \code{g2} here~~
}
  \item{mainGroup}{ integer or character: group which is chosen as based. Possible values are: \code{1},\code{'1'} - group 1, \code{2}, \code{'2'} - group2, \code{'all'} - all samples, \code{'max'} - the bigger group
%%     ~~Describe \code{mainGroup} here~~
}
  \item{plotly}{ logical, should interactive plotly be used.
%%     ~~Describe \code{plotly} here~~
}
 \item{grid}{ logical, should grid be added to the plot.
%%     ~~Describe \code{plotly} here~~
}
 \item{title}{ title of plot
%%     ~~Describe \code{plotly} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ matrix with cellwise outlier information, in range <-1,1>.
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
plot_cellheatmap(data, type = "biweight", g1 = 1:20 , g2 = 21:40, mainGroup = "max",grid = TRUE,plotly = TRUE,title = 'Simulated example')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
