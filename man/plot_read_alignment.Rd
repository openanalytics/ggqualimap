% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qualimap-plots.R
\name{plot_read_alignment}
\alias{plot_read_alignment}
\title{Plot Qualimap read alignment}
\usage{
plot_read_alignment(..., interactive = TRUE, geom = c("jitter", "box"))
}
\arguments{
\item{interactive}{logical, default is \code{TRUE}, which returns an 
\emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
static \code{ggplot2} plot.}

\item{geom}{Possible values are \code{"jitter"} (default) and \code{"box"}.
\code{"jitter"} is only possible for \code{interactive = TRUE}, and is 
usually the preferred option since it provides the lowest ink ratio, and 
contains the least amount of clutter. \code{"box"} is only possible when 
the input contains only one \code{qualimap} object.}

\item{\dots}{The set of \code{qualimap} objects to plot, usually of the 
form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
\code{examples}. The names will be used as title for facets.}
}
\value{
The plot object
}
\description{
This function plots the read alignment stats of samples. If 
info regarding the \code{group} each \code{sample} belongs to is also 
available, then the generated plot will take that into account to colour / 
facet accordingly.
}
\examples{
path = system.file("tests/qualimap-sample", package="ggqualimap")
obj = qualimap(sample_info = file.path(path, "annotation.txt"))

# interactive jitter
plot_read_alignment(sample = obj)
}
\seealso{
\code{\link{qualimap}} \code{\link{plot_bias_profile}}
\code{\link{plot_coverage_profile}} \code{\link{plot_junction_analysis}}
\code{\link{plot_genomic_origin}}
}

