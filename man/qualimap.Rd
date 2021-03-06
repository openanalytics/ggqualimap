% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qualimap-tables.R
\name{qualimap}
\alias{qualimap}
\title{Extract Qualimap summary tables}
\usage{
qualimap(sample_info)
}
\arguments{
\item{sample_info}{Full path to file containing details about 
\code{samples} and their \code{paths}.}
}
\value{
An object of class \code{qualimap} which inherits from 
\code{"data.table"}, with two columns: \code{param} and \code{value}, where 
\code{value} is a list of \code{"data.table"}s.
}
\description{
\code{qualimap} function parses all the summary statistics of 
reports produced by \code{Qualimap} tool and returns a \code{data.table} 
with two columns: \code{param} and \code{value}. 

Each row of the \code{value} column contains the data corresponding to that 
\code{param}, and is itself a \code{data.table}.
}
\details{
The file provided to \code{sample_info} argument should contain 
at least these three columns: 
\itemize{
  \item \code{sample} - contains the \code{sample} name.
  \item \code{type} - one of \code{"summary"}, \code{"coverage_high"}, 
     \code{"coverage_low"} and \code{"coverage_total"} usually. See 
     vignette for an example `sample_info` file.

     For \code{html} files, \code{type} is \code{"summary"} and \code{"txt"} 
     files, one of the coverage types mentioned above should be provided.

     \bold{Note} that there may be other types of summary statistic file 
     generated by \code{Qualimap}. It is possible to provide a custom type 
     for those files, and \code{qualimap()} function would generate the 
     data, but plots are only implemented for summary tables in the 
     \code{html} file and all \code{coverage} files. This means you will 
     have to implement your own functions to plot those summary statistics.

  \item \code{path} - full path to the qualimap summary report, i.e., 
  \code{html} and/or \code{coverage_*.txt} files, for each sample.
  
  If just the file name (without path) is provided, it is assumed that the 
  file is in the same folder as the input file provided to 
  \code{sample_info} argument.
}

It can also optionally contain a \code{group} column. If present, the plots 
generated will take it into account and \code{color} / \code{facet} 
accordingly.
}
\examples{
path = system.file("tests/qualimap-sample", package="ggqualimap")
obj = qualimap(sample_info = file.path(path, "annotation.txt"))
}
\seealso{
\code{\link{plot_read_alignment}} \code{\link{plot_bias_profile}}
\code{\link{plot_coverage_profile}} \code{\link{plot_junction_analysis}}
\code{\link{plot_genomic_origin}}
}

