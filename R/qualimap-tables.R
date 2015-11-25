#' @title Gather all Qualimap summary tables from all samples
#' 
#' @description This function extracts the \code{Input}, \code{Reads alignment}, 
#' \code{Reads genomics origin}, \code{Transcript coverage profile}, \code{Junction 
#' Analysis} and \code{coverage} tables from all samples.
#' 
#' @seealso \code{\link{root_node}}, \code{\link{child_node}}, \code{\link{html_table}}, \code{\link{qualimap_tables}}, \code{\link{tidy}}, \code{\link{extract_coverage_stats}}, , \code{\link{extract_summary_tables}}
#' @param dir Complete path to Qualimap reports for all samples from an experiment
#' @export
qualimap <- function(dir) {

	summary_files  = list.files(dir, pattern="qualimapReport.html$", full.names=TRUE, recursive=TRUE)
	coverage_files = list.files(dir, pattern="coverage_.*total)\\.txt$", full.names=TRUE, recursive=TRUE)
	samples = basename(dirname(summary_files))

	all_samples = lapply(seq_along(samples), function(i) qualimap_tables(samples[i], summary_files[i], coverage_files[i]))
	ans = lapply(seq_along(all_samples[[1L]]), function(i) rbindlist(lapply(all_samples, `[[`, i)))
	setattr(ans, 'names', names(all_samples[[1L]]))
}

#' @title Extract all summary tables from Qualimap report
#' 
#' @description This function extracts the \code{Input}, \code{Reads alignment}, 
#' \code{Reads genomics origin}, \code{Transcript coverage profile}, \code{Junction 
#' Analysis} tables.
#' 
#' @seealso \code{\link{root_node}}, \code{\link{child_node}}, \code{\link{html_table}}, \code{\link{qualimap_tables}}, \code{\link{tidy}}, \code{\link{extract_coverage_stats}}, , \code{\link{extract_summary_tables}}, \code{qualimap}
#' @param sample_name Sample name to be appeneded as a separate column to the summary/coverage tables.
#' @param html_file Complete path to html file from Qualimap.
#' @param coverage_file Complete path to coverage stats from Qualimap.
#' @export
#' @examples 
#' \dontrun{
#' exp_tables = qualimap_tables("qualimapReport.html", "coverage_profile_along_genes_(total).txt")
#' }
qualimap_tables <- function(sample_name = "", html_file, coverage_file) {
	add_sample_col <- function(x) {
		x[, sample_name := sample_name]
		setcolorder(x, c("sample_name", head(names(x), -1L)))
	}
	doc = html_parse(html_file)
	tables = extract_summary_tables(doc)
	coverage = extract_coverage_stats(coverage_file)
	ans = c(tables, `Total coverage` = list(coverage))
	ans = lapply(ans, add_sample_col)
}

#' @title Wrapper for parsing HTML content
#' @description Wrapper for \code{XML::htmlParse}. For internal use only.
#' @seealso \code{\link{root_node}}, \code{\link{child_node}}, \code{\link{html_table}}, \code{\link{qualimap_tables}}, \code{\link{extract_summary_tables}}, \code{\link{tidy}}, \code{\link{extract_coverage_stats}}, \code{qualimap}
#' @param html_file Complete path to html file from Qualimap.
html_parse <- function(html_file) {
	if (!identical("html", tools::file_ext(html_file)))
		stop("Expecting a html file as input.")
	XML::htmlParse(html_file)
}

#' @title Extract all summary tables from Qualimap
#' @description The function extracts just the summary tables from a Qualimap report file on to 
#' a list of \code{data.table}s and tidies them up. For internal use only.
#' @seealso \code{\link{root_node}}, \code{\link{child_node}}, \code{\link{html_table}}, \code{\link{qualimap_tables}}, \code{\link{tidy}}, \code{\link{extract_coverage_stats}}, \code{qualimap}
#' @param doc An object of class \code{HTMLInternalDocument} obtained from \code{html_parse}.
extract_summary_tables <- function(doc) {
	tables = lapply(html_table(doc), setDT)
	# navigate to h3 element for div elements with table-summary class
	# and extract the values. These are the table headings.
	contents = xpathSApply(root_node(doc), "//div[@class='table-summary']/h3")
	setattr(tables, 'names', sapply(contents, xmlValue))
	tidy(tables)
}

#' @title Extract all tables into list of data frames
#' @description Wrapper for \code{XML::readHTMLTable}. For internal use only.
#' @seealso \code{\link{html_parse}}, \code{\link{root_node}}, \code{\link{child_node}}, \code{\link{qualimap_tables}}, \code{\link{extract_summary_tables}}, \code{\link{tidy}}, \code{\link{extract_coverage_stats}}, \code{qualimap}
#' @param x An object of class \code{HTMLInternalDocument} obtained from \code{html_parse}. 
html_table <- function(x) {
	if (!"HTMLInternalDocument" %chin% class(x))
		stop("Expecting object of class 'HTMLInternalDocument' as input.")
	XML::readHTMLTable(x, header=FALSE)
}

#' @title Cleanup all summary tables
#' @description For internal use only.
#' @seealso \code{\link{root_node}}, \code{\link{child_node}}, \code{\link{html_table}}, \code{\link{qualimap_tables}}, \code{\link{extract_summary_tables}}, \code{\link{extract_coverage_stats}}, \code{qualimap}
#' @param html_file Complete path to html file from Qualimap.
tidy <- function(x) {
	origin    = x[["Reads genomic origin"]]
	coverage  = x[["Transcript coverage profile"]]
	junction  = x[["Junction analysis"]]
	alignment = x[["Reads alignment"]]
	input     = x[["Input"]]

	alignment = rbind(alignment, junction[1L])
	setnames(alignment, c("Type", "Count"))
	alignment[, Type := gsub(":$", "", Type)][, Count := as.integer(gsub(",", "", as.character(Count)))]

	setnames(origin, c("Region", "Read %"))
	origin[, Region := gsub(":$", "", Region)][, `Read %` := as.numeric(gsub(".*/[ ]*(.*)%$", "\\1", `Read %`))]

	setnames(coverage, c("Position", "Value"))
	coverage[, Position := gsub(":$", "", Position)][, Value := as.numeric(as.character(Value))]

	junction = junction[-1L]
	setnames(junction, c("Junction", "Read %"))
	junction[, Junction := gsub(":$", "", Junction)][, `Read %` := as.numeric(gsub("%$", "", `Read %`))]

	setnames(input, c("Detail", "Value"))
	list(Input = input, 
		`Reads genomic origin` = origin, 
		 `Transcript coverage profile` = coverage,
		 `Junction analysis` = junction,
		 `Reads alignment` = alignment)
}

#' @title Wrapper to extract top-level XML node
#' @description Wrapper for \code{XML::xmlRoot}. For internal use only.
#' @seealso \code{\link{html_parse}}, \code{\link{child_node}}, \code{\link{html_table}}, \code{\link{qualimap_tables}}, \code{\link{extract_summary_tables}}, \code{\link{tidy}}, \code{\link{extract_coverage_stats}}, \code{qualimap}
#' @param doc An object of class \code{HTMLInternalDocument} obtained from \code{html_parse}. 
root_node <- function(doc) {
	if (!"HTMLInternalDocument" %chin% class(doc))
		stop("Expecting object of class 'HTMLInternalDocument' as input.")
	XML::xmlRoot(doc)
}

#' @title Wrapper to extract sub-nodes within XMLNode objects
#' @description Wrapper for \code{XML::xmlChildren}. For internal use only.
#' @seealso \code{\link{html_parse}}, \code{\link{root_node}}, \code{\link{html_table}}, \code{\link{qualimap_tables}}, \code{\link{extract_summary_tables}}, \code{\link{tidy}}, \code{\link{extract_coverage_stats}}, \code{qualimap}
#' @param x An object of class \code{XMLInternalElementNode} obtained from \code{root_node}. 
child_node <- function(x) {
	if (!"XMLInternalElementNode" %chin% class(x))
		stop("Expecting object of class 'XMLInternalElementNode' as input.")
	XML::xmlChildren(x)
}

#' @title Load Qualimap coverage stats
#' @description Loads the coverage stats. For internal use only.
#' @seealso \code{\link{html_parse}}, \code{\link{root_node}}, \code{\link{html_table}}, \code{\link{qualimap_tables}}, \code{\link{extract_summary_tables}}, \code{\link{tidy}}, \code{qualimap}
#' @param coverage_file A complete path to the text file containing coverage stats.
extract_coverage_stats <- function(coverage_file) {
	fread(coverage_file)
}
