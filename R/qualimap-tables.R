#' @title Gather all Qualimap summary tables from all samples
#' 
#' @description This function extracts the \code{Input}, \code{Reads alignment}, 
#' \code{Reads genomics origin}, \code{Transcript coverage profile}, \code{Junction 
#' Analysis} and \code{coverage} tables from all samples.
#' 
#' @seealso \code{\link{root_node}}, \code{\link{child_node}}, \code{\link{html_table}}, \code{\link{qualimap_tables}}, \code{\link{tidy}}, \code{\link{extract_coverage_stats}}, , \code{\link{extract_summary_tables}}
#' @param dir Complete path to Qualimap reports for all samples from an experiment
#' @param groups Default is \code{NULL}. If the plots need to be grouped by samples, then a two 
#' column \code{data.table} has to be provided with names \code{sample} and \code{group}. The 
#' column \code{group} should have identical values for all samples belonging to the same group.
#' @export
qualimap <- function(dir, groups=NULL) {

	summary_files  = list.files(dir, pattern="qualimapReport.html$", full.names=TRUE, recursive=TRUE)
	samples = basename(dirname(summary_files))

	all_samples = lapply(seq_along(samples), function(i) 
		            qualimap_tables(samples[i], summary_files[i], groups))
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
#' @param groups Default is \code{NULL}. If the plots need to be grouped by samples, then a two 
#' column \code{data.table} has to be provided with names \code{sample} and \code{group}. The 
#' column \code{group} should have identical values for all samples belonging to the same group.
#' @export
#' @examples 
#' \dontrun{
#' exp_tables = qualimap_tables("qualimapReport.html", "coverage_profile_along_genes_(total).txt")
#' }
qualimap_tables <- function(sample_name = "", html_file, groups=NULL) {
	add_sample_and_group_cols <- function(x) {
		x[, sample_name := sample_name
		][, c("sample_group", "pairs") := split_sample(sample_name)]
	    if (is.data.table(groups)) {
      		x[groups, group := factor(i.group), on="sample_group"]
      	} else {
      		x[, group := NULL]
      	}
      	movecols = c("group", "sample_group", "sample_name", "pairs")
		setcolorder( x, c( movecols, head(names(x), -length(movecols)) ) )
	}
	doc = html_parse(html_file)
	tables = extract_summary_tables(doc)
	coverage_files = list.files(dirname(html_file), pattern="^coverage_.*", , full.names=TRUE, recursive=TRUE)
	coverage = lapply(coverage_files, extract_coverage_stats)
	setattr(coverage, 'names', paste(gsub(".*[(](.*)[)].*$", "\\1", coverage_files), "coverage", sep=" "))
	ans = c(tables, coverage)
	if (is.data.table(groups)) {
  		if (!all(c("sample", "group") %in% names(groups)))
   			stop("Columns 'sample' and 'group' should be present in argument 'groups'.")
		groups = copy(groups)[, sample_group := as.character(sample)]
	}
	ans = lapply(ans, add_sample_and_group_cols)
    ans
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

	setnames(origin, c("Region", "Read %"))
	origin[, Count := as.integer(gsub(",", "", gsub(" .*$", "\\1", as.character(`Read %`))))]
	origin[, Region := gsub(":$", "", Region)][, `Read %` := as.numeric(gsub(".*/[ ]*(.*)%$", "\\1", `Read %`))]

	alignment = rbind(alignment, junction[1L])
	setnames(alignment, c("Type", "Count"))
	alignment[, Type := gsub(":$", "", Type)][, Count := as.numeric(gsub(",", "", as.character(Count)))]
	alignment[, Percentage := Count/sum(origin$Count)]

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

split_sample <- function(sample_name) {
	splits = tstrsplit(sample_name, "_(?=[12]$)", perl=TRUE)
	if (length(splits) == 1L)
		splits = c(splits, list(1L))
	if (length(splits) != 2L)
		stop("'sample_name' argument should be of the form <samplename>_<pair> for paired end and <samplename> for single end Fastq files, with no '_' elsewhere.")
	splits[[2]] = as.integer(splits[[2]])
	splits
	# x[, c("sample_name", "pairs") := splits]
	# setcolorder(x, c("sample_name", "pairs", head(names(x), -2L)))
}

