#' @title Extract Qualimap summary tables
#' 
#' @description \code{qualimap} function parses all the summary statistics of 
#' reports produced by \code{Qualimap} tool and returns a \code{data.table} 
#' with two columns: \code{param} and \code{value}. 
#' 
#' Each row of the \code{value} column contains the data corresponding to that 
#' \code{param}, and is itself a \code{data.table}.
#' 
#' @param sample_info Full path to file containing details about 
#' \code{samples} and their \code{paths}.
#' 
#' @details The file provided to \code{sample_info} argument should contain 
#' at least these three columns: 
#' \itemize{
#'   \item \code{sample} - contains the \code{sample} name.
#'   \item \code{type} - one of \code{"summary"}, \code{"coverage_high"}, 
#'      \code{"coverage_low"} and \code{"coverage_total"} usually. See 
#'      vignette for an example `sample_info` file.
#' 
#'      For \code{html} files, \code{type} is \code{"summary"} and \code{"txt"} 
#'      files, one of the coverage types mentioned above should be provided.
#' 
#'      \bold{Note} that there may be other types of summary statistic file 
#'      generated by \code{Qualimap}. It is possible to provide a custom type 
#'      for those files, and \code{qualimap()} function would generate the 
#'      data, but plots are only implemented for summary tables in the 
#'      \code{html} file and all \code{coverage} files. This means you will 
#'      have to implement your own functions to plot those summary statistics.
#' 
#'   \item \code{path} - full path to the qualimap summary report, i.e., 
#'   \code{html} and/or \code{coverage_*.txt} files, for each sample.
#'   
#'   If just the file name (without path) is provided, it is assumed that the 
#'   file is in the same folder as the input file provided to 
#'   \code{sample_info} argument.
#' }
#' 
#' It can also optionally contain a \code{group} column. If present, the plots 
#' generated will take it into account and \code{color} / \code{facet} 
#' accordingly.
#' 
#' @return An object of class \code{qualimap} which inherits from 
#' \code{"data.table"}, with two columns: \code{param} and \code{value}, where 
#' \code{value} is a list of \code{"data.table"}s.
#' @seealso \code{\link{plot_read_alignment}} \code{\link{plot_bias_profile}}
#' \code{\link{plot_coverage_profile}} \code{\link{plot_junction_analysis}}
#' \code{\link{plot_genomic_origin}}
#' @examples
#' path = system.file("tests/qualimap-sample", package="ggqualimap")
#' obj = qualimap(sample_info = file.path(path, "annotation.txt"))
#' @export
qualimap <- function(sample_info) {

    group=path=type=NULL
    info = fread(sample_info, colClasses=list(character=c("sample")))
    cols = c("sample", "type", "path")
    rest = setdiff(cols, names(info))
    if (length(rest)) stop("Columns [", paste(rest, collapse=","), 
        "] not found in file provided to argument 'sample_info'.")
    if (is.null(info[["group"]])) info[, "group" := ""]
    samples = info[["sample"]]
    # set paths correctly, if necessary
    idx = which(info[["path"]] == basename(info[["path"]]))
    if (length(idx)) 
        info[(idx), "path" := file.path(dirname(sample_info), path)][]
    info[, "tables" := list(list(qualimap_tables(path, sample, type, group))), 
                by=1:nrow(info)]
    merge_tables <- function(ll) {
        nm = names(ll[[1L]])
        fans = lapply(seq_along(nm), function(i) 
                rbindlist(lapply(ll, `[[`, i)))
        setDT(list(param = nm, value = fans))
    }
    ans = info[, merge_tables(tables), by="type"
             ][, "type" := NULL][]
    setattr(ans, 'class', c('qualimap', class(ans)))
}

#' @title Extract all summary tables for a particular sample
#' 
#' @description \bold{NOTE:} For internal use only.
#' 
#' This is the workhorse that powers \code{qualimap}. Given a \code{sample}, 
#' \code{type} and \code{path}, \code{extract_summary_tables} extracts all 
#' summary statistics and returns them as a list of \code{data.table}s.
#' 
#' @return A list of data.tables.
#' @param file The file corresponding to \emph{that \code{sample}}.
#' @param sample Sample name.
#' @param type  \code{"summary"} or one of \code{"coverage_.*"} values. See 
#' \code{vignette} for details.
#' @param group Group / condition this sample belongs to.
qualimap_tables <- function(file, sample, type, group) {

    ext = tools::file_ext(file) # html or txt / non-html?
    if (!ext %in% c("html", "txt")) {
        stop("File extensions must be either 'html' or 'txt'.")
    }
    ans = switch(ext, 
        html = html_summary_tables(file, sample, group),
        txt  = txt_summary_tables(file, sample, group, type)
    )
    ans
}

# @title Extract all summary tables from html file
html_summary_tables <- function(file, sample, group) {

    doc = XML::htmlParse(file) # parse html content
    tbl = XML::readHTMLTable(doc, header=FALSE) # extract as list of DFs
    tbl = lapply(tbl, setDT) # convert each element to DT
    # have extracted tables, need to extract what type of table to set as names
    hdr = XML::xpathSApply(XML::xmlRoot(doc),"//div[@class='table-summary']/h3")
    hdr = sapply(hdr, xmlValue)
    hdr = gsub("[ ]+", "_", tolower(hdr))
    setattr(tbl, "names", hdr)
    setattr(tbl, 'class', 'html')
    ans = tidy(tbl)
    nm = c("sample", "pairs", "group")
    nm_vals = c(split_sample(sample), list(group))
    ans = lapply(ans, function(x) {
        x[, c(nm) := nm_vals][]
        setcolorder(x, c(nm, setdiff(names(x), nm)))
    })
}

txt_summary_tables <- function(file, sample, group, type) {
    ans = list(data.table::fread(file))
    setattr(ans, "names", type)
    setattr(ans, "class", "txt")
    ans = tidy(ans)
    nm = c("sample", "pairs", "group")
    nm_vals = c(split_sample(sample), list(group))
    ans = lapply(ans, function(x) {
        x[, c(nm) := nm_vals][]
        setcolorder(x, c(nm, setdiff(names(x), nm)))
    })
}

# @title Cleanup all summary tables
tidy <-function(x) {
    UseMethod("tidy")
}

tidy.default <- function(x) {
    stop("No default 'tidy' method found for object of class ", class(x))
}

tidy.html <- function(x) {

    param=alignment=read_count=V2=region=percentage=bias=value=NULL
    # spaces to underscores, lower case, remove ":"
    fix_text <- function(x) {
        x = gsub("[ ]+", "_", tolower(x))
        gsub(":$", "", x)
    }
    fix_numbers <- function(x, type=c("int", "numeric")) {
        type = match.arg(type)
        x = gsub(",", "", x)
        ans = switch(type, 
            int = as.integer(x), 
            numeric = as.numeric(x)
        )
        ans
    }
    fix_percent <- function(x) {
        as.numeric(gsub("%$", "", x))
    }

    params = c("input", "reads_alignment", "reads_genomic_origin", 
                "transcript_coverage_profile", "junction_analysis")

    if (params[1L] %in% names(x)) {
        input = .subset2(x, params[1L])
        stopifnot(is.data.table(input))
        setnames(input, c("param", "value"))
        input[, "param" := fix_text(param)]
    }
    if (params[2L] %in% names(x)) {
        align = .subset2(x, params[2L])
        setnames(align, c("alignment", "read_count"))
        align[, "alignment" := fix_text(alignment)
            ][, "read_count" := fix_numbers(read_count)]
    }
    if (params[3L] %in% names(x)) {
        geno = .subset2(x, params[3L])
        geno[, c("read_count", "percentage") := tstrsplit(V2, "[ ]*/[ ]*")]
        set(geno, j="V2", value=NULL)
        setnames(geno, c("region", "read_count", "percentage"))
        geno[, "region" := fix_text(region)
           ][, "read_count" := fix_numbers(read_count)
           ][, "percentage" := fix_percent(percentage)]
    }
    if (params[4L] %in% names(x)) {
        cov = .subset2(x, params[4L])
        setnames(cov, c("bias", "value"))
        cov[, "bias" := fix_text(bias)
            ][, "value" := fix_numbers(value, "numeric")]
    }
    if (params[5L] %in% names(x)) {
        junc = .subset2(x, params[5L])
        junc_cnt = data.table(alignment = "at_junction", 
                                read_count = fix_numbers(junc[1, V2]))
        junc = junc[-1L, ]
        setnames(junc, c("junc_type", "percentage"))
        junc[, "percentage" := fix_percent(percentage)]
    }
    align = rbind(align, junc_cnt)
    align[, "percentage" := read_count / sum(geno[["read_count"]])][]
    ans = list(input, align, geno, cov, junc)
    setattr(ans, "names", params)
}

tidy.txt <- function(x) {
    fix_text <- function(x) {
        x = gsub("[#]", "", tolower(x))
        gsub("[ ]+", "_", x)
    }
    ans = lapply(x, function(th) setnames(th, fix_text(names(th))))
}

split_sample <- function(sample_name) {
    splits = tstrsplit(sample_name, "_(?=[12]$)", perl=TRUE)
    if (length(splits) == 1L)
        splits = c(splits, list(1L))
    if (length(splits) != 2L) {
        stop("'sample_name' argument should be of the form <samplename>_<pair> 
                for paired end and <samplename> for single end Fastq files, 
                with no '_' elsewhere.")
    }
    splits[[2]] = as.integer(splits[[2]])
    splits
}

# #' @title Wrapper to extract sub-nodes within XMLNode objects
# child_node <- function(x) {
#     if (!inherits(x, "XMLInternalElementNode"))
#         stop("Expecting object of class 'XMLInternalElementNode' as input.")
#     XML::xmlChildren(x)
# }
