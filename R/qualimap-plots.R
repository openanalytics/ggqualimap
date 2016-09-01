#' @title Plot Qualimap read alignment 
#' 
#' @description This function plots the read alignment stats of samples. If 
#' info regarding the \code{group} each \code{sample} belongs to is also 
#' available, then the generated plot will take that into account to colour / 
#' facet accordingly.
#' 
#' @param \dots The set of \code{qualimap} objects to plot, usually of the 
#' form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
#' \code{examples}. The names will be used as title for facets.
#' @param interactive logical, default is \code{TRUE}, which returns an 
#' \emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
#' static \code{ggplot2} plot.
#' @param geom Possible values are \code{"jitter"} (default) and \code{"box"}.
#' \code{"jitter"} is only possible for \code{interactive = TRUE}, and is 
#' usually the preferred option since it provides the lowest ink ratio, and 
#' contains the least amount of clutter. \code{"box"} is only possible when 
#' the input contains only one \code{qualimap} object.
#' 
#' @seealso \code{\link{qualimap}} \code{\link{plot_bias_profile}}
#' \code{\link{plot_coverage_profile}} \code{\link{plot_junction_analysis}}
#' \code{\link{plot_genomic_origin}}
#' @export
#' @examples
#' require(ggqualimap)
#' path = system.file("tests/qualimap-sample", package="ggqualimap")
#' obj = qualimap(sample_info = file.path(path, "annotation.txt"))
#' 
#' # interactive jitter
#' plot_read_alignment(sample = obj)
plot_read_alignment <- function(..., interactive=TRUE, 
            geom=c("jitter", "box")) {

    geom = match.arg(geom)
    ll = list(...)
    if (length(ll) > 1L && geom == "box") {
        warning("box plots are only possible on one qualimap object.", 
            " Setting geom to 'jitter' and 'interactive' to TRUE.")
        geom = "jitter"
        interactive = TRUE
    }
    ra = lapply(ll, function(l) {
            stopifnot(inherits(l, "qualimap"))
            ans = l[param == "reads_alignment", value]
            data.table::rbindlist(ans)
        })
    if (is.null(names(ll)))
        setattr(ra, 'names', paste0("fastqc_obj", seq_along(ll)))
    ra = data.table::rbindlist(ra, idcol=TRUE)
    ra[, "percent" := sprintf("%.4f", percentage)] # for tooltip
    idx = c(grep("^[Aa]ligned pairs", ra[["alignment"]], ignore.case=TRUE), 
            grep("^[tT]otal number of ", ra[["alignment"]], ignore.case=TRUE), 
            grep("[Ss]econdary align", ra[["alignment"]], ignore.case=TRUE))
    if (length(idx)) ra = ra[!(idx)]

    as_factor <- function(x) factor(x, levels=unique(x))
    cols = c("sample", "alignment", "group")
    ra[, (cols) := lapply(.SD, as_factor), .SDcols=cols][]
    aes = list(
            x = "group", 
            y = "percentage", 
            sample_name = "sample", 
            fill = "group", 
            facet_lhs = ".id", 
            facet_rhs = "alignment", 
            tooltip1 = "percent"
        )
    theme = list(
            xlab = "Groups", 
            ylab = "Percentage", 
            title = "Read alignment stats"
        )
    pl = switch(geom, 
            jitter = qualimap_jitter(ra, aes, theme, interactive), 
            box = qualimap_box(ra, aes, theme, interactive)
        )
    pl
}

#' @title Plot Qualimap reads across genomic origin 
#' 
#' @description This function plots reads across genomic origin. If 
#' info regarding the \code{group} each \code{sample} belongs to is also 
#' available, then the generated plot will take that into account to colour / 
#' facet accordingly.
#' 
#' @param \dots The set of \code{qualimap} objects to plot, usually of the 
#' form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
#' \code{examples}. The names will be used as title for facets.
#' @param interactive logical, default is \code{TRUE}, which returns an 
#' \emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
#' static \code{ggplot2} plot.
#' @param geom Possible values are \code{"jitter"} (default) and \code{"bar"}.
#' \code{"jitter"} is only possible for \code{interactive = TRUE}, and is 
#' usually the preferred option since it provides the lowest ink ratio, and 
#' contains the least amount of clutter. \code{"bar"} is only possible when 
#' the input contains only one \code{qualimap} object.
#' 
#' @seealso \code{\link{qualimap}} \code{\link{plot_bias_profile}}
#' \code{\link{plot_coverage_profile}} \code{\link{plot_junction_analysis}}
#' \code{\link{plot_read_alignment}}
#' @export
#' @examples
#' require(ggqualimap)
#' path = system.file("tests/qualimap-sample", package="ggqualimap")
#' obj = qualimap(sample_info = file.path(path, "annotation.txt"))
#' 
#' # interactive jitter
#' plot_genomic_origin(sample = obj)
#' 
#' interactive bar plot
#' plot_genomic_origin(sample = obj, geom="bar")
#' 
#' non-interactive bar plot
#' plot_genomic_origin(sample = obj, geom="bar", interactive=FALSE)
plot_genomic_origin <- function(..., interactive=TRUE, 
                            geom=c("jitter", "bar")) {

    geom = match.arg(geom)
    ll = list(...)
    if (length(ll) > 1L && geom == "bar") {
        warning("bar plots are only possible on one qualimap object.", 
            " Setting geom to 'jitter' and 'interactive' to TRUE.")
        geom = "jitter"
        interactive = TRUE
    }
    go = lapply(ll, function(l) {
            stopifnot(inherits(l, "qualimap"))
            ans = l[param == "reads_genomic_origin", value]
            data.table::rbindlist(ans)
        })
    if (is.null(names(ll)))
        setattr(go, 'names', paste0("fastqc_obj", seq_along(ll)))
    go = data.table::rbindlist(go, idcol=TRUE)
    go[, "percent" := sprintf("%.4f", percentage)] # for tooltip

    as_factor <- function(x) factor(x, levels=unique(x))
    cols = c("sample", "region")
    go[, (cols) := lapply(.SD, as_factor), .SDcols=cols][]
    aes = list(
            x = if (geom == "jitter") "region" else "sample", 
            y = "percentage", 
            sample_name = "sample", 
            fill = if (geom == "jitter") "group" else "region", 
            group = "sample", 
            position = "stack", 
            facet_lhs = ".id", 
            facet_rhs = ".", 
            x_breaks_discrete = "region", 
            tooltip1 = "percent"
        )
    theme = list(
            xlab = "Sample",
            ylab = "Read %", 
            title = "Read across genomic origin"
        )
    pl = switch(geom, 
            jitter = qualimap_jitter2(go, aes, theme, interactive), 
            bar = qualimap_bar(go, aes, theme, interactive)
        )
    pl
}

#' @title Plot Qualimap bias profile
#' 
#' @description This function plots the bias profiles of samples. If 
#' info regarding the \code{group} each \code{sample} belongs to is also 
#' available, then the generated plot will take that into account to colour / 
#' facet accordingly.
#' 
#' @param \dots The set of \code{qualimap} objects to plot, usually of the 
#' form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
#' \code{examples}. The names will be used as title for facets.
#' @param interactive logical, default is \code{TRUE}, which returns an 
#' \emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
#' static \code{ggplot2} plot.
#' @param geom Possible values are \code{"jitter"} (default) and \code{"bar"}.
#' \code{"jitter"} is only possible for \code{interactive = TRUE}, and is 
#' usually the preferred option since it provides the lowest ink ratio, and 
#' contains the least amount of clutter. \code{"bar"} is only possible when 
#' the input contains only one \code{qualimap} object.
#' 
#' @seealso \code{\link{qualimap}} \code{\link{plot_read_alignment}}
#' \code{\link{plot_coverage_profile}} \code{\link{plot_junction_analysis}}
#' \code{\link{plot_genomic_origin}}
#' @export
#' @examples
#' require(ggqualimap)
#' path = system.file("tests/qualimap-sample", package="ggqualimap")
#' obj = qualimap(sample_info = file.path(path, "annotation.txt"))
#' 
#' # interactive jitter
#' plot_bias_profile(sample = obj)
#' interactive bar plot
#' plot_bias_profile(sample = obj, geom="bar")
#' 
#' non-interactive bar plot
#' plot_bias_profile(sample = obj, geom="bar", interactive=FALSE)
plot_bias_profile <- function(..., interactive=TRUE, geom=c("jitter", "bar")) {

    geom = match.arg(geom)
    ll = list(...)
    if (length(ll) > 1L && geom == "bar") {
        warning("bar plots are only possible on one qualimap object.", 
            " Setting geom to 'jitter' and 'interactive' to TRUE.")
        geom = "jitter"
        interactive = TRUE
    }
    bias = lapply(ll, function(l) {
            stopifnot(inherits(l, "qualimap"))
            ans = l[param == "transcript_coverage_profile", value]
            data.table::rbindlist(ans)
        })
    if (is.null(names(ll)))
        setattr(bias, 'names', paste0("fastqc_obj", seq_along(ll)))
    bias = data.table::rbindlist(bias, idcol=TRUE)

    as_factor <- function(x) factor(x, levels=unique(x))
    cols = c("sample", "bias")
    bias[, (cols) := lapply(.SD, as_factor), .SDcols=cols][]
    aes = list(
            x = if (geom == "jitter") "bias" else "sample", 
            y = "value", 
            sample_name = "sample", 
            fill = if (geom == "jitter") "group" else "bias", 
            group = if (geom == "jitter") "sample" else "bias", 
            position = position_dodge(), 
            facet_lhs = ".id", 
            facet_rhs = ".", 
            x_breaks_discrete = "bias", 
            tooltip1 = "value"
        )
    theme = list(
            xlab = "Sample",
            ylab = "Value", 
            title = "Bias profiles"
        )
    pl = switch(geom, 
            jitter = qualimap_jitter2(bias, aes, theme, interactive), 
            bar = qualimap_bar(bias, aes, theme, interactive)
        )
    pl
}

#' @title Plot Qualimap coverage profiles
#' 
#' @description This function plots the coverage profiles of samples. If 
#' info regarding the \code{group} each \code{sample} belongs to is also 
#' available, then the generated plot will take that into account to colour / 
#' facet accordingly.
#' 
#' @param \dots The set of \code{qualimap} objects to plot, usually of the 
#' form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
#' \code{examples}. The names will be used as title for facets.
#' @param interactive logical, default is \code{TRUE}, which returns an 
#' \emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
#' static \code{ggplot2} plot.
#' @param geom Only possible value is \code{"line"}.
#' 
#' @seealso \code{\link{qualimap}} \code{\link{plot_read_alignment}}
#' \code{\link{plot_bias_profile}} \code{\link{plot_junction_analysis}}
#' \code{\link{plot_genomic_origin}}
#' @export
#' @examples
#' require(ggqualimap)
#' path = system.file("tests/qualimap-sample", package="ggqualimap")
#' obj = qualimap(sample_info = file.path(path, "annotation.txt"))
#' 
#' # interactive
#' plot_coverage_profile(sample = obj)
#' 
#' # non-interactive
#' plot_coverage_profile(sample = obj, interactive=FALSE)
plot_coverage_profile <- function(..., interactive=TRUE, geom="line") {

    geom = match.arg(geom)
    ll = list(...)
    # if (length(ll) > 1L && !interactive) {
    #     warning("Static plot is only possible on one qualimap object.", 
    #         " Setting 'interactive' to TRUE.")
    #     interactive = TRUE
    # }
    cov = lapply(ll, function(l) {
            stopifnot(inherits(l, "qualimap"))
            vals = grep("^coverage_.*$", l$param, value=TRUE)
            ans = l[param %in% vals, value]
            if (length(ans)) {
                data.table::setattr(ans, 'names', vals)
                data.table::rbindlist(ans, idcol="coverage_type")
            }
        })
    if (is.null(names(ll)))
        setattr(cov, 'names', paste0("fastqc_obj", seq_along(ll)))
    cov = data.table::rbindlist(cov, idcol=TRUE)

    as_factor <- function(x) factor(x, levels=unique(x))
    cols = c("sample", "coverage_type", ".id")
    cov[, (cols) := lapply(.SD, as_factor), .SDcols=cols][]
    aes = list(
            x = "transcript_position", 
            y = "transcript_coverage_profile", 
            sample_name = "sample", 
            colour = if (interactive && length(ll)>1L) ".id" else "sample", 
            facet_lhs = "coverage_type", 
            facet_rhs = if (length(ll) == 1L) NULL else ".id", 
            tooltip1 = ".id"
        )
    theme = list(
            xlab = "Position of transcript",
            ylab = "Coverage", 
            title = "Coverage profiles"
        )

    p = ggplot(cov, aes_string(x=aes$x, y=aes$y, sample_name=aes$sample_name, 
            tooltip1=aes$tooltip1)) + 
        geom_line(aes_string(colour=aes$colour)) + 
        xlab(theme$xlab) + ylab(theme$ylab) + theme_bw() + 
        facet_wrap(reformulate(aes$facet_lhs, aes$facet_rhs), ncol=1L, 
            scales="free_y") + 
        theme(legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06), 
            legend.position="bottom") + 
        ggtitle(theme$title)
    if (interactive) {
        if (!requireNamespace("plotly"))
            stop("Package 'plotly' is not available.")
        p = plotly::ggplotly(p, tooltip=c("x", "y", "sample_name", "tooltip1"))
    }
    p
}

#' @title Plot Qualimap junction analysis
#' 
#' @description This function plots the junction stats for samples. If 
#' info regarding the \code{group} each \code{sample} belongs to is also 
#' available, then the generated plot will take that into account to colour / 
#' facet accordingly.
#' 
#' @param \dots The set of \code{qualimap} objects to plot, usually of the 
#' form \code{sample_name_1 = obj1}, \code{sample_name_2 = obj2}, etc. See 
#' \code{examples}. The names will be used as title for facets.
#' @param interactive logical, default is \code{TRUE}, which returns an 
#' \emph{interactive} \code{plotly} plot. If \code{FALSE}, it returns a 
#' static \code{ggplot2} plot.
#' @param geom Only possible value is \code{"bar"}.
#' 
#' @seealso \code{\link{qualimap}} \code{\link{plot_read_alignment}}
#' \code{\link{plot_bias_profile}} \code{\link{plot_coverage_profile}}
#' \code{\link{plot_genomic_origin}}
#' @export
#' @examples
#' require(ggqualimap)
#' path = system.file("tests/qualimap-sample", package="ggqualimap")
#' obj = qualimap(sample_info = file.path(path, "annotation.txt"))
#' 
#' # interactive
#' plot_junction_analysis(sample = obj)
#' 
#' # non-interactive
#' plot_junction_analysis(sample = obj, interactive=FALSE)
plot_junction_analysis <- function(..., interactive=TRUE, geom="bar") {

    geom = match.arg(geom)
    ll = list(...)
    # if (length(ll) > 1L && !interactive) {
    #     warning("Static plot is only possible on one qualimap object.", 
    #         " Setting 'interactive' to TRUE.")
    #     interactive = TRUE
    # }
    jun = lapply(ll, function(l) {
            stopifnot(inherits(l, "qualimap"))
            ans = l[param == "junction_analysis", value]
            data.table::rbindlist(ans, idcol="coverage_type")
        })
    if (is.null(names(ll)))
        setattr(jun, 'names', paste0("fastqc_obj", seq_along(ll)))
    jun = data.table::rbindlist(jun, idcol=TRUE)

    as_factor <- function(x) factor(x, levels=unique(x))
    cols = c("sample", "group", ".id", "junc_type")
    jun[, (cols) := lapply(.SD, as_factor), .SDcols=cols][]
    aes = list(
            x = "junc_type", 
            y = "percentage", 
            sample_name = "sample", 
            fill = ".id", 
            facet_lhs = "sample", 
            facet_rhs = if (length(ll) == 1L) NULL else ".id", 
            tooltip1 = ".id"
        )
    theme = list(
            xlab = "Junction",
            ylab = "Percentage", 
            title = "Junction stats"
        )

    p = ggplot(jun, aes_string(x=aes$x, y=aes$y, sample_name=aes$sample_name, 
            tooltip1=aes$tooltip1, fill=aes$fill)) + 
        geom_bar(stat="identity") + 
        xlab(theme$xlab) + ylab(theme$ylab) + theme_bw() + 
        facet_wrap(reformulate(aes$facet_lhs, aes$facet_rhs), ncol=1L, 
            scales="free_y") + 
        theme(axis.text.x=element_text(angle=45, hjust=1), 
                text=element_text(size=16L), 
            legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06), 
            legend.position="bottom") + 
        ggtitle(theme$title)
    if (interactive) {
        if (!requireNamespace("plotly"))
            stop("Package 'plotly' is not available.")
        p = plotly::ggplotly(p, tooltip=c("x", "y", "sample_name", "tooltip1"))
    }
    p
}

## internal geoms ----------------------- 

qualimap_jitter <- function(dt, aes, theme, interactive=TRUE) {

    if (!interactive) {
        warn = paste("When geom='jitter', only interactive", 
                "plots are possible. Setting interactive=TRUE")
        warning(warn)
        interactive=TRUE
    }
    if (!requireNamespace("plotly"))
        stop("Package 'plotly' is not available.")
    p = ggplot(dt, aes_string(x=aes$x, y=aes$y, sample=aes$sample_name, 
                    tooltip1=aes$tooltip1)) + 
        geom_point(position = position_jitter(), aes_string(fill=aes$fill), 
            shape=21L, size=4L) + 
        xlab(theme$xlab) + ylab(theme$ylab) + theme_bw() + 
        (if (length(unique(dt$`.id`)) == 1L) 
            facet_wrap(reformulate(aes$facet_rhs), scales="free") 
        else facet_grid(reformulate(aes$facet_rhs, aes$facet_lhs), 
                scales="free", space="free")) + 
        scale_x_discrete(breaks=NULL) + 
        theme(legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06), 
            legend.position="none") + 
        ggtitle(theme$title)
    plotly::ggplotly(p, tooltip = c("sample", "tooltip1"))
}

qualimap_jitter2 <- function(dt, aes, theme, interactive=TRUE) {

    if (!interactive) {
        warn = paste("When geom='jitter', only interactive", 
                "plots are possible. Setting interactive=TRUE")
        warning(warn)
        interactive=TRUE
    }
    if (!requireNamespace("plotly"))
        stop("Package 'plotly' is not available.")
    p = ggplot(dt, aes_string(x=aes$x, y=aes$y, sample=aes$sample_name, 
                    tooltip1=aes$tooltip1)) + 
        geom_point(position = position_jitter(), aes_string(fill=aes$fill), 
            shape=21L, size=4L) + 
        xlab(theme$xlab) + ylab(theme$ylab) + theme_bw() + 
        facet_grid(reformulate(aes$facet_rhs, aes$facet_lhs), scales="free")  + 
        (if (is.null(aes$x_breaks_discrete))
            scale_x_discrete(breaks = NULL)) + 
        theme(legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06), 
            legend.position="bottom") + 
        ggtitle(theme$title)
    plotly::ggplotly(p, tooltip = c("sample", "x", "tooltip1"))
}

qualimap_box <- function(dt, aes, theme, interactive=TRUE) {

    p = ggplot(dt, aes_string(x=aes$x, y=aes$y)) + 
        geom_boxplot(aes_string(fill=aes$fill)) + 
        ylab(theme$ylab) + theme_bw() + 
        facet_wrap(~ alignment, scales="free") + 
        scale_x_discrete(breaks=NULL) + 
        theme(legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06), 
            legend.position="bottom") + 
        ggtitle(theme$title)
    if (interactive) {
        if (!requireNamespace("plotly"))
            stop("Package 'plotly' is not available.")
        p = plotly::ggplotly(p)
    }
    p
}

qualimap_bar <- function(dt, aes, theme, interactive=TRUE) {

    p = ggplot(dt, aes_string(x=aes$x, y=aes$y, fill=aes$fill, 
                percent_tooltip=aes$percent_tooltip, group=aes$group)) + 
        geom_bar(stat="identity", position=aes$position) + 
        xlab(theme$xlab) + ylab(theme$ylab) + theme_bw() + 
        facet_wrap(~ .id, scales="free") + 
        theme(legend.text=element_text(size = 12L), 
            legend.background=element_blank(),
            legend.key=element_blank(),
            panel.grid.major.y=element_line(
                colour="#333333", size=0.12), 
            panel.grid.major.x=element_blank(), 
            panel.grid.minor=element_line(
                colour="#999999", size=0.06), 
            legend.position="bottom") + 
        ggtitle(theme$title)
    if (interactive) {
        if (!requireNamespace("plotly"))
            stop("Package 'plotly' is not available.")
        p = plotly::ggplotly(p, tooltip=c("group", "percent_tooltip"))
    }
    p
}
