
## Style definitions for the graphs
STYLE_COLOR <- c('EN' = 'green',
                 'Ada. EN' = 'green',
                 'I-LAMM' = 'gray',
                 'I-LAMM (SCAD)' = 'gray',
                 'MM' = 'lightblue',
                 'Ada. MM' = 'lightblue',
                 'PENSE' = 'blue',
                 'Ada. PENSE' = 'blue')

STYLE_SHAPE <- c('EN' = 0,
                 'Ada. EN' = 15,
                 'I-LAMM' = 2,
                 'I-LAMM (SCAD)' = 17,
                 'MM' = 5,
                 'Ada. MM' = 23,
                 'PENSE' = 1,
                 'Ada. PENSE' = 19)

STYLE_LINETYPE <- c('EN' = '22',
                    'Ada. EN' = 'solid',
                    'I-LAMM' = '22',
                    'I-LAMM (SCAD)' = 'solid',
                    'MM' = '22',
                    'Ada. MM' = 'solid',
                    'PENSE' = '22',
                    'Ada. PENSE' = 'solid')

#' Recode the contamination setting to a human readable string
#'
#' @param x contamination proportion
#' @return factor with human readable strings in the correct order
recode_cont <- function (x) {
  cont_str <- factor(if_else(x < .Machine$double.eps,
                             "No contamination",
                             sprintf("%.0f%% contamination", 100 * x)))
  fct_reorder(cont_str, x)
}

#' Recode the estimatino method to a human readable string
#'
#' @param x estimation method identifier
#' @return factor with human readable strings in the correct order
recode_method <- function (x) {
  recode_factor(x,
                adapense = 'Ada. PENSE',
                adamm = 'Ada. MM',
                adaen = 'Ada. EN',
                ilammscad = 'I-LAMM (SCAD)',
                ilamm = 'I-LAMM',
                pense = 'PENSE',
                mm = 'MM',
                en = 'EN',
                .default = NA_character_)
}

#' Recode the residual distribution to a human readable string
#'
#' @param x residual distribution identifier
#' @return factor with human readable strings in the correct order
recode_resid_dist <- function (x) {
  recode_factor(x,
                `norm()` = 'Normal residuals',
                `stable(alpha = 1.33, beta = 0)` = 'Stable(1.33) residuals',
                `cauchy()` = 'Cauchy residuals')
}

#' Colorblind-Friendly Color Palette
#'
#' A list of 9 colors which are somewhat easy to distinguish even with mild forms of common issues
#' in color-perception.
#'
#' Although the colors are "friendly" for many different types of
#' color-blindness, it's not guaranteed that the contrast between the colors
#' is great for anyone!
#' Therefore, these colors should not be the only differentiating factor between elements
#' in a graph.
#'
#' The following colors are available:
#'
#' * `lightorange`
#' * `lightblue`
#' * `green`
#' * `yellow`
#' * `blue`
#' * `orange`
#' * `cyan`
#' * `gray`
#' * `lightgray`
#'
#' @param ... colors to choose (possibly named) (see examples).
#'   If none are given, all colors are returned.
#' @return a named vector with color codes
#'
#' @examples
#' colorblind_cols(air = 'blue', fruit = 'orange')
colorblind_cols <- function (...) {
  colors <- c(lightorange = '#E69F00',
              lightblue = '#56B4E9',
              green = '#009E73',
              yellow = '#F0E442',
              blue = '#0072B2',
              orange = '#D55E00',
              cyan = '#CC79A7',
              gray = '#333333',
              lightgray = '#CCCCCC')
  requested <- as.character(...)
  if (length(requested) == 0L) {
    colors
  } else {
    colors <- colors[requested]
    colors[!is.na(colors)]
  }
}

#' Create a Colorblind-Friendly Color Palette
#'
#' Similar to [colorblind_pal()], but returns a function
#'
#' @param values a named character vector of chosen colors.
#' @return a function which takes the number of colors and returns a vector
#'   of that length with the chosen colors.
colorblind_pal <- function (values) {
  palette <- colorblind_cols()
  if (is.null(values)) {
    values <- seq_along(palette)
  }

  color_specs <- grepl('#[a-fA-F0-9]{6}', values)

  pal <- unname(values)
  pal[!color_specs] <- unname(palette[values[!color_specs]])

  if (anyNA(pal)) {
    stop(sprintf("Invalid color(s) selected Available colors are: %s",
                 paste(names(palette), collapse = ', ')))
  }

  if (!is.null(names(values))) {
    names(pal) <- names(values)
  }

  function (n) {
    if (n > length(pal)) {
      stop(sprintf("Requested more colors (%d) than available (%d)", n, length(pal)))
    }
    pal
  }
}

#' Colorblind-Friendly Color Palette to use with ggplot2
#' See [colorblind_pal()] for available colors.
scale_color_colorblind <- function(..., values = NULL) {
  if (!require('ggplot2', quietly = TRUE)) {
    stop("ggplot2 package not available")
  }
  discrete_scale('colour', 'colorblind', colorblind_pal(values), ...)
}

#' Colorblind-Friendly Color Palette to use with ggplot2
#' See [colorblind_pal()] for available colors.
scale_fill_colorblind <- function(..., values = NULL) {
  if (!require('ggplot2', quietly = TRUE)) {
    stop("ggplot2 package not available")
  }

  discrete_scale('fill', 'colorblind', colorblind_pal(values), ...)
}

#' ggplot2 Theme Used For Plots In The Manuscript
#'
#' Theme based on the [ggplot2][theme_bw()] theme.
#'
#' @param base_family,base_size see [ggplot2][theme_bw()] for the meaning of these
#'   parameters.
ggplot_theme <- function(base_size = 12, base_family = '') {
  if (!require('ggplot2', quietly = TRUE)) {
    stop("Package `ggplot2` not available")
  }

  theme_bw(base_size = base_size, base_family = base_family) +
    theme(plot.title = element_text(size = rel(0.8), hjust = 0.5, face = 'bold'),
          plot.margin = margin(0.2, 0.4, 0, 0.2, 'lines'),
          panel.background = element_rect(fill = 'transparent', color = NA),
          plot.background = element_rect(fill = 'transparent', color = NA),
          legend.title = element_text(size = rel(0.75)),
          legend.text = element_text(size = rel(0.75)),
          legend.margin = margin(),
          axis.title = element_text(size = rel(0.75)),
          axis.text = element_text(size = rel(0.7)),
          axis.title.y = element_text(vjust = 1),
          axis.title.x = element_text(vjust = 0),
          panel.grid.major = element_line(color = 'gray30', size = rel(0.5), linetype='dotted'),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = '#ffffff', color = 'gray50', size = 0.3),
          strip.text = element_text(size = rel(0.75)),
          panel.border = element_rect(color = 'gray50', size = 0.3),
          legend.background = element_blank())
}

