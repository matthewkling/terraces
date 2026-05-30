# Method registry
#
# Each disaggregation method is registered with a function that performs
# the disaggregation. The registry supports two method types:
#
#   "prefilter": method has a finite-support inverse kernel computed
#       analytically from a round-trip kernel. v0.1 includes bilinear;
#       cubic planned for v0.2.
#
#   "iterative": method runs an iterative algorithm directly on the
#       fine grid. Planned for v0.2 (pycnophylactic / Tobler).
#
# The registry abstraction is method-type-agnostic: a method just needs
# a `disagg_fn(coarse, fact, ...)` that returns a fine SpatRaster.
# Method-specific parameters (radius, max_iter, etc.) are passed via `...`
# and the method's `disagg_fn` validates them.

.ces_methods <- new.env(parent = emptyenv())

#' Register a disaggregation method with `terraces`
#'
#' @param name Character, method name.
#' @param disagg_fn Function `function(coarse, fact, ...)` returning a
#'   fine SpatRaster.
#' @param type Character, one of "prefilter" or "iterative". Determines
#'   what auxiliary functions are expected; prefilter methods also need
#'   `roundtrip_fn` and `default_radius`.
#' @param roundtrip_fn (prefilter methods only) Function returning the
#'   round-trip kernel as a square odd-sized matrix.
#' @param default_radius (prefilter methods only) Function returning the
#'   default inverse-kernel radius given a disagg factor.
#' @param description Optional description shown by `ces_list_methods()`.
#'
#' @return Invisibly the method name.
#' @export
ces_register_method <- function(name, disagg_fn,
                                type = c("prefilter", "iterative"),
                                roundtrip_fn = NULL,
                                default_radius = NULL,
                                description = "") {
      type <- match.arg(type)
      stopifnot(is.character(name), length(name) == 1L)
      stopifnot(is.function(disagg_fn))

      if (type == "prefilter") {
            if (is.null(roundtrip_fn) || is.null(default_radius)) {
                  stop("Prefilter methods require both `roundtrip_fn` and ",
                       "`default_radius`.")
            }
            # quick sanity-check
            K_test <- roundtrip_fn(2L)
            if (!is.matrix(K_test) || nrow(K_test) != ncol(K_test) ||
                nrow(K_test) %% 2L != 1L) {
                  stop("`roundtrip_fn` must return a square matrix with odd side length")
            }
            if (abs(sum(K_test) - 1) > 1e-10) {
                  stop("Round-trip kernel must sum to 1; got ", sum(K_test))
            }
      }

      .ces_methods[[name]] <- list(
            disagg_fn      = disagg_fn,
            type           = type,
            roundtrip      = roundtrip_fn,
            default_radius = default_radius,
            description    = description
      )
      invisible(name)
}

#' List registered disaggregation methods
#' @return Data frame with columns name, type, description.
#' @export
ces_list_methods <- function() {
      names <- ls(.ces_methods)
      if (length(names) == 0L) {
            return(data.frame(name = character(), type = character(),
                              description = character(), stringsAsFactors = FALSE))
      }
      data.frame(
            name = names,
            type = vapply(names, function(n) .ces_methods[[n]]$type, character(1)),
            description = vapply(names, function(n) .ces_methods[[n]]$description,
                                 character(1)),
            stringsAsFactors = FALSE
      )
}

# Internal accessor with friendly error
.get_method <- function(name) {
      if (!exists(name, envir = .ces_methods, inherits = FALSE)) {
            avail <- paste(ls(.ces_methods), collapse = ", ")
            stop("Unknown method '", name, "'. Available methods: ", avail)
      }
      get(name, envir = .ces_methods, inherits = FALSE)
}
