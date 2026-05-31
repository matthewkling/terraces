# Keys cubic disaggregation via prefilter
#
# Cubic version of disagg_bl(). The math is identical in structure:
# precompute an inverse kernel that undoes the 5x5 round-trip operator,
# apply it as a focal pass on the coarse raster, then do standard cubic
# interpolation.
#
# Note: terra::disagg() doesn't expose cubic, so we use terra::resample()
# instead, with the target fine grid built to nest perfectly within the
# coarse cells. Verified empirically that terra::resample(method = "cubic")
# uses Keys with a = -0.5, which is the convention our kernel inverts.

#' Cubic disaggregation with conservation of empirical statistics
#'
#' Disaggregates `coarse` to a finer resolution by integer factor `fact`
#' using Keys cubic convolution interpolation, with a prefilter that makes
#' the result aggregation-consistent: aggregating it back by block mean
#' recovers `coarse` exactly (up to a small truncation error governed by
#' `radius`).
#'
#' Compared to [disagg_bl()]:
#' \itemize{
#'   \item Cubic produces smoother output (C¹ continuity at coarse cell
#'     boundaries) at modestly higher computational cost.
#'   \item Cubic can place within-block extrema in cell interiors (bilinear
#'     can only place them at coarse cell vertices).
#'   \item Cubic has negative side lobes and may overshoot input range
#'     near sharp gradients. For data with hard physical bounds
#'     (precipitation, fractional cover), this may matter.
#' }
#'
#' @param coarse SpatRaster. Multi-layer rasters are supported.
#' @param fact Integer disagg factor (>= 2).
#' @param radius Integer half-width of the inverse kernel. `NULL` uses a
#'   method-specific default (slightly larger than bilinear's because
#'   cubic's round-trip operator is more diffusive).
#' @param max_radius_frac Numeric. Same semantics as in [disagg_bl()];
#'   upper bound on kernel radius as a fraction of the coarse raster's
#'   smaller dimension.
#' @param na_fill Character. Boundary/NA handling mode, same options as in
#'   [disagg_bl()]: `"auto"` (default), `"reflect"`, or `"fill"`.
#'
#' @return Fine SpatRaster.
#' @export
disagg_cub <- function(coarse, fact, radius = NULL,
                       max_radius_frac = 1/3,
                       na_fill = c("auto", "reflect", "fill")) {
      na_fill <- match.arg(na_fill)
      .validate_inputs(coarse, fact, max_radius_frac)
      fact <- as.integer(fact)

      if (is.null(radius)) {
            radius <- .default_radius("cubic", fact)
      }
      radius <- as.integer(radius)

      min_dim <- min(nrow(coarse), ncol(coarse))
      r_hard <- (min_dim - 1L) %/% 2L
      r_soft <- floor(max_radius_frac * min_dim)
      r_max  <- min(r_hard, r_soft)
      if (r_max < 3L) {
            stop("Coarse raster too small (", nrow(coarse), " x ", ncol(coarse),
                 ") for cubic; need at least ~9 cells per side.")
      }
      if (radius > r_max) {
            warning("Inverse-kernel radius reduced from ", radius, " to ", r_max,
                    " to fit coarse raster (", nrow(coarse), " x ", ncol(coarse),
                    "); round-trip accuracy will be lower than requested.")
            radius <- r_max
      }

      K_inv <- kernel("cubic", fact, radius)
      apply_kernel(coarse, K_inv, fact, method = "cubic", na_fill = na_fill)
}
