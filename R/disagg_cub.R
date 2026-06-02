# Pre-sharpened cubic disaggregation
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

#' Pre-sharpened cubic disaggregation
#'
#' Disaggregates `coarse` to a finer resolution by integer factor `fact`
#' using Keys cubic convolution interpolation, with a pre-sharpening filter
#' that makes the result mass-preserving: aggregating it back by block
#' mean recovers the input coarse values up to a small truncation error
#' governed by `radius`.
#'
#' Compared to [disagg_bl()]:
#' \itemize{
#'   \item Cubic produces smoother output (no kinks at cell centers, where
#'     bilinear has them) at modestly higher computational cost.
#'   \item Cubic has negative side lobes and may overshoot input range
#'     near sharp gradients. For data with hard physical bounds
#'     (precipitation, fractional cover), this may matter.
#'   \item At the same `radius`, cubic's truncation error is somewhat
#'     larger than bilinear's because its inverse kernel decays more
#'     slowly. The default `radius` is correspondingly larger to
#'     compensate.
#' }
#'
#' @param coarse SpatRaster. Multi-layer rasters are supported.
#' @param fact Integer disagg factor (>= 2).
#' @param radius Integer half-width of the inverse kernel. If `NULL`
#'   (the default), `radius` is set to `max(9, fact + 4)`, slightly
#'   larger than bilinear's default because cubic's inverse kernel
#'   decays more slowly. Users rarely need to change this. Larger
#'   values reduce truncation error further at a small one-time cost
#'   (the kernel is cached per `fact`); smaller values are faster on
#'   small rasters but allow larger residual error. Automatically
#'   reduced for small rasters; see `max_radius_frac`.
#' @param max_radius_frac Numeric. Same semantics as in [disagg_bl()];
#'   upper bound on kernel radius as a fraction of the coarse raster's
#'   smaller dimension.
#' @param na_fill Character. Boundary/NA handling mode, same options as in
#'   [disagg_bl()]: `"auto"` (default), `"reflect"`, or `"fill"`.
#'
#' @return Fine SpatRaster.
#'
#' @details
#' **Mass preservation is approximate.** Like [disagg_bl()], the inverse
#' kernel is a finite approximation to an ideal infinite operator. The
#' truncation error is concentrated near the raster boundary and decays
#' with larger `radius`. Cubic's error at equal radius is somewhat larger
#' than bilinear's because its inverse kernel has alternating-sign rings
#' that decay more slowly than bilinear's geometric decay. See
#' [edge_effects()] to visualize which fine cells fall in the
#' edge-affected zone for given inputs.
#'
#' **Implementation note.** [terra::disagg()] does not expose cubic
#' interpolation directly, so `disagg_cub()` uses [terra::resample()]
#' internally for the interpolation step, with a target fine grid that
#' nests perfectly within the coarse cells. We have verified empirically
#' that `terra::resample(method = "cubic")` uses Keys cubic with
#' `a = -0.5`, which is the convention our pre-sharpening kernel inverts.
#'
#' @references
#' Keys, R. G. (1981). Cubic convolution interpolation for digital
#' image processing. *IEEE Trans. Acoust. Speech Signal Process.*
#' 29(6), 1153–1160.
#'
#' @examples
#' library(terra)
#' coarse <- rast(nrows = 30, ncols = 30, vals = runif(900))
#' fine_cub <- disagg_cub(coarse, fact = 5)
#' back     <- aggregate(fine_cub, fact = 5, fun = "mean")
#' max(abs(values(coarse) - values(back)))    # small; reduced by larger radius
#'
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
