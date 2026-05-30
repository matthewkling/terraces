# Disaggregation functions
#
# Two-tier API:
#   ces_disagg() - top-level dispatcher; pick method by name
#   ces_disagg_bl() - method-specific entry point (bilinear)
#   (planned: ces_disagg_cub, ces_disagg_pyc)

#' Bilinear disaggregation with conservation of empirical statistics
#'
#' Disaggregates `coarse` to a finer resolution by integer factor `fact`,
#' such that aggregating the result back by block mean recovers `coarse`
#' exactly (up to a small truncation error governed by `radius`). The fine
#' raster is a true bilinear surface — same smoothness as standard bilinear
#' disagg, just debiased.
#'
#' Internally:
#' 1. Computes (or looks up cached) the bilinear inverse kernel for `fact`.
#' 2. Applies the kernel as a focal pass on `coarse`, yielding an adjusted
#'    coarse raster whose extremes are amplified to compensate for bilinear's
#'    regression-to-the-mean bias.
#' 3. Bilinear-disaggregates the adjusted coarse raster.
#'
#' @param coarse SpatRaster. Multi-layer rasters are supported (the kernel
#'   is applied to each layer in turn).
#' @param fact Integer disagg factor (>= 2).
#' @param radius Integer half-width of the inverse kernel. If `NULL`, a
#'   sensible method-specific default is used (typically 7-11). May be
#'   automatically reduced for small rasters; see `max_radius_frac`.
#' @param max_radius_frac Numeric. Upper bound on radius as a fraction of
#'   the coarse raster's smaller dimension. Defaults to `1/3` to keep the
#'   focal-boundary band (where reflective extension is approximate) from
#'   dominating the interior. Allowed range: (0, 0.5].
#'
#' @return Fine SpatRaster. Layer names preserved.
#'
#' @details
#' **Boundary handling.** The focal pass uses reflective extension at raster
#' edges. Block-mean preservation is slightly approximate in a boundary band
#' of width `radius`; for continental rasters this is invisible. Future
#' versions may add an "exact" boundary mode via a small global solve.
#'
#' **Numerical accuracy.** Round-trip RMSE is controlled by the `tail_max`
#' attribute of the cached kernel (see [ces_kernel()]). Default radii target
#' ~1e-5 to ~1e-7 relative round-trip error.
#'
#' @references
#' Tobler, W. R. (1979). Smooth pycnophylactic interpolation for geographical
#' regions. *J. Am. Stat. Assoc.* 74(367), 519-530.
#'
#' Unser, M. (1999). Splines: A perfect fit for signal and image processing.
#' *IEEE Signal Processing Magazine* 16(6), 22-38.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' coarse <- rast(nrows = 30, ncols = 30, vals = runif(900))
#'
#' # standard bilinear: biased
#' fine_std <- disagg(coarse, fact = 5, method = "bilinear")
#' back_std <- aggregate(fine_std, fact = 5, fun = "mean")
#' max(abs(values(coarse) - values(back_std)))   # nonzero
#'
#' # conservation-enforcing bilinear: empirical stats preserved
#' fine_ces <- ces_disagg_bl(coarse, fact = 5)
#' back_ces <- aggregate(fine_ces, fact = 5, fun = "mean")
#' max(abs(values(coarse) - values(back_ces)))   # ~ machine epsilon
#' }
#' @export
#' @param na_fill Boundary/NA handling mode passed to [ces_apply_kernel()].
#'   `"auto"` (default) selects `"reflect"` for rasters with no NAs and
#'   `"fill"` for rasters with NAs.
ces_disagg_bl <- function(coarse, fact, radius = NULL,
                          max_radius_frac = 1/3,
                          na_fill = c("auto", "reflect", "fill")) {
      na_fill <- match.arg(na_fill)
      .validate_inputs(coarse, fact, max_radius_frac)
      fact <- as.integer(fact)

      if (is.null(radius)) {
            radius <- .get_method("bilinear")$default_radius(fact)
      }
      radius <- as.integer(radius)

      min_dim <- min(nrow(coarse), ncol(coarse))
      r_hard <- (min_dim - 1L) %/% 2L
      r_soft <- floor(max_radius_frac * min_dim)
      r_max  <- min(r_hard, r_soft)
      if (r_max < 2L) {
            stop("Coarse raster too small (", nrow(coarse), " x ", ncol(coarse),
                 ") for this method; need at least ~6 cells per side.")
      }
      if (radius > r_max) {
            warning("Inverse-kernel radius reduced from ", radius, " to ", r_max,
                    " to fit coarse raster (", nrow(coarse), " x ", ncol(coarse),
                    "); round-trip accuracy will be lower than requested.")
            radius <- r_max
      }

      K_inv <- ces_kernel("bilinear", fact, radius)
      ces_apply_kernel(coarse, K_inv, fact, method = "bilinear", na_fill = na_fill)
}

#' Disaggregate a raster with conservation of empirical statistics
#'
#' Top-level dispatcher: pick a method by name. Method-specific arguments
#' pass through via `...`. For most uses, calling the method-specific
#' function (e.g. [ces_disagg_bl()]) directly gives clearer signatures
#' and per-method documentation.
#'
#' @param coarse SpatRaster.
#' @param fact Integer disagg factor.
#' @param method Character, one of the registered methods. See
#'   [ces_list_methods()]. Default `"bilinear"`.
#' @param ... Method-specific arguments forwarded to the chosen method's
#'   disagg function.
#'
#' @return Fine SpatRaster.
#' @export
ces_disagg <- function(coarse, fact, method = "bilinear", ...) {
      m <- .get_method(method)
      m$disagg_fn(coarse, fact, ...)
}

#' Apply a precomputed inverse kernel to perform CES disaggregation
#'
#' Power-user function: skip the kernel-lookup step when disaggregating
#' many rasters with the same geometry and method.
#'
#' @param coarse SpatRaster.
#' @param K_inv Inverse kernel from [ces_kernel()].
#' @param fact Integer disagg factor; should match the kernel.
#' @param method Character, the prefilter method the kernel was built for.
#' @param na_fill Character, boundary/NA handling mode. One of:
#'   * `"auto"` (default): use `"reflect"` if `coarse` has no NA values,
#'     otherwise use `"fill"`.
#'   * `"reflect"`: use terra's reflective extension at boundaries.
#'     Fastest and most accurate for round-trip preservation, but
#'     propagates NAs (any NA in the coarse raster yields an
#'     `(2*radius+1)`-cell wide NA region in the fine output).
#'   * `"fill"`: pad the coarse raster with nearest-valid values before
#'     the focal pass, then crop back. Robust to NAs, but slightly
#'     degrades round-trip preservation in the boundary band.
#' @return Fine SpatRaster.
#' @export
ces_apply_kernel <- function(coarse, K_inv, fact, method,
                             na_fill = c("auto", "reflect", "fill")) {
      na_fill <- match.arg(na_fill)
      ka <- attributes(K_inv)
      if (!is.null(ka$fact)   && ka$fact   != as.integer(fact))
            warning("Kernel was built for fact = ", ka$fact, " but called with ", fact)
      if (!is.null(ka$method) && ka$method != method)
            warning("Kernel was built for method = '", ka$method, "' but called ",
                    "with '", method, "'")

      if (na_fill == "auto") {
            na_fill <- if (any(is.na(terra::values(coarse)))) "fill" else "reflect"
      }

      if (na_fill == "fill") {
            radius <- (nrow(K_inv) - 1L) %/% 2L
            coarse_to_focal <- .pad_with_nearest(coarse, radius)
      } else {
            coarse_to_focal <- coarse
      }

      coarse_adj <- terra::focal(coarse_to_focal, w = K_inv, fun = "sum",
                                 na.policy = "omit", expand = TRUE)

      # crop back to original coarse extent (no-op if no padding was applied)
      coarse_adj <- terra::crop(coarse_adj, terra::ext(coarse))

      fine <- terra::disagg(coarse_adj, fact = fact, method = "bilinear")
      names(fine) <- names(coarse)
      fine
}

# Pad a raster by `pad` cells on all sides, filling padded cells (and any
# pre-existing NAs in the original) by repeated nearest-valid-cell
# propagation. This makes the focal pass behave well at edges and around
# interior NAs.
.pad_with_nearest <- function(r, pad) {
      if (pad < 1L) return(r)

      res <- terra::res(r)
      e <- terra::ext(r)
      e_ext <- terra::ext(
            e[1] - pad * res[1],
            e[2] + pad * res[1],
            e[3] - pad * res[2],
            e[4] + pad * res[2]
      )
      r_ext <- terra::extend(r, e_ext)

      # Iteratively fill NAs by propagating values inward from valid cells.
      # Each pass: for any NA cell with at least one valid neighbor in a 3x3
      # window, replace it with the mean of its valid neighbors.
      #
      # We loop until no NA cells remain or progress stalls. The number of
      # iterations needed is roughly the maximum distance from any NA cell
      # to the nearest valid cell in the original raster.
      w <- matrix(1, 3, 3)
      prev_n_na <- Inf
      for (i in seq_len(max(pad + 5L, 50L))) {
            vals <- terra::values(r_ext)
            n_na <- sum(is.na(vals))
            if (n_na == 0L) break
            if (n_na >= prev_n_na) break  # stalled, no further progress
            prev_n_na <- n_na

            filled <- terra::focal(r_ext, w = w, fun = "mean",
                                   na.policy = "only", na.rm = TRUE)
            # focal can produce NaN for cells with all-NA neighborhoods; convert
            # back to NA so the next iteration can attempt them once their
            # neighbors are filled.
            f_vals <- terra::values(filled)
            f_vals[is.nan(f_vals)] <- NA
            terra::values(filled) <- f_vals

            r_ext <- filled
      }
      r_ext
}

.validate_inputs <- function(coarse, fact, max_radius_frac) {
      if (!inherits(coarse, "SpatRaster"))
            stop("`coarse` must be a SpatRaster")
      if (length(fact) != 1L || fact < 2L || fact != as.integer(fact))
            stop("`fact` must be an integer >= 2")
      if (!is.numeric(max_radius_frac) || max_radius_frac <= 0 ||
          max_radius_frac > 0.5)
            stop("`max_radius_frac` must be in (0, 0.5]")
}

# Register the built-in methods on package load
.register_builtins <- function() {
      ces_register_method(
            name           = "bilinear",
            type           = "prefilter",
            disagg_fn      = ces_disagg_bl,
            roundtrip_fn   = .roundtrip_bilinear,
            default_radius = .default_radius_bilinear,
            description    = "Bilinear with conservation of empirical statistics"
      )
      ces_register_method(
            name        = "pycnophylactic",
            type        = "iterative",
            disagg_fn   = ces_disagg_pyc,
            description = "Tobler's pycnophylactic interpolation (iterative)"
      )
}

.onLoad <- function(libname, pkgname) {
      .register_builtins()
}
