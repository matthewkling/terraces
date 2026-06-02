#' Apply a precomputed inverse kernel to perform aggregation-consistent disaggregation
#'
#' Power-user function: skip the kernel-lookup step when disaggregating
#' many rasters with the same geometry and method.
#'
#' @param coarse SpatRaster.
#' @param K_inv Inverse kernel from [kernel()].
#' @param fact Integer disagg factor; should match the kernel.
#' @param method Character, the prefilter method the kernel was built for.
#'   Determines the interpolation step: `"bilinear"` uses [terra::disagg()];
#'   `"cubic"` uses [terra::resample()] with `method = "cubic"`.
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
#' @keywords internal
#' @noRd
apply_kernel <- function(coarse, K_inv, fact, method,
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

      # Dispatch interpolation step on method
      fine <- switch(method,
                     "bilinear" = terra::disagg(coarse_adj, fact = fact, method = "bilinear"),
                     "cubic"    = {
                           # terra::disagg doesn't expose cubic, so use resample to a
                           # nested-target fine grid
                           fine_template <- terra::disagg(coarse_adj, fact = fact,
                                                          method = "near")
                           terra::resample(coarse_adj, fine_template, method = "cubic")
                     },
                     stop("apply_kernel: unsupported method '", method, "'")
      )

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
