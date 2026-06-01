# Pycnophylactic (Tobler 1979) disaggregation
#
# Iterative method: smooth the fine field with a Laplacian kernel, then
# restore each coarse block's mean by adding a per-block constant. Repeat
# until convergence. Produces a fine raster that is smooth in the Laplacian
# sense and whose block means equal the input coarse data exactly.
#
# References:
#   Tobler, W. R. (1979). Smooth pycnophylactic interpolation for
#     geographical regions. J. Am. Stat. Assoc. 74(367), 519-530.
#   Brunsdon's `pycno` R package implements the polygon-source version
#     using the same algorithmic core.

#' Pycnophylactic disaggregation (Tobler 1979)
#'
#' Iteratively smooths the fine raster while restoring each coarse block's
#' mean at every iteration. The result is a Laplacian-smooth surface whose
#' block means equal the input coarse data exactly.
#'
#' Compared to [disagg_bl()]:
#' \itemize{
#'   \item Tobler produces smoother output (no kinks at cell centers,
#'     where bilinear has them).
#'   \item Tobler is iterative and slower (~10-100x).
#'   \item Tobler tends to produce less ringing/overshoot near sharp
#'     coarse-scale gradients than the cubic prefilter method.
#' }
#'
#' @param coarse SpatRaster.
#' @param fact Integer disagg factor (>= 2).
#' @param max_iter Integer, maximum iterations. Default 100.
#' @param tol Numeric, convergence tolerance as a fraction of the coarse
#'   raster's value range. Iteration stops when the largest per-cell
#'   change between successive iterations falls below `tol * range(coarse)`.
#'   Default 1e-4.
#' @param smoother Character, smoothing kernel. One of:
#'   * `"laplacian_9"` (default): full 3x3 mean, including diagonals.
#'   * `"laplacian_5"`: 5-point stencil (cardinal neighbors + center).
#' @param variant Character, mean-restoration variant. One of:
#'   * `"additive"` (default): add a per-block constant. Works for any data.
#'   * `"multiplicative"`: scale within-block values. Strictly positive
#'     data only; will warn and fall back to additive if any block mean
#'     is non-positive.
#' @param na_fill Character, NA handling. Same semantics as [disagg_bl()];
#'   `"auto"` (default) uses `"fill"` when the coarse raster has NAs,
#'   `"reflect"` otherwise.
#' @param initial Character, choice of initial fine raster before iteration
#'   begins. One of:
#'   * `"near"` (default): nearest-neighbor disagg of `coarse`. Each fine
#'     cell starts at its parent coarse cell's value. Matches the original
#'     Tobler (1979) formulation.
#'   * `"bilinear"`: standard bilinear disagg of `coarse`. Smoother starting
#'     point that typically converges faster and produces less block-edge
#'     blockiness at high disagg factors. Final result still satisfies the
#'     same convergence criterion; only the iteration trajectory differs.
#' @param verbose Logical, print iteration progress. Default `FALSE`.
#' @return Fine SpatRaster.
#'
#' @references
#' Tobler, W. R. (1979). Smooth pycnophylactic interpolation for
#' geographical regions. *J. Am. Stat. Assoc.* 74(367), 519-530.
#'
#' @export
disagg_pyc <- function(coarse, fact,
                       max_iter = 100L,
                       tol = 1e-4,
                       smoother = c("laplacian_9", "laplacian_5"),
                       variant = c("additive", "multiplicative"),
                       initial = c("near", "bilinear"),
                       na_fill = c("auto", "reflect", "fill"),
                       verbose = FALSE) {
      smoother <- match.arg(smoother)
      variant  <- match.arg(variant)
      initial  <- match.arg(initial)
      na_fill  <- match.arg(na_fill)
      .validate_inputs_pyc(coarse, fact, max_iter, tol)
      fact <- as.integer(fact)
      max_iter <- as.integer(max_iter)

      if (na_fill == "auto") {
            na_fill <- if (any(is.na(terra::values(coarse)))) "fill" else "reflect"
      }

      # Handle NAs in coarse via padding if requested
      if (na_fill == "fill") {
            coarse_in <- .pad_with_nearest(coarse, fact)  # pad by one block
            crop_back <- TRUE
      } else {
            coarse_in <- coarse
            crop_back <- FALSE
      }

      # Multiplicative variant validation
      if (variant == "multiplicative") {
            coarse_vals <- terra::values(coarse_in)
            if (any(coarse_vals <= 0, na.rm = TRUE)) {
                  warning("Multiplicative variant requires strictly positive coarse ",
                          "values; some are <= 0. Falling back to additive.")
                  variant <- "additive"
            }
      }

      rng <- diff(range(terra::values(coarse_in), na.rm = TRUE))
      if (rng == 0) rng <- 1
      abs_tol <- tol * rng

      # Initial guess for the fine raster
      fine <- terra::disagg(coarse_in, fact = fact, method = initial)

      # Smoother kernel
      w <- switch(smoother,
                  "laplacian_9" = matrix(1/9, 3, 3),
                  "laplacian_5" = { k <- matrix(0, 3, 3)
                  k[2, ] <- 1/5; k[, 2] <- 1/5
                  k[2, 2] <- 1/5  # center counted once
                  k })

      prev_fine <- fine
      converged <- FALSE
      for (it in seq_len(max_iter)) {
            # Smoothing step
            fine <- terra::focal(fine, w = w, fun = "sum",
                                 na.policy = "omit", expand = TRUE)

            # Block-mean restoration
            block_mean <- terra::aggregate(fine, fact = fact, fun = "mean")
            if (variant == "additive") {
                  correction_coarse <- coarse_in - block_mean
                  correction_fine <- terra::disagg(correction_coarse, fact = fact,
                                                   method = "near")
                  fine <- fine + correction_fine
            } else {  # multiplicative
                  ratio_coarse <- coarse_in / block_mean
                  ratio_fine <- terra::disagg(ratio_coarse, fact = fact, method = "near")
                  fine <- fine * ratio_fine
            }

            # Convergence check
            delta <- max(abs(terra::values(fine) - terra::values(prev_fine)),
                         na.rm = TRUE)
            if (verbose) {
                  cat(sprintf("  pyc iter %3d:  max|delta| = %.3e  (tol = %.3e)\n",
                              it, delta, abs_tol))
            }
            if (delta < abs_tol) {
                  converged <- TRUE
                  break
            }
            prev_fine <- fine
      }

      if (!converged) {
            warning("Pycnophylactic iteration did not converge in ", max_iter,
                    " iterations (last delta = ", signif(delta, 3),
                    ", tolerance = ", signif(abs_tol, 3),
                    "). Result may not perfectly preserve block means.")
      }

      # Crop back to original extent if we padded
      if (crop_back) {
            fine_target_ext <- terra::ext(terra::disagg(coarse, fact = fact,
                                                        method = "near"))
            fine <- terra::crop(fine, fine_target_ext)
      }

      names(fine) <- names(coarse)
      attr(fine, "iterations") <- it
      attr(fine, "converged") <- converged
      fine
}

.validate_inputs_pyc <- function(coarse, fact, max_iter, tol) {
      if (!inherits(coarse, "SpatRaster"))
            stop("`coarse` must be a SpatRaster")
      if (length(fact) != 1L || fact < 2L || fact != as.integer(fact))
            stop("`fact` must be an integer >= 2")
      if (length(max_iter) != 1L || max_iter < 1L || max_iter != as.integer(max_iter))
            stop("`max_iter` must be a positive integer")
      if (!is.numeric(tol) || length(tol) != 1L || tol <= 0)
            stop("`tol` must be a positive numeric scalar")
}
