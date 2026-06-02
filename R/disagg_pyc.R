# Pycnophylactic (Tobler 1979) disaggregation
#
# Iterative method: smooth the fine field with a Laplacian kernel, then
# restore each coarse block's mean by adding a per-block constant. Repeat
# until convergence. Produces a fine raster that is smooth in the Laplacian
# sense and whose block means equal the input coarse data exactly (by
# construction — mass preservation is enforced at every iteration).
#
# References:
#   Tobler, W. R. (1979). Smooth pycnophylactic interpolation for
#     geographical regions. J. Am. Stat. Assoc. 74(367), 519-530.
#   The `pycno` R package implements the polygon-source version
#     using the same algorithmic core.

#' Pycnophylactic disaggregation (Tobler 1979)
#'
#' Iteratively smooths the fine raster while restoring each coarse block's
#' mean at every iteration. The result is a Laplacian-smooth surface whose
#' block means equal the input coarse values exactly — mass preservation
#' is enforced at every iteration by construction, not approximated.
#'
#' Compared to the pre-sharpening methods ([disagg_bl()], [disagg_cub()]):
#' \itemize{
#'   \item Pycnophylactic produces smoother output (no kinks at cell
#'     centers, where bilinear has them).
#'   \item Pycnophylactic is *exactly* mass-preserving by construction;
#'     the pre-sharpening methods are approximately so, with a small
#'     truncation error concentrated near the raster boundary. For
#'     small rasters or boundary-sensitive applications, this matters.
#'   \item Pycnophylactic tends to produce less ringing near sharp
#'     coarse-scale gradients than [disagg_cub()].
#'   \item Pycnophylactic is iterative and substantially slower
#'     (~10-50x the cost of [disagg_bl()]).
#' }
#'
#' **Caveat at large disaggregation factors.** At factors much larger
#' than the smoother's bandwidth (say, `fact` > ~20 with the default 3x3
#' kernel), the iteration's fixed point can show subtle within-block
#' flattening: high-frequency content within a coarse block becomes
#' piecewise-constant at block scale, limiting how much curvature the
#' method can express within a block. This is a structural property of
#' the iterative algorithm on regular grids, not a convergence issue —
#' running more iterations doesn't remove it. Consider [disagg_cub()]
#' for high-factor use cases. Using `init = "bilinear"` can also reduce
#' the effect somewhat.
#'
#' @param coarse SpatRaster.
#' @param fact Integer disagg factor (>= 2).
#' @param max_iter Integer, maximum iterations. Default 100.
#' @param tol Numeric, convergence tolerance as a fraction of the coarse
#'   raster's value range. Iteration stops when the largest per-cell
#'   change between successive iterations falls below `tol * range(coarse)`.
#'   Default 1e-4. Controls only how close the result is to its smooth
#'   fixed point; mass preservation is enforced exactly at every iteration
#'   regardless of convergence.
#' @param smoother Character, smoothing kernel. One of:
#'   * `"laplacian_9"` (default): full 3x3 mean, including diagonals.
#'   * `"laplacian_5"`: 5-point stencil (cardinal neighbors + center).
#' @param variant Character, mean-restoration variant. One of:
#'   * `"additive"` (default): add a per-block constant. Works for any data.
#'   * `"multiplicative"`: scale within-block values. Strictly positive
#'     data only; will warn and fall back to additive if any block mean
#'     is non-positive.
#' @param init Character, initial fine raster before iteration begins. One of:
#'   * `"near"` (default): nearest-neighbor disagg of `coarse`. Each fine
#'     cell starts at its parent coarse cell's value. Matches the original
#'     Tobler (1979) formulation.
#'   * `"bilinear"`: standard bilinear disagg of `coarse`. A smoother
#'     starting point that typically converges faster and produces less
#'     within-block flattening at high disaggregation factors. The final
#'     mass-preservation property is unchanged; only the iteration
#'     trajectory and fixed point differ.
#' @param na_fill Character, NA handling. Same semantics as [disagg_bl()];
#'   `"auto"` (default) uses `"fill"` when the coarse raster has NAs,
#'   `"reflect"` otherwise.
#' @param verbose Logical, print iteration progress. Default `FALSE`.
#'
#' @return Fine SpatRaster, with attributes `iterations` (number actually
#'   used) and `converged` (logical).
#'
#' @references
#' Tobler, W. R. (1979). Smooth pycnophylactic interpolation for
#' geographical regions. *J. Am. Stat. Assoc.* 74(367), 519-530.
#'
#' @seealso [disagg_bl()] and [disagg_cub()] for faster, non-iterative
#'   pre-sharpening alternatives. The `pycno` package
#'   (`pycno::pycno`) implements the polygon-source version of the
#'   underlying algorithm.
#'
#' @examples
#' library(terra)
#' coarse <- rast(matrix(runif(900), 30, 30))
#' fine_pyc <- disagg_pyc(coarse, fact = 5)
#' back     <- aggregate(fine_pyc, fact = 5, fun = "mean")
#' max(abs(values(coarse) - values(back)))  # ~ machine precision
#' @export
disagg_pyc <- function(coarse, fact,
                       max_iter = 100L,
                       tol = 1e-4,
                       smoother = c("laplacian_9", "laplacian_5"),
                       variant = c("additive", "multiplicative"),
                       init = c("near", "bilinear"),
                       na_fill = c("auto", "reflect", "fill"),
                       verbose = FALSE) {
      smoother <- match.arg(smoother)
      variant  <- match.arg(variant)
      init  <- match.arg(init)
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

      # Suppress terra progress bars during iteration (many small ops would
      # otherwise produce flickering empty bars). Restored on function exit.
      old_progress <- terra::terraOptions(print = FALSE)$progress
      terra::terraOptions(progress = 0)
      on.exit(terra::terraOptions(progress = old_progress), add = TRUE)

      # Initial guess for the fine raster
      fine <- terra::disagg(coarse_in, fact = fact, method = init)

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
                    "). Block means are still preserved exactly; ",
                    "the result may simply be less smooth than requested.")
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
