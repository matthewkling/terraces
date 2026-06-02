
#' Pycnophylactic disaggregation (Tobler 1979)
#'
#' Iteratively smooths the fine raster while restoring each coarse block's
#' mean at every iteration. The result is a smooth surface whose block
#' means equal the input coarse values exactly — mass preservation is
#' enforced at every iteration by construction, not approximated.
#'
#' For composite `fact` values, the disaggregation is cascaded through
#' a sequence of smaller stages by default, which is substantially
#' faster than running at the full factor in one shot and produces
#' smoother results.
#'
#' Compared to the pre-sharpening methods ([disagg_bl()], [disagg_cub()]):
#' \itemize{
#'   \item Pycnophylactic and cubic produce smoother output than bilinear
#'     (no kinks at cell centers).
#'   \item Pycnophylactic is *exactly* mass-preserving by construction;
#'     the pre-sharpening methods are approximately so, with a small
#'     truncation error concentrated near the raster boundary. For
#'     small rasters or boundary-sensitive applications, this matters.
#'   \item Pycnophylactic tends to produce less ringing near sharp
#'     coarse-scale gradients than [disagg_cub()].
#'   \item Pycnophylactic is iterative and substantially slower
#'     [disagg_bl()], though cascading is typically much faster than
#'     one-shot pycno.
#' }
#'
#' @param coarse SpatRaster.
#' @param fact Integer disagg factor (>= 2).
#' @param cascade Disaggregation stages. One of:
#'   * `NULL` (default): auto-factorize `fact` into prime factors,
#'     sorted descending. For example, `fact = 48` becomes
#'     `c(3, 2, 2, 2, 2)`. Each stage's per-stage disaggregation is
#'     small, keeping the iteration in its comfortable regime where
#'     the 3x3 smoother is comparable in bandwidth to the per-stage
#'     block size. This is the recommended default.
#'   * `fact` (a single integer equal to `fact`): single-stage
#'     disaggregation. This is the original Tobler (1979) formulation;
#'     at large `fact` it converges slowly and may show within-block
#'     flattening (see "Caveat" below).
#'   * An integer vector with `prod(cascade) == fact`: explicit
#'     factorization. Useful for fine-tuning, e.g. `cascade = c(4, 4, 3)`
#'     for `fact = 48`.
#' @param max_iter Integer, maximum iterations per stage. Default 100.
#' @param tol Numeric, convergence tolerance as a fraction of the coarse
#'   raster's value range. Iteration stops when the largest per-cell
#'   change between successive iterations falls below `tol * range`.
#'   Default 1e-4. Controls only how close the result is to its smooth
#'   fixed point; mass preservation is enforced exactly at every iteration
#'   regardless of convergence.
#' @param smoother Character, smoothing kernel. One of:
#'   * `"laplacian_9"` (default): full 3x3 mean, including diagonals.
#'   * `"laplacian_5"`: 5-point pattern (cardinal neighbors + center).
#' @param variant Character, mean-restoration variant. One of:
#'   * `"additive"` (default): add a per-block constant. Works for any data.
#'   * `"multiplicative"`: scale within-block values. Strictly positive
#'     data only; will warn and fall back to additive if any block mean
#'     is non-positive.
#' @param init Character, initial fine raster at each stage. One of:
#'   * `"near"` (default): nearest-neighbor disagg. Matches the original
#'     Tobler formulation.
#'   * `"bilinear"`: standard bilinear disagg. A smoother starting point
#'     that typically converges faster.
#' @param na_fill Character, NA handling. Same semantics as [disagg_bl()];
#'   `"auto"` (default) uses `"fill"` when the coarse raster has NAs,
#'   `"reflect"` otherwise.
#' @param verbose Logical, print iteration progress. Default `FALSE`.
#'
#' @return Fine SpatRaster, with attributes `cascade` (the stage sequence
#'   actually used), `iterations` (integer vector, one per stage), and
#'   `converged` (logical vector, one per stage).
#'
#' @section Caveat for single-stage runs:
#' When you explicitly set `cascade = fact` (or `fact` is prime, leaving
#' no cascade options), the algorithm runs as a single pycnophylactic
#' iteration. At large `fact` (> ~20 with the default 3x3 smoother), the
#' single-stage fixed point can show subtle within-block flattening:
#' high-frequency content within a coarse block becomes piecewise-constant
#' at block scale, because the smoother cannot bridge across blocks in a
#' single stage. This is a structural property of the iterative algorithm
#' on regular grids, not a convergence issue — running more iterations
#' doesn't remove it.
#'
#' The default cascade behavior mostly eliminates this issue by keeping
#' each stage in a regime where the smoother's bandwidth matches the
#' per-stage block size. For `fact` with a large prime factor (e.g.
#' `fact = 22` = 2 x 11), the 11-stage will still exhibit the artifact;
#' consider [disagg_cub()] in that case.
#'
#' @references
#' Tobler, W. R. (1979). Smooth pycnophylactic interpolation for
#' geographical regions. *J. Am. Stat. Assoc.* 74(367), 519-530.
#'
#' @seealso [disagg_bl()] and [disagg_cub()] for faster, non-iterative
#'   pre-sharpening alternatives. The `pycno` package (`pycno::pycno`)
#'   implements the polygon-source version of the underlying algorithm.
#'
#' @examples
#' library(terra)
#' coarse <- rast(nrows = 30, ncols = 30, vals = runif(900))
#' fine_pyc <- disagg_pyc(coarse, fact = 6) # cascade = c(3, 2)
#' back <- aggregate(fine_pyc, fact = 6, fun = "mean")
#' max(abs(values(coarse) - values(back))) # ~ tolerance level
#'
#' # Single-stage -- original Tobler 1979 algorithm
#' fine_tobler <- disagg_pyc(coarse, fact = 6, cascade = 6)
#' @export
disagg_pyc <- function(coarse, fact,
                       cascade = NULL,
                       max_iter = 100L,
                       tol = 1e-4,
                       smoother = c("laplacian_9", "laplacian_5"),
                       variant = c("additive", "multiplicative"),
                       init = c("near", "bilinear"),
                       na_fill = c("auto", "reflect", "fill"),
                       verbose = FALSE) {
      smoother <- match.arg(smoother)
      variant  <- match.arg(variant)
      init     <- match.arg(init)
      na_fill  <- match.arg(na_fill)
      .validate_inputs_pyc(coarse, fact, max_iter, tol)
      fact <- as.integer(fact)
      max_iter <- as.integer(max_iter)

      cascade <- .resolve_cascade(cascade, fact)

      if (verbose && length(cascade) > 1L) {
            cat(sprintf("Cascading: %s = %d\n",
                        paste(cascade, collapse = " x "), fact))
      }

      if (na_fill == "auto") {
            na_fill <- if (any(is.na(terra::values(coarse)))) "fill" else "reflect"
      }

      # NA padding (once, on the original coarse)
      if (na_fill == "fill") {
            coarse_work <- .pad_with_nearest(coarse, fact)
            crop_back <- TRUE
      } else {
            coarse_work <- coarse
            crop_back <- FALSE
      }

      # Suppress terra progress bars during iteration
      old_progress <- terra::terraOptions(print = FALSE)$progress
      terra::terraOptions(progress = 0)
      on.exit(terra::terraOptions(progress = old_progress), add = TRUE)

      # Run cascade stages
      current <- coarse_work
      iterations <- integer(length(cascade))
      converged  <- logical(length(cascade))
      for (k in seq_along(cascade)) {
            stage_fact <- cascade[k]
            if (verbose) {
                  cat(sprintf("Stage %d/%d (fact = %d):\n",
                              k, length(cascade), stage_fact))
            }
            stage <- .disagg_pyc_stage(
                  coarse_stage = current,
                  fact         = stage_fact,
                  max_iter     = max_iter,
                  tol          = tol,
                  smoother     = smoother,
                  variant      = variant,
                  init         = init,
                  verbose      = verbose
            )
            current        <- stage$fine
            iterations[k]  <- stage$iterations
            converged[k]   <- stage$converged
      }

      fine <- current

      # Crop back to original fine extent if we padded
      if (crop_back) {
            fine_target_ext <- terra::ext(terra::disagg(coarse, fact = fact,
                                                        method = "near"))
            fine <- terra::crop(fine, fine_target_ext)
      }

      names(fine) <- names(coarse)
      attr(fine, "cascade")    <- cascade
      attr(fine, "iterations") <- iterations
      attr(fine, "converged")  <- converged
      fine
}

# Run one pycno stage (smooth + restore until convergence) at the given
# disagg factor. Returns list(fine, iterations, converged).
.disagg_pyc_stage <- function(coarse_stage, fact, max_iter, tol,
                              smoother, variant, init, verbose) {
      # Multiplicative variant validation
      if (variant == "multiplicative") {
            coarse_vals <- terra::values(coarse_stage)
            if (any(coarse_vals <= 0, na.rm = TRUE)) {
                  warning("Multiplicative variant requires strictly positive ",
                          "coarse values; some are <= 0. Falling back to additive.")
                  variant <- "additive"
            }
      }

      rng <- diff(range(terra::values(coarse_stage), na.rm = TRUE))
      if (rng == 0) rng <- 1
      abs_tol <- tol * rng

      # Initial guess
      fine <- terra::disagg(coarse_stage, fact = fact, method = init)

      # Smoother kernel
      w <- switch(smoother,
                  "laplacian_9" = matrix(1/9, 3, 3),
                  "laplacian_5" = {
                        k <- matrix(0, 3, 3)
                        k[2, ] <- 1/5
                        k[, 2] <- 1/5
                        k[2, 2] <- 1/5
                        k
                  })

      prev_fine <- fine
      converged <- FALSE
      for (it in seq_len(max_iter)) {
            fine <- terra::focal(fine, w = w, fun = "sum",
                                 na.policy = "omit", expand = TRUE)
            block_mean <- terra::aggregate(fine, fact = fact, fun = "mean")
            if (variant == "additive") {
                  correction_coarse <- coarse_stage - block_mean
                  correction_fine <- terra::disagg(correction_coarse,
                                                   fact = fact, method = "near")
                  fine <- fine + correction_fine
            } else {
                  ratio_coarse <- coarse_stage / block_mean
                  ratio_fine <- terra::disagg(ratio_coarse, fact = fact,
                                              method = "near")
                  fine <- fine * ratio_fine
            }

            delta <- max(abs(terra::values(fine) - terra::values(prev_fine)),
                         na.rm = TRUE)
            if (verbose) {
                  cat(sprintf("  iter %3d:  max|delta| = %.3e  (tol = %.3e)\n",
                              it, delta, abs_tol))
            }
            if (delta < abs_tol) {
                  converged <- TRUE
                  break
            }
            prev_fine <- fine
      }

      if (!converged) {
            warning(sprintf(
                  "Pycnophylactic iteration (stage fact = %d) did not converge in %d iterations (last delta = %s, tolerance = %s). Block means are still preserved exactly; the result may simply be less smooth than requested.",
                  fact, max_iter, signif(delta, 3), signif(abs_tol, 3)
            ))
      }

      list(fine = fine, iterations = it, converged = converged)
}

# Resolve the cascade argument to a vector of disagg factors.
#   NULL    -> reverse prime factorization of `fact`
#   integer -> validated and returned as-is
.resolve_cascade <- function(cascade, fact) {
      if (is.null(cascade)) {
            return(.prime_factors_desc(fact))
      }
      if (!is.numeric(cascade) || any(is.na(cascade))) {
            stop("`cascade` must be NULL or an integer vector with all elements >= 2")
      }
      cascade_int <- as.integer(cascade)
      if (any(cascade != cascade_int) || any(cascade_int < 2L)) {
            stop("`cascade` must be NULL or an integer vector with all elements >= 2")
      }
      if (prod(cascade_int) != fact) {
            stop(sprintf("prod(cascade) = %d does not equal fact = %d",
                         as.integer(prod(cascade_int)), fact))
      }
      cascade_int
}

# Reverse prime factorization: returns prime factors of n in decreasing
# order. Examples:
#   2  -> 2
#   12 -> c(3, 2, 2)
#   48 -> c(3, 2, 2, 2, 2)
#   100 -> c(5, 5, 2, 2)
#   11 -> 11  (prime)
.prime_factors_desc <- function(n) {
      n <- as.integer(n)
      factors <- integer(0)
      d <- 2L
      while (d * d <= n) {
            while (n %% d == 0L) {
                  factors <- c(factors, d)
                  n <- n %/% d
            }
            d <- d + 1L
      }
      if (n > 1L) factors <- c(factors, n)
      sort(factors, decreasing = TRUE)
}

.validate_inputs_pyc <- function(coarse, fact, max_iter, tol) {
      if (!inherits(coarse, "SpatRaster"))
            stop("`coarse` must be a SpatRaster")
      if (length(fact) != 1L || fact < 2L || fact != as.integer(fact))
            stop("`fact` must be an integer >= 2")
      if (length(max_iter) != 1L || max_iter < 1L ||
          max_iter != as.integer(max_iter))
            stop("`max_iter` must be a positive integer")
      if (!is.numeric(tol) || length(tol) != 1L || tol <= 0)
            stop("`tol` must be a positive numeric scalar")
}
