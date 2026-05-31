# Round-trip kernel for bilinear disaggregation
#
# K is the linear operator on coarse cells representing "disaggregate by
# bilinear, then aggregate back by block mean." Derived analytically by
# enumerating fine-cell positions within one representative coarse block
# and computing bilinear's weights pulling from the surrounding 3x3
# coarse stencil. Depends only on the disagg factor.

.roundtrip_bilinear <- function(fact) {
      k <- as.integer(fact)
      K <- matrix(0, 3, 3)
      centers <- (seq_len(k) - 0.5) / k

      for (a in centers) for (b in centers) {
            if (a < 0.5) { ix <- c(-1, 0); wx <- c(0.5 - a, 0.5 + a) }
            else         { ix <- c( 0, 1); wx <- c(1.5 - a, a - 0.5) }
            if (b < 0.5) { iy <- c(-1, 0); wy <- c(0.5 - b, 0.5 + b) }
            else         { iy <- c( 0, 1); wy <- c(1.5 - b, b - 0.5) }
            for (pi in 1:2) for (pj in 1:2) {
                  K[iy[pi] + 2L, ix[pj] + 2L] <-
                        K[iy[pi] + 2L, ix[pj] + 2L] + wx[pj] * wy[pi]
            }
      }
      K / (k * k)
}

# Default inverse-kernel radius for bilinear; scales gently with fact.
# At fact = 5, radius = 7 gives ~1e-5 truncation error.
.default_radius_bilinear <- function(fact) max(7L, as.integer(fact) + 2L)

# Round-trip kernel for Keys cubic convolution disaggregation (a = -0.5).
#
# Keys cubic is the convention used by terra::resample(method = "cubic")
# (verified empirically; matches impulse response to machine precision).
# Each fine cell pulls from a 4x4 stencil of coarse cells; the round-trip
# operator is a 5x5 kernel.

.keys_cubic_1d <- function(t, a = -0.5) {
      at <- abs(t)
      ifelse(at < 1,
             (a + 2) * at^3 - (a + 3) * at^2 + 1,
             ifelse(at < 2,
                    a * at^3 - 5 * a * at^2 + 8 * a * at - 4 * a,
                    0))
}

.roundtrip_cubic <- function(fact) {
      k <- as.integer(fact)
      K <- matrix(0, 5, 5)
      centers <- (seq_len(k) - 0.5) / k

      for (u in centers) for (v in centers) {
            # Which 4 coarse cells does cubic pull from in each axis?
            if (u < 0.5) { ix <- c(-2, -1, 0, 1) } else { ix <- c(-1, 0, 1, 2) }
            if (v < 0.5) { iy <- c(-2, -1, 0, 1) } else { iy <- c(-1, 0, 1, 2) }

            # Distances from fine cell to each surrounding coarse cell center
            dx <- u - (ix + 0.5); dy <- v - (iy + 0.5)
            wx <- .keys_cubic_1d(dx); wy <- .keys_cubic_1d(dy)
            wx <- wx / sum(wx); wy <- wy / sum(wy)   # safety renormalization

            # Accumulate into the 5x5 stencil
            for (pi in 1:4) for (pj in 1:4) {
                  K[iy[pi] + 3L, ix[pj] + 3L] <-
                        K[iy[pi] + 3L, ix[pj] + 3L] + wx[pj] * wy[pi]
            }
      }
      K / (k * k)
}

# Cubic round-trip is more diffusive than bilinear (center weight ~0.71
# vs 0.58) and has negative side lobes, so the inverse kernel decays more
# slowly. Use a slightly larger default radius.
.default_radius_cubic <- function(fact) max(9L, as.integer(fact) + 4L)

#' Round-trip kernel for a disaggregation method
#'
#' For prefilter methods, returns the analytical kernel representing
#' "disagg then aggregate" as a small convolution on the coarse grid.
#' Useful for diagnostics and for understanding the bias structure of
#' standard interpolation.
#'
#' @param method Character. Currently `"bilinear"` (3x3 kernel) and
#'   `"cubic"` (5x5 kernel, Keys convention with a = -0.5) are supported
#'   as prefilter methods. (Pycnophylactic disaggregation is iterative
#'   and does not have a round-trip kernel in this sense.)
#' @param fact Integer disagg factor.
#' @return A square odd-sized numeric matrix.
#' @export
roundtrip_kernel <- function(method, fact) {
      switch(method,
             "bilinear" = .roundtrip_bilinear(as.integer(fact)),
             "cubic"    = .roundtrip_cubic(as.integer(fact)),
             stop("Unknown prefilter method: '", method, "'. ",
                  "Currently supported: 'bilinear', 'cubic'.")
      )
}

# Default kernel radius for a given method/fact. Used by kernel() when
# the user passes radius = NULL.
.default_radius <- function(method, fact) {
      switch(method,
             "bilinear" = .default_radius_bilinear(fact),
             "cubic"    = .default_radius_cubic(fact),
             stop("Unknown method: '", method, "'")
      )
}

#' List disaggregation methods available in terraces
#'
#' @return A data frame with columns `name`, `type`, and `description`.
#' @export
list_methods <- function() {
      data.frame(
            name = c("bilinear", "cubic", "pycnophylactic"),
            type = c("prefilter", "prefilter", "iterative"),
            description = c(
                  "Bilinear via prefilter (disagg_bl)",
                  "Keys cubic via prefilter (disagg_cub)",
                  "Tobler 1979 pycnophylactic, iterative (disagg_pyc)"
            ),
            stringsAsFactors = FALSE
      )
}
