# Round-trip kernel for bilinear disaggregation
#
# The round-trip kernel K is the linear operator on coarse cells representing
# "disaggregate by bilinear, then aggregate back by block mean." Derived
# analytically by enumerating fine-cell positions within one representative
# coarse block and computing bilinear's weights pulling from the surrounding
# 3x3 coarse stencil. The result depends only on the disagg factor.

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

# Default inverse-kernel radius: scales gently with disagg factor.
# At fact = 5, radius = 7 gives ~1e-5 truncation error.
.default_radius_bilinear <- function(fact) max(7L, as.integer(fact) + 2L)

#' Round-trip kernel for a registered disaggregation method
#'
#' For prefilter methods, returns the analytical kernel that represents
#' "disagg then aggregate" as a small convolution on the coarse grid.
#' Useful for diagnostics and for understanding the bias structure of
#' standard interpolation.
#'
#' @param method Character, registered method name.
#' @param fact Integer disagg factor.
#' @return A square odd-sized numeric matrix.
#' @export
ces_roundtrip_kernel <- function(method, fact) {
      m <- .get_method(method)
      if (is.null(m$roundtrip)) {
            stop("Method '", method, "' is not a prefilter method and has no ",
                 "round-trip kernel.")
      }
      m$roundtrip(as.integer(fact))
}
