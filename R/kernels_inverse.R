# Inverse kernel computation
#
# For prefilter methods, the inverse kernel K_inv satisfies K * K_inv = delta
# (on the kernel's support window). Computed by direct linear-system solve;
# cached in a session-level environment keyed on (method, fact, radius).

.kernel_cache <- new.env(parent = emptyenv())
.cache_key <- function(method, fact, radius) {
      paste(method, fact, radius, sep = "/")
}

# Methods that have a prefilter (i.e., support kernel()/roundtrip_kernel()).
# Pycnophylactic is iterative and is not in this list.
.prefilter_methods <- c("bilinear")

#' Compute the inverse kernel for a prefilter disaggregation method
#'
#' The inverse kernel is the focal filter applied to the coarse raster
#' before standard interpolation; it predistorts the coarse data so that
#' subsequent disagg produces a fine raster whose block means equal the
#' input coarse values.
#'
#' Kernels depend only on the method and disagg factor (not the data) and
#' are cached for the R session.
#'
#' @param method Character, prefilter method name. Currently `"bilinear"`.
#' @param fact Integer disagg factor.
#' @param radius Integer half-width of the kernel. Larger = more accurate
#'   at higher one-time cost. `NULL` uses a method-specific default.
#' @return A square matrix of side `2 * radius + 1`, with attributes
#'   `method`, `fact`, `radius`, `tail_max` (largest absolute value at
#'   the kernel boundary; an upper bound on truncation error), and
#'   `conv_max_err` (largest deviation of K * K_inv from the delta on
#'   the interior; should be machine epsilon).
#' @export
kernel <- function(method, fact, radius = NULL) {
      if (!method %in% .prefilter_methods) {
            stop("Method '", method, "' is not a prefilter method. ",
                 "Prefilter methods: ", paste(.prefilter_methods, collapse = ", "), ".")
      }
      fact <- as.integer(fact)
      if (is.null(radius)) radius <- .default_radius(method, fact)
      radius <- as.integer(radius)
      if (radius < 1L) stop("radius must be at least 1")

      key <- .cache_key(method, fact, radius)
      if (exists(key, envir = .kernel_cache, inherits = FALSE)) {
            return(get(key, envir = .kernel_cache, inherits = FALSE))
      }

      K <- roundtrip_kernel(method, fact)
      K_inv <- .compute_inverse(K, radius)
      attr(K_inv, "method") <- method
      attr(K_inv, "fact")   <- fact
      attr(K_inv, "radius") <- radius

      assign(key, K_inv, envir = .kernel_cache)
      K_inv
}

# Solve the finite-window deconvolution K * K_inv = delta.
.compute_inverse <- function(K, radius) {
      r <- as.integer(radius); sz <- 2L * r + 1L; n <- sz * sz
      rK <- (nrow(K) - 1L) / 2L

      K_off <- expand.grid(da = -rK:rK, db = -rK:rK)
      K_off$w <- mapply(function(da, db) K[db + rK + 1L, da + rK + 1L],
                        K_off$da, K_off$db)

      idx <- function(di, dj) (di + r) * sz + (dj + r) + 1L

      M <- matrix(0, n, n)
      e <- numeric(n); e[idx(0L, 0L)] <- 1

      for (di in -r:r) for (dj in -r:r) {
            row <- idx(di, dj)
            for (q in seq_len(nrow(K_off))) {
                  a <- K_off$da[q]; b <- K_off$db[q]; w <- K_off$w[q]
                  ii <- di - a; jj <- dj - b
                  if (abs(ii) <= r && abs(jj) <= r) {
                        M[row, idx(ii, jj)] <- M[row, idx(ii, jj)] + w
                  }
            }
      }
      x <- solve(M, e)

      K_inv <- matrix(0, sz, sz)
      for (di in -r:r) for (dj in -r:r) {
            K_inv[di + r + 1L, dj + r + 1L] <- x[idx(di, dj)]
      }

      conv <- .convolve_small(K, K_inv)
      delta <- matrix(0, sz, sz); delta[r + 1L, r + 1L] <- 1
      attr(K_inv, "conv_max_err") <- max(abs(conv - delta))
      attr(K_inv, "tail_max") <- max(abs(K_inv[1L, ]), abs(K_inv[sz, ]),
                                     abs(K_inv[, 1L]), abs(K_inv[, sz]))
      K_inv
}

.convolve_small <- function(K, K_inv) {
      sz <- nrow(K_inv); r <- (sz - 1L) / 2L; rK <- (nrow(K) - 1L) / 2L
      out <- matrix(0, sz, sz)
      for (di in -r:r) for (dj in -r:r) {
            s <- 0
            for (a in -rK:rK) for (b in -rK:rK) {
                  ii <- di - a; jj <- dj - b
                  if (abs(ii) <= r && abs(jj) <= r) {
                        s <- s + K[b + rK + 1L, a + rK + 1L] * K_inv[ii + r + 1L, jj + r + 1L]
                  }
            }
            out[di + r + 1L, dj + r + 1L] <- s
      }
      out
}

# Internal: clear cache (for tests).
.clear_kernel_cache <- function() {
      rm(list = ls(.kernel_cache), envir = .kernel_cache)
      invisible(NULL)
}
