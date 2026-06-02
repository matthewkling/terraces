library(terra)

# Round-trip kernels ------------------------------------------------------

test_that("bilinear round-trip kernel sums to 1 across factors", {
      for (fact in c(2L, 3L, 5L, 10L)) {
            K <- roundtrip_kernel("bilinear", fact)
            expect_equal(sum(K), 1, tolerance = 1e-12, info = paste("fact =", fact))
            expect_equal(rowSums(K), colSums(K), tolerance = 1e-12)
      }
})

test_that("cubic round-trip kernel sums to 1 across factors", {
      for (fact in c(3L, 5L, 10L)) {
            K <- roundtrip_kernel("cubic", fact)
            expect_equal(sum(K), 1, tolerance = 1e-12, info = paste("fact =", fact))
            expect_equal(rowSums(K), colSums(K), tolerance = 1e-12)
      }
})

test_that("cubic round-trip kernel matches terra::resample empirically", {
      # Verifies that our analytical Keys cubic kernel (a = -0.5) matches
      # what terra::resample(method = "cubic") actually computes.
      fact <- 5L
      K_pred <- roundtrip_kernel("cubic", fact)

      # Build empirical round-trip kernel by impulse response.
      n_pad <- 4L
      n <- 2L * n_pad + 1L
      coarse <- rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
      v <- numeric(ncell(coarse))
      v[n_pad * n + n_pad + 1L] <- 1
      values(coarse) <- v

      fine_template <- disagg(coarse, fact = fact, method = "near")
      fine <- resample(coarse, fine_template, method = "cubic")
      back <- aggregate(fine, fact = fact, fun = "mean")

      mat <- matrix(values(back), nrow = n, byrow = TRUE)
      center <- n_pad + 1L
      K_emp <- mat[(center - 2L):(center + 2L), (center - 2L):(center + 2L)]

      expect_lt(max(abs(K_emp - K_pred)), 1e-7)
})


# Inverse kernels ---------------------------------------------------------

test_that("inverse kernel inverts round-trip on the interior", {
      K_inv <- kernel("bilinear", fact = 5, radius = 9)
      expect_lt(attr(K_inv, "conv_max_err"), 1e-10)
      K_inv_c <- kernel("cubic", fact = 5, radius = 9)
      expect_lt(attr(K_inv_c, "conv_max_err"), 1e-10)
})

test_that("kernels are cached and reused", {
      .clear_kernel_cache()
      expect_length(ls(.kernel_cache), 0)
      K1 <- kernel("bilinear", 5, 7)
      expect_length(ls(.kernel_cache), 1)
      K2 <- kernel("bilinear", 5, 7)
      expect_identical(K1, K2)
      expect_length(ls(.kernel_cache), 1)
})


# Bilinear disaggregation -------------------------------------------------

test_that("disagg_bl preserves block means to tol on smooth data", {
      # Smooth synthetic field (representative of real use cases) on a raster
      # large enough to fit the default radius without clamping.
      n <- 30; fact <- 5
      r <- rast(nrows = n, ncols = n, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      xy <- xyFromCell(r, 1:ncell(r))
      values(r) <- with(as.data.frame(xy),
                        sin(6 * x) * cos(5 * y) + 0.5 * sin(8 * x * y))

      fine <- disagg_bl(r, fact = fact, radius = 9)
      back <- aggregate(fine, fact = fact, fun = "mean")
      err <- values(r) - values(back)
      # Round-trip preservation: ~0.5% of data range at radius 9, fact 5,
      # on smooth data. Larger fact and rougher data yield larger error.
      data_range <- diff(range(values(r)))
      expect_lt(max(abs(err), na.rm = TRUE), 0.01 * data_range)
})

test_that("disagg_bl preserves block means on noisy data", {
      # Uncorrelated noise is the hardest case for the prefilter
      # (worst-case high-frequency content). Tolerance is correspondingly loose.
      set.seed(42)
      n <- 30; fact <- 5
      r <- rast(nrows = n, ncols = n, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                vals = runif(n * n, -5, 5))
      fine <- disagg_bl(r, fact = fact, radius = 9)
      back <- aggregate(fine, fact = fact, fun = "mean")
      err <- values(r) - values(back)
      data_range <- diff(range(values(r)))
      expect_lt(max(abs(err), na.rm = TRUE), 0.05 * data_range)
})

test_that("disagg_bl beats plain bilinear on round-trip preservation", {
      # The defining property of the method: round-trip error should be
      # substantially smaller than plain bilinear's on the same data.
      set.seed(42)
      r <- rast(nrows = 30, ncols = 30, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                vals = runif(900, -5, 5))

      fine_aces <- disagg_bl(r, fact = 5, radius = 9)
      fine_bl  <- terra::disagg(r, fact = 5, method = "bilinear")

      back_ces <- aggregate(fine_aces, fact = 5, fun = "mean")
      back_bl  <- aggregate(fine_bl,  fact = 5, fun = "mean")

      err_ces <- max(abs(values(r) - values(back_ces)), na.rm = TRUE)
      err_bl  <- max(abs(values(r) - values(back_bl)),  na.rm = TRUE)

      # disagg_bl should be at least an order of magnitude better.
      expect_lt(err_ces, err_bl / 10)
})

test_that("multi-layer rasters disaggregate correctly with bilinear", {
      set.seed(1)
      r <- c(
            rast(nrows = 25, ncols = 25, vals = runif(625)),
            rast(nrows = 25, ncols = 25, vals = runif(625))
      )
      names(r) <- c("a", "b")
      fine <- disagg_bl(r, fact = 5, radius = 7)
      expect_equal(nlyr(fine), 2L)
      expect_equal(names(fine), c("a", "b"))
})

test_that("small raster triggers radius clamping with warning", {
      r <- rast(nrows = 8, ncols = 8, vals = runif(64))
      expect_warning(
            disagg_bl(r, fact = 3, radius = 7),
            "radius reduced"
      )
})

test_that("too-small raster errors cleanly", {
      r <- rast(nrows = 4, ncols = 4, vals = runif(16))
      expect_error(disagg_bl(r, fact = 2), "too small")
})


# Cubic disaggregation ----------------------------------------------------

test_that("disagg_cub preserves block means on smooth data", {
      n <- 30; fact <- 5
      r <- rast(nrows = n, ncols = n, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      xy <- xyFromCell(r, 1:ncell(r))
      values(r) <- with(as.data.frame(xy),
                        sin(6 * x) * cos(5 * y) + 0.5 * sin(8 * x * y))

      fine <- disagg_cub(r, fact = fact, radius = 9)
      back <- aggregate(fine, fact = fact, fun = "mean")
      err <- values(r) - values(back)
      data_range <- diff(range(values(r)))
      # Cubic is more diffusive than bilinear, so truncation error at the same
      # radius is somewhat larger. Allow ~1% of data range.
      expect_lt(max(abs(err), na.rm = TRUE), 0.02 * data_range)
})

test_that("disagg_cub uses cubic interpolation end-to-end", {
      # Regression test for a previous bug where apply_kernel hardcoded
      # bilinear interpolation regardless of method. If the bug were
      # present, the cubic pre-sharpening kernel followed by bilinear
      # interpolation would not be mass-preserving, and the impulse
      # would not be recovered.
      n <- 15L
      coarse <- rast(nrows = n, ncols = n)
      v <- numeric(ncell(coarse))
      v[(n %/% 2L) * n + (n %/% 2L) + 1L] <- 1
      values(coarse) <- v

      fine <- disagg_cub(coarse, fact = 5L, radius = 5)
      back <- aggregate(fine, fact = 5L, fun = "mean")

      # Impulse should be recovered to within cubic's truncation error.
      expect_lt(max(abs(values(coarse) - values(back)), na.rm = TRUE), 0.02)
})

test_that("multi-layer rasters disaggregate correctly with cubic", {
      set.seed(1)
      r <- c(
            rast(nrows = 25, ncols = 25, vals = runif(625)),
            rast(nrows = 25, ncols = 25, vals = runif(625))
      )
      names(r) <- c("a", "b")
      fine <- disagg_cub(r, fact = 5, radius = 8)
      expect_equal(nlyr(fine), 2L)
      expect_equal(names(fine), c("a", "b"))
})


# Pycnophylactic disaggregation -------------------------------------------

test_that("disagg_pyc preserves block means to near machine precision", {
      # Pycno enforces block-mean equality at every iteration by construction;
      # the only error is floating-point rounding in terra's focal/aggregate.
      set.seed(42)
      r <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                vals = runif(400))
      fine <- disagg_pyc(r, fact = 5, max_iter = 50)
      back <- aggregate(fine, fact = 5, fun = "mean")
      expect_lt(max(abs(values(r) - values(back)), na.rm = TRUE), 1e-10)
})

test_that("disagg_pyc auto-cascades composite fact by default", {
      set.seed(1)
      r <- rast(nrows = 20, ncols = 20, vals = runif(400))
      fine <- disagg_pyc(r, fact = 6, max_iter = 50)
      expect_equal(attr(fine, "cascade"), c(3L, 2L))
      expect_length(attr(fine, "iterations"), 2L)
      expect_length(attr(fine, "converged"), 2L)
})

test_that("disagg_pyc handles prime fact as single stage", {
      set.seed(1)
      r <- rast(nrows = 20, ncols = 20, vals = runif(400))
      fine <- disagg_pyc(r, fact = 5, max_iter = 50)
      expect_equal(attr(fine, "cascade"), 5L)
      expect_length(attr(fine, "iterations"), 1L)
})

test_that("disagg_pyc accepts explicit cascade", {
      set.seed(1)
      r <- rast(nrows = 20, ncols = 20, vals = runif(400))
      fine <- disagg_pyc(r, fact = 12, cascade = c(4, 3), max_iter = 50)
      expect_equal(attr(fine, "cascade"), c(4L, 3L))
})

test_that("disagg_pyc single-stage via cascade = fact", {
      set.seed(1)
      r <- rast(nrows = 20, ncols = 20, vals = runif(400))
      fine <- disagg_pyc(r, fact = 6, cascade = 6, max_iter = 50)
      expect_equal(attr(fine, "cascade"), 6L)
      expect_length(attr(fine, "iterations"), 1L)
})

test_that("disagg_pyc rejects invalid cascade", {
      r <- rast(nrows = 20, ncols = 20, vals = runif(400))
      expect_error(disagg_pyc(r, fact = 6, cascade = c(2, 4)),
                   "does not equal fact")
      expect_error(disagg_pyc(r, fact = 6, cascade = c(1, 6)),
                   ">= 2")
      expect_error(disagg_pyc(r, fact = 6, cascade = c(2.5, 3)),
                   ">= 2")
})

test_that("disagg_pyc cascade preserves block means as exactly as single-stage", {
      set.seed(1)
      r <- rast(nrows = 20, ncols = 20, vals = runif(400))
      fine <- disagg_pyc(r, fact = 6, max_iter = 50)
      err <- max(abs(values(r) - values(aggregate(fine, 6, fun = "mean"))),
                 na.rm = TRUE)
      expect_lt(err, 1e-10)
})

test_that("multi-layer rasters disaggregate correctly with pycnophylactic", {
      set.seed(1)
      r <- c(
            rast(nrows = 20, ncols = 20, vals = runif(400)),
            rast(nrows = 20, ncols = 20, vals = runif(400))
      )
      names(r) <- c("a", "b")
      fine <- disagg_pyc(r, fact = 4, max_iter = 30)
      expect_equal(nlyr(fine), 2L)
      expect_equal(names(fine), c("a", "b"))
})

test_that("disagg_pyc verbose mode prints progress without erroring", {
      set.seed(1)
      r <- rast(nrows = 20, ncols = 20, vals = runif(400))
      expect_output(
            disagg_pyc(r, fact = 6, verbose = TRUE, max_iter = 30),
            "Stage|iter|Cascading"
      )
})


# Edge effects ------------------------------------------------------------

test_that("edge_effects returns expected geometry and value range", {
      r <- rast(nrows = 20, ncols = 20, vals = runif(400))
      m <- edge_effects(r, fact = 5, method = "bilinear")
      expect_s4_class(m, "SpatRaster")
      expect_equal(nrow(m), 100L)
      expect_equal(ncol(m), 100L)
      expect_equal(as.vector(ext(m)), as.vector(ext(r)))
      expect_true(all(values(m) %in% c(0L, 1L, 2L)))
})

test_that("edge_effects has interior zeros and boundary non-zeros", {
      r <- rast(nrows = 40, ncols = 40, vals = runif(1600))
      m <- edge_effects(r, fact = 5, method = "cubic")
      v <- values(m)
      expect_true(any(v == 0L))   # interior is clean
      expect_true(any(v > 0L))    # boundary is affected
})

test_that("edge_effects cubic zone is at least as wide as bilinear", {
      # Both the default radius and the interp_reach are larger for cubic.
      r <- rast(nrows = 40, ncols = 40, vals = runif(1600))
      m_bl  <- edge_effects(r, fact = 5, method = "bilinear")
      m_cub <- edge_effects(r, fact = 5, method = "cubic")
      expect_gte(mean(values(m_cub) > 0L), mean(values(m_bl) > 0L))
})

test_that("edge_effects rejects unsupported methods", {
      r <- rast(nrows = 20, ncols = 20, vals = runif(400))
      expect_error(edge_effects(r, fact = 5, method = "pycnophylactic"))
      expect_error(edge_effects(r, fact = 5, method = "near"))
})


# Method registry ---------------------------------------------------------

test_that("list_methods includes all methods", {
      m <- list_methods()
      expect_true("bilinear" %in% m$name)
      expect_true("cubic" %in% m$name)
      expect_true("pycnophylactic" %in% m$name)
      expect_true("prefilter" %in% m$type)
      expect_true("iterative" %in% m$type)
})
