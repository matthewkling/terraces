library(terra)

test_that("bilinear round-trip kernel sums to 1 across factors", {
      for (fact in c(2L, 3L, 5L, 10L)) {
            K <- roundtrip_kernel("bilinear", fact)
            expect_equal(sum(K), 1, tolerance = 1e-12, info = paste("fact =", fact))
            expect_equal(rowSums(K), colSums(K), tolerance = 1e-12)
      }
})

test_that("inverse kernel inverts round-trip on the interior", {
      K_inv <- kernel("bilinear", fact = 5, radius = 9)
      expect_lt(attr(K_inv, "conv_max_err"), 1e-10)
})

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

test_that("kernels are cached and reused", {
      .clear_kernel_cache()
      expect_length(ls(.kernel_cache), 0)
      K1 <- kernel("bilinear", 5, 7)
      expect_length(ls(.kernel_cache), 1)
      K2 <- kernel("bilinear", 5, 7)
      expect_identical(K1, K2)
      expect_length(ls(.kernel_cache), 1)
})

test_that("multi-layer rasters disaggregate correctly", {
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

test_that("disagg_pyc preserves block means", {
      set.seed(42)
      n <- 30; fact <- 5
      r <- rast(nrows = n, ncols = n, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      xy <- xyFromCell(r, 1:ncell(r))
      values(r) <- with(as.data.frame(xy),
                        sin(6 * x) * cos(5 * y) + 0.5 * sin(8 * x * y))

      fine <- disagg_pyc(r, fact = fact, max_iter = 200, tol = 1e-5)
      back <- aggregate(fine, fact = fact, fun = "mean")
      err <- values(r) - values(back)
      data_range <- diff(range(values(r)))
      # Tobler preserves block means very precisely by construction (each
      # iteration restores them exactly).
      expect_lt(max(abs(err), na.rm = TRUE), 1e-3 * data_range)
})

test_that("list_methods includes both methods", {
      m <- list_methods()
      expect_true("bilinear" %in% m$name)
      expect_true("pycnophylactic" %in% m$name)
      expect_true("prefilter" %in% m$type)
      expect_true("iterative" %in% m$type)
})
