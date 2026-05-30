library(testthat)
library(terra)

test_that("bilinear round-trip kernel sums to 1 across factors", {
      for (fact in c(2L, 3L, 5L, 10L)) {
            K <- ces_roundtrip_kernel("bilinear", fact)
            expect_equal(sum(K), 1, tolerance = 1e-12, info = paste("fact =", fact))
            expect_equal(rowSums(K), colSums(K), tolerance = 1e-12)
      }
})

test_that("inverse kernel inverts round-trip on the interior", {
      K_inv <- ces_kernel("bilinear", fact = 5, radius = 9)
      expect_lt(attr(K_inv, "conv_max_err"), 1e-10)
})

test_that("ces_disagg_bl preserves block means to tol on smooth data", {
      # Smooth synthetic field (representative of real use cases) on a raster
      # large enough to fit the default radius without clamping.
      n <- 30; fact <- 5
      r <- rast(nrows = n, ncols = n, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      xy <- xyFromCell(r, 1:ncell(r))
      values(r) <- with(as.data.frame(xy),
                        sin(6 * x) * cos(5 * y) + 0.5 * sin(8 * x * y))

      fine <- ces_disagg_bl(r, fact = fact, radius = 9)
      back <- aggregate(fine, fact = fact, fun = "mean")
      err <- values(r) - values(back)
      # Round-trip preservation: ~0.5% of data range at radius 9, fact 5,
      # on smooth data. Larger fact and rougher data yield larger error.
      data_range <- diff(range(values(r)))
      expect_lt(max(abs(err), na.rm = TRUE), 0.01 * data_range)
})

test_that("ces_disagg_bl preserves block means on noisy data", {
      # Uncorrelated noise is the hardest case for the prefilter
      # (worst-case high-frequency content). Tolerance is correspondingly loose.
      set.seed(42)
      n <- 30; fact <- 5
      r <- rast(nrows = n, ncols = n, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                vals = runif(n * n, -5, 5))
      fine <- ces_disagg_bl(r, fact = fact, radius = 9)
      back <- aggregate(fine, fact = fact, fun = "mean")
      err <- values(r) - values(back)
      data_range <- diff(range(values(r)))
      expect_lt(max(abs(err), na.rm = TRUE), 0.05 * data_range)
})

test_that("ces_disagg_bl beats plain bilinear on round-trip preservation", {
      # The defining property of the method: CES round-trip error should be
      # substantially smaller than plain bilinear's on the same data.
      set.seed(42)
      r <- rast(nrows = 30, ncols = 30, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                vals = runif(900, -5, 5))

      fine_ces <- ces_disagg_bl(r, fact = 5, radius = 9)
      fine_bl  <- terra::disagg(r, fact = 5, method = "bilinear")

      back_ces <- aggregate(fine_ces, fact = 5, fun = "mean")
      back_bl  <- aggregate(fine_bl,  fact = 5, fun = "mean")

      err_ces <- max(abs(values(r) - values(back_ces)), na.rm = TRUE)
      err_bl  <- max(abs(values(r) - values(back_bl)),  na.rm = TRUE)

      # CES should be at least an order of magnitude better.
      expect_lt(err_ces, err_bl / 10)
})

test_that("kernels are cached and reused", {
      .ces_clear_cache()
      expect_length(ls(.ces_cache), 0)
      K1 <- ces_kernel("bilinear", 5, 7)
      expect_length(ls(.ces_cache), 1)
      K2 <- ces_kernel("bilinear", 5, 7)
      expect_identical(K1, K2)
      expect_length(ls(.ces_cache), 1)
})

test_that("multi-layer rasters disaggregate correctly", {
      set.seed(1)
      r <- c(
            rast(nrows = 25, ncols = 25, vals = runif(625)),
            rast(nrows = 25, ncols = 25, vals = runif(625))
      )
      names(r) <- c("a", "b")
      fine <- ces_disagg_bl(r, fact = 5, radius = 7)
      expect_equal(nlyr(fine), 2L)
      expect_equal(names(fine), c("a", "b"))
})

test_that("ces_disagg dispatches to method-specific function", {
      set.seed(1)
      r <- rast(nrows = 25, ncols = 25, vals = runif(625))
      out_direct   <- ces_disagg_bl(r, fact = 5, radius = 7)
      out_dispatch <- ces_disagg(r, fact = 5, method = "bilinear", radius = 7)
      expect_equal(values(out_direct), values(out_dispatch))
})

test_that("small raster triggers radius clamping with warning", {
      r <- rast(nrows = 8, ncols = 8, vals = runif(64))
      expect_warning(
            ces_disagg_bl(r, fact = 3, radius = 7),
            "radius reduced"
      )
})

test_that("too-small raster errors cleanly", {
      r <- rast(nrows = 4, ncols = 4, vals = runif(16))
      expect_error(ces_disagg_bl(r, fact = 2), "too small")
})

test_that("methods registry lists registered methods", {
      m <- ces_list_methods()
      expect_true("bilinear" %in% m$name)
})

test_that("ces_disagg_pyc preserves block means", {
      set.seed(42)
      n <- 30; fact <- 5
      r <- rast(nrows = n, ncols = n, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      xy <- xyFromCell(r, 1:ncell(r))
      values(r) <- with(as.data.frame(xy),
                        sin(6 * x) * cos(5 * y) + 0.5 * sin(8 * x * y))

      fine <- ces_disagg_pyc(r, fact = fact, max_iter = 200, tol = 1e-5)
      back <- aggregate(fine, fact = fact, fun = "mean")
      err <- values(r) - values(back)
      data_range <- diff(range(values(r)))
      # Tobler preserves block means very precisely by construction (each
      # iteration restores them exactly).
      expect_lt(max(abs(err), na.rm = TRUE), 1e-3 * data_range)
})

test_that("pycnophylactic method is registered", {
      m <- ces_list_methods()
      expect_true("pycnophylactic" %in% m$name)
      expect_true("iterative" %in% m$type)
})

test_that("ces_disagg dispatches to pycnophylactic correctly", {
      set.seed(1)
      r <- rast(nrows = 25, ncols = 25, vals = runif(625))
      out_direct   <- ces_disagg_pyc(r, fact = 5, max_iter = 20, tol = 1e-3)
      out_dispatch <- ces_disagg(r, fact = 5, method = "pycnophylactic",
                                 max_iter = 20, tol = 1e-3)
      expect_equal(values(out_direct), values(out_dispatch))
})

test_that("ces_register_method validates round-trip kernel", {
      expect_error(
            ces_register_method(
                  "bad_test",
                  disagg_fn      = function(coarse, fact, ...) coarse,
                  type           = "prefilter",
                  roundtrip_fn   = function(fact) matrix(0.5, 3, 3),  # doesn't sum to 1
                  default_radius = function(fact) 3L
            ),
            "sum to 1"
      )
})
