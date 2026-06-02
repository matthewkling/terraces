#' Visualize edge-effect regions in a disaggregated raster
#'
#' Returns a fine-resolution mask showing which cells are affected by
#' boundary effects in [disagg_bl()] / [disagg_cub()] output. Useful
#' for masking out unreliable regions in boundary-sensitive applications,
#' or for inspecting the affected zone before processing.
#'
#' Edge effects come from two sources, encoded as values 0, 1, 2:
#' \itemize{
#'   \item `0`: Interior cells with no edge effects. Mass preservation
#'     and interpolation accuracy are both at full precision here.
#'   \item `1`: Cells where the pre-sharpening focal pass used
#'     reflective extension (directly or via the interpolation kernell
#'     reaching into the reflected band). Mass preservation is slightly
#'     degraded but interpolation is full-order.
#'   \item `2`: Cells additionally affected by terra's interpolation
#'     boundary handling, where the cubic or bilinear kernel would
#'     extend outside the raster and terra falls back to a lower-order
#'     scheme. This is a narrow band within ~0.5 coarse-cell-widths of
#'     the edge for bilinear, ~1.5 for cubic.
#' }
#'
#' The bands shrink to nothing in the interior, so for large rasters
#' the bulk of the output is `0`. For small rasters, the bands may
#' overlap and cover most or all of the output; the function still
#' returns a sensible mask in that case.
#'
#' @param coarse SpatRaster; the coarse input that would be passed to
#'   [disagg_bl()] or [disagg_cub()]. Only its geometry is used.
#' @param fact Integer disagg factor (>= 2).
#' @param method One of `"bilinear"` or `"cubic"`. Must match the
#'   method you intend to use.
#' @param radius Pre-sharpening kernel radius. If `NULL`, the default
#'   for the method is used (matches [disagg_bl()] / [disagg_cub()]
#'   defaults). For accurate masks, pass the same value you would pass
#'   to the disagg function.
#'
#' @return A single-layer fine-resolution SpatRaster with integer
#'   values 0, 1, or 2. Same geometry as the output of
#'   `disagg_*(coarse, fact)`. Attributes `method`, `radius`, and
#'   `interp_reach` are attached for inspection.
#'
#' @examples
#' library(terra)
#' coarse <- rast(matrix(runif(30 * 30), 30, 30))
#'
#' # See where edge effects would land for a fact=5 cubic disagg
#' m <- edge_effects(coarse, fact = 5, method = "cubic")
#' plot(m)
#'
#' # Use the mask to drop edge-affected cells from a disagg result
#' fine <- disagg_cub(coarse, fact = 5)
#' fine_clean <- mask(fine, m, maskvalues = c(1, 2))
#'
#' @export
edge_effects <- function(coarse, fact,
                         method = c("bilinear", "cubic"),
                         radius = NULL) {
      method <- match.arg(method)
      if (!inherits(coarse, "SpatRaster")) {
            stop("`coarse` must be a SpatRaster")
      }
      if (length(fact) != 1L || fact < 2L || fact != as.integer(fact)) {
            stop("`fact` must be an integer >= 2")
      }
      fact <- as.integer(fact)

      # Resolve radius (matches disagg_* defaults and clamping)
      if (is.null(radius)) {
            radius <- if (method == "bilinear") {
                  .default_radius_bilinear(fact)
            } else {
                  .default_radius_cubic(fact)
            }
      }
      radius <- as.integer(radius)

      min_dim <- min(nrow(coarse), ncol(coarse))
      r_hard <- (min_dim - 1L) %/% 2L
      r_soft <- floor((1/3) * min_dim)
      r_max <- min(r_hard, r_soft)
      if (radius > r_max) radius <- r_max

      # Stencil reach (in coarse-cell-widths). Theoretical worst case:
      # bilinear needs a 2x2 stencil → fits when d_edge >= 0.5
      # cubic needs a 4x4 stencil → fits when d_edge >= 1.5
      interp_reach <- if (method == "bilinear") 0.5 else 1.5

      # Fine grid geometry
      fine_template <- terra::disagg(coarse, fact = fact, method = "near")

      ext_c <- terra::ext(coarse)
      res_c <- terra::res(coarse)
      xy_fine <- terra::xyFromCell(fine_template,
                                   seq_len(terra::ncell(fine_template)))

      # Distance from each fine cell center to nearest edge, in coarse-cell widths
      d_edge <- pmin(
            (xy_fine[, 1] - terra::xmin(ext_c)) / res_c[1],
            (terra::xmax(ext_c) - xy_fine[, 1]) / res_c[1],
            (xy_fine[, 2] - terra::ymin(ext_c)) / res_c[2],
            (terra::ymax(ext_c) - xy_fine[, 2]) / res_c[2]
      )

      # Classify cells
      mask_values <- rep(0L, length(d_edge))
      mask_values[d_edge < radius + interp_reach] <- 1L
      mask_values[d_edge < interp_reach]          <- 2L

      result <- fine_template
      terra::values(result) <- mask_values
      names(result) <- "edge_effects"

      attr(result, "method")       <- method
      attr(result, "radius")       <- radius
      attr(result, "interp_reach") <- interp_reach

      result
}
