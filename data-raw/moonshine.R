
# Generates vignettes/articles/moonshine.tif from topoclimate.pred::moonshine.
# Run manually when the source data changes; otherwise the tif is static.

library(terra)
truth <- rast(topoclimate.pred::moonshine)
terra::writeRaster(truth,
                   "vignettes/articles/moonshine.tif",
                   overwrite = TRUE,
                   gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
