###############################################################################
##                                                                           ##
##   FOREST STRUCTURAL COMPLEXITY — APPLY THE PUBLISHED MODEL (Script 2B)     ##
##                                                                           ##
##   Produces a wall-to-wall Forest Structural Complexity Index (FSCI) map    ##
##   from airborne lidar, using the published model of Light et al. (Priest  ##
##   River Experimental Forest study, northern Idaho mixed-conifer forest).  ##
##                                                                           ##
##   -----------------------------------------------------------------       ##
##   APPLICABILITY                                                            ##
##   This model was trained at Priest River Experimental Forest (northern    ##
##   Idaho mixed-conifer). It is most reliable in structurally and           ##
##   ecologically comparable forests. This script does NOT mask predictions; ##
##   instead it outputs a separate "confidence" layer (Mahalanobis distance  ##
##   of each pixel's lidar metrics from the model's training data). Pixels    ##
##   with high Mahalanobis distance fall outside the range the model was      ##
##   trained on and should be interpreted with caution. Review that layer    ##
##   before drawing conclusions.                                             ##
##   -----------------------------------------------------------------       ##
##                                                                           ##
##   WHAT YOU PROVIDE:                                                        ##
##     - a folder of airborne lidar tiles (.laz/.las)                        ##
##     - a digital terrain model (DTM) covering the area                     ##
##     - a study-area boundary shapefile                                     ##
##     - the trained model file  (FSCI_model.rds)                            ##
##                                                                           ##
##   WHAT THIS PRODUCES (GeoTIFFs + a summary):                              ##
##     - FSCI prediction (0-100)                                             ##
##     - prediction uncertainty (half-width of 95% prediction interval)      ##
##     - lower / upper 95% prediction bounds                                 ##
##     - confidence layer (Mahalanobis distance; higher = less reliable)     ##
##                                                                           ##
###############################################################################


###############################################################################
# 0. SETUP                                                                     #
###############################################################################

required_packages <- c("lidR", "terra", "sf")
to_install <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_packages, library, character.only = TRUE))

# lidRmetrics is not on CRAN; install from GitHub if missing.
if (!("lidRmetrics" %in% installed.packages()[, "Package"])) {
  if (!("remotes" %in% installed.packages()[, "Package"])) install.packages("remotes")
  remotes::install_github("ptompalski/lidRmetrics")
}
library(lidRmetrics)
# NOTE: the zentropy metric can differ slightly across lidRmetrics versions.
# To reproduce a specific published map exactly, match the package version used
# to create it (see the reproducibility note in the README).


###############################################################################
# 1. CONFIGURATION — edit these paths and settings, then run top to bottom    #
###############################################################################

## -- Inputs ----------------------------------------------------------------
lidar_folder   <- "path/to/lidar_tiles_folder"      # folder of .laz/.las tiles
dtm_path       <- "path/to/your_dtm.tif"            # DTM raster covering the area
boundary_path  <- "path/to/study_area_boundary.shp" # study-area boundary shapefile
model_path     <- "path/to/FSCI_model.rds"          # the trained model (shipped)

## -- Output ----------------------------------------------------------------
output_dir     <- "path/to/output_folder"           # where results are written
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## -- Spatial settings ------------------------------------------------------
target_crs     <- "EPSG:26911"    # coordinate system for outputs (NAD83 UTM 11N)
grid_res       <- 30              # prediction grid cell size, metres (DO NOT change:
                                  # the model was trained at 30 m)
height_min     <- 0               # height filter, metres (as in training)
height_max     <- 60              # height filter, metres (as in training)

## -- Model training statistics (embedded; do not edit) ---------------------
## These reproduce the Mahalanobis confidence layer exactly. From Light et al.
train_center <- c(pz_above_30 = 9.9927863,
                  lad_cv      = 0.9766774,
                  zentropy    = 0.8114259)

train_cov <- matrix(
  c(112.4039128, -1.4762245,  0.4533199,
     -1.4762245,  0.1922812, -0.0011702,
      0.4533199, -0.0011702,  0.0235231),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("pz_above_30","lad_cv","zentropy"),
                  c("pz_above_30","lad_cv","zentropy"))
)

# Chi-squared cutoff (3 predictors) marking the edge of the training range.
# Pixels above this are flagged as low-confidence in the summary.
mahal_threshold <- qchisq(0.999, df = 3)


###############################################################################
# 2. METRIC FUNCTION — computes the three model predictors                     #
#    (identical to the definitions used to train the model)                    #
###############################################################################

fsci_metrics <- function(Z) {
  if (length(Z) == 0) {
    return(list(pz_above_30 = NA_real_, lad_cv = NA_real_, zentropy = NA_real_))
  }
  list(
    pz_above_30 = (sum(Z > 30) / length(Z)) * 100,   # PERCENT of returns above 30 m
    lad_cv      = lidRmetrics::metrics_lad(Z)$lad_cv,
    zentropy    = lidRmetrics::metrics_dispersion(Z)$zentropy
  )
}


###############################################################################
# 3. READ LIDAR, CLIP, NORMALIZE, COMPUTE METRICS                              #
###############################################################################

# Load the trained model
best_mod <- readRDS(model_path)

# Read lidar tiles as a catalog
ctg <- readLAScatalog(lidar_folder)
if (is.na(crs(ctg))) crs(ctg) <- target_crs

# Read and align the boundary
boundary <- sf::st_read(boundary_path, quiet = TRUE)
boundary <- sf::st_transform(boundary, crs = sf::st_crs(ctg))

# -- Normalize heights against the supplied DTM, KEEPING A CATALOG ----------
# The catalog is normalized tile-by-tile to temporary files (an output
# template is required for catalog-level processing). This preserves the
# original pipeline's structure (catalog + read filter), which is important
# for exact reproduction of the distribution-based zentropy metric.
dtm <- terra::rast(dtm_path)
opt_output_files(ctg) <- paste0(tempdir(), "/norm_{ORIGINALFILENAME}")
ctg_norm <- normalize_height(ctg, dtm)

# -- Height filter as a READ filter on the normalized catalog ---------------
# Matches the original "-keep_z 0 60" behaviour (applied at read time during
# metric computation), rather than filtering points already in memory.
opt_filter(ctg_norm) <- paste("-keep_z", height_min, height_max)

# Optionally restrict processing to the study-area boundary. The read filter
# plus the metric computation are applied across the catalog.
# (Metrics are computed wall-to-wall; clip the result to the boundary below.)

# -- Compute the three predictor metrics on a 30 m grid ---------------------
opt_output_files(ctg_norm) <- ""   # metrics returned in memory as a raster
metrics <- pixel_metrics(ctg_norm, ~fsci_metrics(Z), res = grid_res)

# Clip metrics to the study-area boundary
metrics <- terra::mask(metrics, terra::vect(boundary))
metrics <- terra::crop(metrics, terra::vect(boundary))


###############################################################################
# 4. PREDICT FSCI + PREDICTION INTERVALS                                       #
###############################################################################

# To a data frame of pixel coordinates + predictors
metrics_df <- as.data.frame(metrics, xy = TRUE, na.rm = TRUE)

# Point prediction
metrics_df$predicted <- predict(best_mod, newdata = metrics_df)

# 95% prediction interval (exact, using the shipped model object)
pred_int <- predict(best_mod, newdata = metrics_df,
                    interval = "prediction", level = 0.95)
metrics_df$lower       <- pred_int[, "lwr"]
metrics_df$upper       <- pred_int[, "upr"]
metrics_df$uncertainty <- (metrics_df$upper - metrics_df$lower) / 2


###############################################################################
# 5. CONFIDENCE LAYER — Mahalanobis distance (NOT a mask)                       #
###############################################################################

raster_preds <- metrics_df[, c("pz_above_30", "lad_cv", "zentropy")]
metrics_df$mahal_dist <- mahalanobis(raster_preds,
                                     center = train_center[c("pz_above_30","lad_cv","zentropy")],
                                     cov    = train_cov)
metrics_df$low_confidence <- metrics_df$mahal_dist > mahal_threshold


###############################################################################
# 6. RASTERIZE AND WRITE OUTPUTS                                               #
###############################################################################

make_raster <- function(df, col) {
  rast(df[, c("x", "y", col)], type = "xyz", crs = target_crs)
}

pred_r   <- make_raster(metrics_df, "predicted")
unc_r    <- make_raster(metrics_df, "uncertainty")
lower_r  <- make_raster(metrics_df, "lower")
upper_r  <- make_raster(metrics_df, "upper")
mahal_r  <- make_raster(metrics_df, "mahal_dist")

writeRaster(pred_r,  file.path(output_dir, "FSCI_prediction.tif"),   overwrite = TRUE)
writeRaster(unc_r,   file.path(output_dir, "FSCI_uncertainty.tif"),  overwrite = TRUE)
writeRaster(lower_r, file.path(output_dir, "FSCI_lower95.tif"),      overwrite = TRUE)
writeRaster(upper_r, file.path(output_dir, "FSCI_upper95.tif"),      overwrite = TRUE)
writeRaster(mahal_r, file.path(output_dir, "FSCI_confidence_mahalanobis.tif"),
            overwrite = TRUE)


###############################################################################
# 7. SUMMARY                                                                   #
###############################################################################

n_total <- nrow(metrics_df)
n_low   <- sum(metrics_df$low_confidence, na.rm = TRUE)

cat("\n================ FSCI MAPPING COMPLETE ================\n")
cat("Predicted pixels:            ", n_total, "\n")
cat("Mean FSCI:                   ", round(mean(metrics_df$predicted, na.rm = TRUE), 1), "\n")
cat("FSCI range:                  ",
    round(min(metrics_df$predicted, na.rm = TRUE), 1), "to",
    round(max(metrics_df$predicted, na.rm = TRUE), 1), "\n")
cat("Low-confidence pixels:       ", n_low,
    sprintf("(%.1f%% of area — outside training range)\n", 100 * n_low / n_total))
cat("Mean prediction uncertainty: ",
    round(mean(metrics_df$uncertainty, na.rm = TRUE), 2), "(+/- FSCI units)\n")
cat("Outputs written to:          ", output_dir, "\n")
cat("======================================================\n")
cat("\nNOTE: predictions are NOT masked. Consult FSCI_confidence_mahalanobis.tif;\n")
cat("higher values indicate pixels outside the model's training range, where\n")
cat("FSCI should be interpreted with caution.\n")
