###############################################################################
##                                                                           ##
##   FOREST STRUCTURAL COMPLEXITY — BUILD YOUR OWN MODEL (Script 2A)          ##
##                                                                           ##
##   The full framework of Light et al. (Priest River Experimental Forest    ##
##   study). Constructs a field-derived Forest Structural Complexity Index    ##
##   (FSCI) from field variables, then builds and validates a model that      ##
##   predicts that index from airborne-lidar metrics.                         ##
##                                                                           ##
##   USE THIS when you have BOTH field-measured structural data and lidar     ##
##   metrics at the same plots and want a model tailored to your own forest.  ##
##   If instead you want to apply the published PREF model directly, use      ##
##   Script 2B.                                                              ##
##                                                                           ##
##   WORKFLOW (runs top to bottom after editing the CONFIG block):           ##
##     1. Build the FSCI index from your chosen field variables              ##
##     2. Split into training / validation plots                             ##
##     3. Filter and reduce candidate lidar metrics (RF + LASSO)             ##
##     4. Detect informative interactions (Friedman's H)                     ##
##     5. Best-subsets regression -> a table of candidate models             ##
##     6. YOU review the table and choose a model (parsimony + performance)  ##
##     7. Fit, validate, and save the final model                            ##
##                                                                           ##
##   NOTE: this framework involves several modelling steps and can take       ##
##   time (the interaction step is the slowest). That is expected.           ##
##                                                                           ##
###############################################################################


###############################################################################
# 0. SETUP                                                                     #
###############################################################################

required_packages <- c("dplyr", "caret", "randomForest", "glmnet",
                       "iml", "leaps")
to_install <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_packages, library, character.only = TRUE))


###############################################################################
# 1. CONFIGURATION                                                             #
###############################################################################

## -- Input / output --------------------------------------------------------
## One CSV, one row per plot, containing: plot ID, a grouping column that marks
## training vs. validation plots, the field variables used to build the index,
## and the candidate lidar metrics.
input_csv   <- "path/to/your_plot_data.csv"          # <-- EDIT
output_dir  <- "path/to/output_folder"               # <-- EDIT
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## -- Column roles ----------------------------------------------------------
plot_id_col <- "Plot"          # <-- plot identifier column

## Column that identifies the plot group. Training and validation are selected
## from this column by the values you give below.
group_col       <- "Type"       # e.g. "FIA" (training) vs "Perm" (validation)
train_group     <- "FIA"        # value marking TRAINING plots
valid_group     <- "Perm"       # value marking VALIDATION plots

## Optional: a second validation set (e.g. a later measurement year). If you
## do not have one, leave valid2_col = NULL.
valid2_col      <- "Meas_year"  # column distinguishing the 2nd validation set
valid2_value    <- 2019         # value marking the 2nd validation set (or NULL)
##   Set valid2_col <- NULL if you have only one validation set.

## -- FSCI index field variables (from Script 1) ----------------------------
## The field variables your Script-1 analysis retained. The index is the mean
## of these variables after each is min-max normalized to 0-1, scaled to 0-100.
index_field_vars <- c("field_var_1", "field_var_2", "field_var_3")  # <-- EDIT

## Reference endpoints for normalization. "data" uses the min and max observed
## across all plots. "manual" uses values you supply (representative min/max).
index_endpoints  <- "data"      # "data" or "manual"

## If index_endpoints = "manual", supply named min (young) and max (old) vectors
## covering every variable in index_field_vars:
index_min <- c(field_var_1 = NA, field_var_2 = NA, field_var_3 = NA)  # <-- EDIT if manual
index_max <- c(field_var_1 = NA, field_var_2 = NA, field_var_3 = NA)  # <-- EDIT if manual

## -- Candidate lidar metrics -----------------------------------------------
## The lidar metric columns to consider as predictors of the index.
lidar_metrics <- c("metric_1", "metric_2", "metric_3")  # <-- EDIT
## You may append specific complexity-related metrics you consider important
## even if they would otherwise be dropped by the correlation filter, e.g.:
##   force_keep_metrics <- c("rumple.mean", "zentropy")
force_keep_metrics <- character(0)

## -- Method settings (defaults reproduce the published framework) -----------
cor_cutoff  <- 0.85    # drop lidar metrics correlated above this
H_threshold <- 0.01    # Friedman's H threshold for keeping an interaction
nvmax       <- 25      # max model size for best-subsets
seed        <- 123     # reproducibility

## -- FINAL MODEL SELECTION -------------------------------------------------
## After you run once and inspect the candidate table (printed to console and
## saved to candidate_models.csv), set your chosen predictors here and re-run
## from section 8 onward. Leave NULL on the first pass.
chosen_predictors <- NULL       # e.g. c("pz_above_30", "zentropy", "lad_cv")


###############################################################################
# 2. LOAD DATA + BUILD THE FSCI INDEX                                          #
###############################################################################

set.seed(seed)
dat <- read.csv(input_csv, header = TRUE, stringsAsFactors = FALSE)

# Checks
needed <- c(plot_id_col, group_col, index_field_vars, lidar_metrics)
missing <- setdiff(needed, names(dat))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

# Reference endpoints
if (index_endpoints == "data") {
  ref_min <- sapply(dat[index_field_vars], min, na.rm = TRUE)
  ref_max <- sapply(dat[index_field_vars], max, na.rm = TRUE)
} else {
  ref_min <- index_min[index_field_vars]
  ref_max <- index_max[index_field_vars]
  if (any(is.na(ref_min)) || any(is.na(ref_max)))
    stop("index_endpoints = 'manual' but index_min/index_max contain NA.")
}

# Min-max normalize each field variable to 0-1, then FSCI = 100/n * sum
norm_mat <- sapply(index_field_vars, function(v) {
  abs((dat[[v]] - ref_min[v]) / (ref_max[v] - ref_min[v]))
})
dat$FSCI <- (100 / length(index_field_vars)) * rowSums(norm_mat)

cat("FSCI index constructed from", length(index_field_vars),
    "field variables. Range:",
    round(min(dat$FSCI, na.rm = TRUE), 1), "-",
    round(max(dat$FSCI, na.rm = TRUE), 1), "\n\n")


###############################################################################
# 3. TRAIN / VALIDATION SPLIT                                                  #
###############################################################################

keep_cols <- c(plot_id_col, "FSCI", lidar_metrics)

train_df <- dat[dat[[group_col]] == train_group, keep_cols, drop = FALSE]
train_df <- train_df[complete.cases(train_df), ]

valid_df <- dat[dat[[group_col]] == valid_group, keep_cols, drop = FALSE]
if (!is.null(valid2_col)) {
  valid_df  <- valid_df[dat[[valid2_col]][match(valid_df[[plot_id_col]],
                          dat[[plot_id_col]])] != valid2_value | is.na(
                          dat[[valid2_col]][match(valid_df[[plot_id_col]],
                          dat[[plot_id_col]])]), , drop = FALSE]
}
valid_df <- valid_df[complete.cases(valid_df), ]

valid2_df <- NULL
if (!is.null(valid2_col)) {
  v2 <- dat[dat[[group_col]] == valid_group &
              dat[[valid2_col]] == valid2_value, keep_cols, drop = FALSE]
  valid2_df <- v2[complete.cases(v2), ]
}


###############################################################################
# 4. FILTER CANDIDATE LIDAR METRICS                                            #
###############################################################################

lidar_cols <- lidar_metrics

# Drop zero-variance metrics
zv <- lidar_cols[sapply(train_df[lidar_cols], function(x) sd(x, na.rm = TRUE) == 0)]
lidar_cols <- setdiff(lidar_cols, zv)

# Drop highly correlated metrics (above cor_cutoff)
cor_mat  <- cor(train_df[lidar_cols], use = "pairwise.complete.obs")
high_cor <- caret::findCorrelation(cor_mat, cutoff = cor_cutoff)
if (length(high_cor) > 0) lidar_cols <- lidar_cols[-high_cor]

# Add back any metrics you chose to force-keep
lidar_cols <- unique(c(lidar_cols, force_keep_metrics))


###############################################################################
# 5. RANDOM FOREST -> TOP PREDICTORS                                           #
###############################################################################

set.seed(seed)
rf_form <- as.formula(paste("FSCI ~", paste(lidar_cols, collapse = " + ")))

rf_full <- randomForest(rf_form, data = train_df, ntree = 2000, importance = TRUE)
varImp_full <- caret::varImp(rf_full, scale = FALSE)

# Optimal number of predictors by training R^2 across top-k
imp_order <- rownames(varImp_full)[order(-varImp_full$Overall)]
r2_by_k <- sapply(seq_len(min(8, length(imp_order))), function(k) {
  fm <- as.formula(paste("FSCI ~", paste(imp_order[1:k], collapse = " + ")))
  set.seed(seed)
  m  <- randomForest(fm, data = train_df, ntree = 1000)
  cor(train_df$FSCI, predict(m, train_df))^2
})
best_k   <- which.max(r2_by_k)
top_vars <- imp_order[1:best_k]

# Reduced RF (used by the interaction step below)
red_form   <- as.formula(paste("FSCI ~", paste(top_vars, collapse = " + ")))
set.seed(seed)
rf_reduced <- randomForest(red_form, data = train_df, ntree = 1000, importance = TRUE)


###############################################################################
# 6. LASSO -> ADDITIONAL MAIN EFFECTS                                          #
###############################################################################

set.seed(seed)
x        <- model.matrix(FSCI ~ . - 1, train_df[, c("FSCI", lidar_cols)])
x_scaled <- scale(x)
y        <- train_df$FSCI

lasso_mod <- glmnet::cv.glmnet(x_scaled, y, alpha = 1, standardize = FALSE)
lasso_co  <- coef(lasso_mod, s = "lambda.min")
lasso_hit <- rownames(lasso_co)[as.numeric(lasso_co) != 0]
lasso_hit <- setdiff(lasso_hit, "(Intercept)")

# Main-effect pool = RF top vars + LASSO-selected vars
main_effects <- unique(c(top_vars, lasso_hit))
main_effects <- intersect(main_effects, lidar_cols)   # keep valid columns only


###############################################################################
# 7. FRIEDMAN'S H -> INFORMATIVE INTERACTIONS  (slowest step)                   #
###############################################################################

pred_iml <- iml::Predictor$new(model = rf_reduced,
                               data  = train_df[, top_vars, drop = FALSE],
                               y     = train_df$FSCI)

pairwise_H <- function(pred, v1, v2, grid_size = 20) {
  f1  <- iml::FeatureEffect$new(pred, feature = v1, method = "pdp",
                                grid.size = grid_size)$results
  f2  <- iml::FeatureEffect$new(pred, feature = v2, method = "pdp",
                                grid.size = grid_size)$results
  f12 <- iml::FeatureEffect$new(pred, feature = c(v1, v2), method = "pdp",
                                grid.size = grid_size)$results
  a <- f1$.value  - mean(f1$.value)
  b <- f2$.value  - mean(f2$.value)
  ab <- f12$.value - mean(f12$.value)
  sqrt(max((var(ab) - (var(a) + var(b))) / var(ab), 0))
}

inter_vars <- main_effects
interactions <- list()
if (length(inter_vars) >= 2) {
  set.seed(seed)
  for (i in 1:(length(inter_vars) - 1)) {
    for (j in (i + 1):length(inter_vars)) {
      Hval <- tryCatch(pairwise_H(pred_iml, inter_vars[i], inter_vars[j]),
                       error = function(e) NA_real_)
      interactions[[length(interactions) + 1]] <-
        data.frame(var1 = inter_vars[i], var2 = inter_vars[j], H = Hval)
    }
  }
}
interaction_df <- if (length(interactions)) do.call(rbind, interactions) else
  data.frame(var1 = character(), var2 = character(), H = numeric())

selected_interactions <- interaction_df[!is.na(interaction_df$H) &
                                          interaction_df$H > H_threshold, , drop = FALSE]
if (nrow(selected_interactions) > 0) {
  selected_interactions$inter_name <- paste0(selected_interactions$var1, "_X_",
                                             selected_interactions$var2)
}

# Add interaction columns to each dataset
add_interactions <- function(df, inter_df) {
  if (nrow(inter_df) == 0) return(df)
  for (k in seq_len(nrow(inter_df))) {
    df[[inter_df$inter_name[k]]] <- df[[inter_df$var1[k]]] * df[[inter_df$var2[k]]]
  }
  df
}
train_df2  <- add_interactions(train_df,  selected_interactions)
valid_df2  <- add_interactions(valid_df,  selected_interactions)
valid2_df2 <- if (!is.null(valid2_df)) add_interactions(valid2_df, selected_interactions) else NULL


###############################################################################
# 8. BEST-SUBSETS -> CANDIDATE MODEL TABLE                                     #
###############################################################################

all_predictors <- c(main_effects,
                     if (nrow(selected_interactions)) selected_interactions$inter_name)

full_form <- as.formula(paste("FSCI ~", paste(all_predictors, collapse = " + ")))
best_sub  <- leaps::regsubsets(full_form, data = train_df2, nbest = 1,
                               nvmax = min(nvmax, length(all_predictors)),
                               method = "exhaustive")
bs_sum    <- summary(best_sub)

# Validation performance helper
val_perf <- function(vars, vdf) {
  if (is.null(vdf)) return(c(R2 = NA, RMSE = NA))
  fm <- as.formula(paste("FSCI ~", paste(vars, collapse = " + ")))
  m  <- lm(fm, data = train_df2)
  p  <- predict(m, newdata = vdf)
  c(R2 = cor(vdf$FSCI, p)^2, RMSE = sqrt(mean((vdf$FSCI - p)^2)))
}

# Build one row per candidate model size
candidate_table <- do.call(rbind, lapply(seq_len(nrow(bs_sum$which)), function(sz) {
  vars <- names(which(bs_sum$which[sz, ]))[-1]   # drop intercept
  v1 <- val_perf(vars, valid_df2)
  v2 <- val_perf(vars, valid2_df2)
  data.frame(
    size        = sz,
    predictors  = paste(vars, collapse = " + "),
    adjR2       = round(bs_sum$adjr2[sz], 4),
    BIC         = round(bs_sum$bic[sz], 2),
    valid1_R2   = round(v1["R2"], 3),
    valid1_RMSE = round(v1["RMSE"], 3),
    valid2_R2   = round(v2["R2"], 3),
    valid2_RMSE = round(v2["RMSE"], 3),
    stringsAsFactors = FALSE, row.names = NULL
  )
}))

cat("\n================ CANDIDATE MODELS ================\n")
cat("Review this table and choose a model balancing parsimony (fewer\n")
cat("predictors, lower BIC) with performance (higher adjR2 and validation R2).\n")
cat("Then set `chosen_predictors` in the CONFIG block and re-run from here.\n\n")
print(candidate_table, row.names = FALSE)
write.csv(candidate_table, file.path(output_dir, "candidate_models.csv"),
          row.names = FALSE)


###############################################################################
# 9. FIT + VALIDATE + SAVE THE CHOSEN MODEL                                     #
#    (runs only after you set `chosen_predictors` in CONFIG)                    #
###############################################################################

if (!is.null(chosen_predictors)) {

  final_form <- as.formula(paste("FSCI ~", paste(chosen_predictors, collapse = " + ")))
  final_model <- lm(final_form, data = train_df2)

  cat("\n================ FINAL MODEL ================\n")
  print(summary(final_model))

  report <- function(df, label) {
    if (is.null(df)) return(invisible())
    p <- predict(final_model, newdata = df)
    cat(sprintf("%-22s R2 = %.3f   RMSE = %.3f\n",
                label, cor(df$FSCI, p)^2, sqrt(mean((df$FSCI - p)^2))))
  }
  cat("\nPerformance:\n")
  report(train_df2,  "Training:")
  report(valid_df2,  "Validation (set 1):")
  report(valid2_df2, "Validation (set 2):")

  saveRDS(final_model, file.path(output_dir, "FSCI_model.rds"))
  cat("\nFinal model saved to:", file.path(output_dir, "FSCI_model.rds"), "\n")
  cat("Use this .rds with Script 2B to map FSCI across your area.\n")

} else {
  cat("\n(No final model fit yet: set `chosen_predictors` in CONFIG and re-run.)\n")
}
