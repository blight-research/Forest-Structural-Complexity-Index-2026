###############################################################################
##                                                                           ##
##   FOREST STRUCTURAL COMPLEXITY — STEP 1                                    ##
##   Identifying field variables that distinguish canopy-complexity groups   ##
##                                                                           ##
##   Adapted from the framework of Light et al. (Priest River Experimental   ##
##   Forest study). This script clusters plots on lidar-derived canopy       ##
##   complexity metrics, then identifies which field-measured structural     ##
##   variables best distinguish those groups.                                ##
##                                                                           ##
##   -----------------------------------------------------------------       ##
##   APPLICABILITY NOTE                                                       ##
##   This framework was developed in northern Idaho mixed-conifer forest.    ##
##   The *method* is transferable, but results (which variables are          ##
##   selected) are specific to the forest and data you supply. Use your      ##
##   own field and lidar data; do not assume the PREF variable set applies   ##
##   to a structurally or ecologically different forest.                     ##
##   -----------------------------------------------------------------       ##
##                                                                           ##
##   WHAT YOU PROVIDE:  one CSV containing, for each plot, the lidar         ##
##   canopy-complexity metrics AND the candidate field variables.           ##
##                                                                           ##
##   WHAT THIS PRODUCES: a table of all candidate field variables with       ##
##   their selection statistics, flagging which were retained. Written to    ##
##   CSV in the output folder.                                               ##
##                                                                           ##
###############################################################################


###############################################################################
# 0. SETUP — install/load required packages                                   #
###############################################################################

required_packages <- c(
  "dplyr", "tidyr", "cluster", "factoextra", "clusterCrit",
  "mclust", "randomForest"
)

to_install <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(42)   # ensures reproducible clustering criteria, RF, and bootstrap


###############################################################################
# 1. CONFIGURATION — the only section most users need to edit                 #
###############################################################################

## -- 1a. Input file --------------------------------------------------------
## One CSV with one row per plot. Must contain: a plot ID column, the lidar
## complexity metric columns, and the candidate field-variable columns.
input_csv <- "path/to/your_plot_data.csv"          # <-- EDIT

## -- 1b. Output folder -----------------------------------------------------
output_dir <- "path/to/output_folder"              # <-- EDIT
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## -- 1c. Column names ------------------------------------------------------
plot_id_col      <- "Plot"                          # <-- EDIT: plot identifier

## Lidar canopy-complexity metrics used to FORM the groups.
## (PREF used three metrics after Kane et al. 2010: canopy roughness,
##  95th-percentile height, and canopy density.)
complexity_metrics <- c("rumple.mean", "zq95", "canopydensity2")   # <-- EDIT

## The metric used to ORDER groups from least to most complex
## (group 1 = lowest, group N = highest). Must be one of complexity_metrics.
ordering_metric <- "zq95"                           # <-- EDIT if desired

## Candidate FIELD variables to be evaluated as group discriminators.
## List every field variable column name you want considered.
field_variables <- c(
  "field_var_1", "field_var_2", "field_var_3"       # <-- EDIT: your columns
)

## -- 1d. Optional row filter ----------------------------------------------
## If your CSV contains plots you want to exclude, set a filter here.
## Leave as NULL to use all rows. Example: filter_expr <- quote(Type == "FIA")
filter_expr <- NULL                                 # <-- EDIT or leave NULL

## -- 1e. Method settings (defaults reproduce the published framework) ------
k_min            <- 2      # smallest cluster number considered
k_max            <- 15     # largest cluster number considered
n_boot           <- 200    # bootstrap iterations for stability selection
topk_proportion  <- 0.20   # "top 20%" definition for stability frequency
topk_threshold   <- 0.20   # min. probability of appearing in top-K to retain


###############################################################################
# 2. LOAD AND PREPARE DATA                                                     #
###############################################################################

dat <- read.csv(input_csv, header = TRUE, stringsAsFactors = FALSE)

# Optional row filter
if (!is.null(filter_expr)) {
  dat <- dat[eval(filter_expr, envir = dat), , drop = FALSE]
}

# Basic checks so problems surface early with clear messages
needed_cols <- c(plot_id_col, complexity_metrics, field_variables)
missing_cols <- setdiff(needed_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("These columns are missing from the input CSV: ",
       paste(missing_cols, collapse = ", "))
}
if (!ordering_metric %in% complexity_metrics) {
  stop("`ordering_metric` must be one of `complexity_metrics`.")
}

dat[[plot_id_col]] <- as.factor(dat[[plot_id_col]])


###############################################################################
# 3. HIERARCHICAL CLUSTERING ON LIDAR COMPLEXITY METRICS                       #
###############################################################################

# Standardize the complexity metrics so none dominates by scale
clust_input <- scale(dat[, complexity_metrics, drop = FALSE])

# Ward's method on Euclidean distance
clust_dist <- dist(clust_input)
clust_tree <- hclust(clust_dist, method = "ward.D2")


###############################################################################
# 4. CHOOSE NUMBER OF CLUSTERS (ensemble of four criteria, averaged)           #
###############################################################################

# (1) Silhouette
sil_scores <- sapply(k_min:k_max, function(k) {
  cl <- cutree(clust_tree, k)
  if (length(unique(cl)) < 2) return(NA_real_)
  mean(cluster::silhouette(cl, clust_dist)[, 3])
})
k_sil <- (k_min:k_max)[which.max(sil_scores)]

# (2) Gap statistic
gap       <- clusGap(clust_input, FUN = hcut, K.max = k_max, B = 500)
k_gap     <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")

# (3) Calinski-Harabasz
ch_scores <- sapply(k_min:k_max, function(k) {
  intCriteria(as.matrix(clust_input),
              as.integer(cutree(clust_tree, k)),
              "calinski_harabasz")[[1]]
})
k_ch      <- (k_min:k_max)[which.max(ch_scores)]

# (4) Model-based (Mclust)
mc        <- Mclust(clust_input, G = 1:k_max, verbose = FALSE)
k_mclust  <- mc$G

# Ensemble: average of the four, rounded
k_selection <- data.frame(
  criterion = c("Gap", "Silhouette", "Calinski-Harabasz", "Mclust"),
  k         = c(k_gap, k_sil, k_ch, k_mclust)
)
n_clusters <- round(mean(k_selection$k))

cat("\nCluster-number selection by criterion:\n")
print(k_selection)
cat("Chosen number of clusters (rounded mean):", n_clusters, "\n\n")


###############################################################################
# 5. ASSIGN GROUPS AND ORDER THEM BY COMPLEXITY                                #
###############################################################################

raw_groups <- cutree(clust_tree, k = n_clusters)

# Order groups by mean of the ordering metric: 1 = least complex ... N = most.
group_means <- tapply(dat[[ordering_metric]], raw_groups, mean, na.rm = TRUE)
group_rank  <- rank(group_means, ties.method = "first")   # raw group -> rank
dat$complexity_group <- factor(group_rank[as.character(raw_groups)],
                               levels = sort(unique(group_rank)))

cat("Plots per complexity group (1 = least complex):\n")
print(table(dat$complexity_group))
cat("\n")


###############################################################################
# 6. KRUSKAL-WALLIS FILTER — which field variables differ across groups?       #
###############################################################################

kw_results <- do.call(rbind, lapply(field_variables, function(v) {
  test <- kruskal.test(dat[[v]] ~ dat$complexity_group)
  data.frame(variable = v,
             H       = unname(test$statistic),
             df      = unname(test$parameter),
             p_value = test$p.value,
             stringsAsFactors = FALSE)
}))
kw_results$adj_p_value <- p.adjust(kw_results$p_value, method = "fdr")

significant_vars <- kw_results$variable[kw_results$adj_p_value <= 0.05]

if (length(significant_vars) < 2) {
  stop("Fewer than two field variables passed the Kruskal-Wallis filter; ",
       "cannot proceed to Random Forest. Check your data or variable list.")
}


###############################################################################
# 7. RANDOM FOREST + BOOTSTRAP STABILITY SELECTION                             #
###############################################################################

# Data for the RF: significant field variables + the group response
rf_data <- dat[, c(significant_vars), drop = FALSE]
rf_data$complexity_group <- droplevels(dat$complexity_group)
rf_data <- rf_data[complete.cases(rf_data), , drop = FALSE]

predictors <- setdiff(names(rf_data), "complexity_group")
top_k      <- max(1, round(length(predictors) * topk_proportion))

# Bootstrap: stratified resampling, record importance each iteration
boot_importance <- vector("list", n_boot)
pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
for (i in seq_len(n_boot)) {
  boot_data <- rf_data %>%
    group_by(complexity_group) %>%
    sample_frac(1, replace = TRUE) %>%
    ungroup()

  rf_boot <- randomForest(
    complexity_group ~ .,
    data       = boot_data,
    importance = TRUE,
    ntree      = 500,
    mtry       = floor(sqrt(length(predictors)))
  )
  boot_importance[[i]] <- randomForest::importance(rf_boot, type = 1)[, 1]
  setTxtProgressBar(pb, i)
}
close(pb)

# Assemble importance matrix (variables x iterations)
imp_mat   <- do.call(cbind, boot_importance)
var_names <- names(boot_importance[[1]])
rownames(imp_mat) <- var_names

# Stability statistics
mean_imp <- rowMeans(imp_mat, na.rm = TRUE)
cv_imp   <- apply(imp_mat, 1, function(x) {
  mu <- mean(x); if (mu == 0) NA else sd(x) / mu
})

# Frequency each variable lands in the top-K by importance across iterations
topk_freq <- setNames(numeric(length(var_names)), var_names)
for (i in seq_len(ncol(imp_mat))) {
  top_vars <- names(sort(imp_mat[, i], decreasing = TRUE))[seq_len(top_k)]
  topk_freq[top_vars] <- topk_freq[top_vars] + 1
}
topk_prob <- topk_freq / ncol(imp_mat)


###############################################################################
# 8. FINAL SELECTION AND OUTPUT TABLE                                          #
###############################################################################

# Retention criterion: above-median mean importance AND top-K probability
# at or above the threshold (this is the framework's stability rule).
results <- data.frame(
  variable        = var_names,
  mean_importance = as.numeric(mean_imp),
  CV              = as.numeric(cv_imp),
  topK_prob       = as.numeric(topk_prob),
  stringsAsFactors = FALSE
)
results$kruskal_adj_p <- kw_results$adj_p_value[match(results$variable,
                                                      kw_results$variable)]
results$selected <- with(results,
                         mean_importance > median(mean_importance, na.rm = TRUE) &
                         topK_prob >= topk_threshold)

# Field variables that did NOT pass the Kruskal-Wallis filter appear here too,
# marked as not selected, so the table shows the full candidate pool.
not_tested <- setdiff(field_variables, results$variable)
if (length(not_tested) > 0) {
  filler <- data.frame(
    variable        = not_tested,
    mean_importance = NA_real_,
    CV              = NA_real_,
    topK_prob       = NA_real_,
    kruskal_adj_p   = kw_results$adj_p_value[match(not_tested, kw_results$variable)],
    selected        = FALSE,
    stringsAsFactors = FALSE
  )
  results <- rbind(results, filler)
}

# Order: selected first, then by importance
results <- results[order(-results$selected, -results$mean_importance,
                         na.last = TRUE), ]

selected_vars <- results$variable[results$selected]

cat("\n=== SELECTED FIELD VARIABLES ===\n")
print(selected_vars)
cat("\nFull selection table:\n")
print(results, row.names = FALSE)

# Write outputs
write.csv(results,
          file.path(output_dir, "field_variable_selection.csv"),
          row.names = FALSE)
write.csv(k_selection,
          file.path(output_dir, "cluster_number_selection.csv"),
          row.names = FALSE)

# Also save the per-plot group assignments (useful for Step 2)
plot_groups <- data.frame(
  Plot            = dat[[plot_id_col]],
  complexity_group = dat$complexity_group
)
write.csv(plot_groups,
          file.path(output_dir, "plot_complexity_groups.csv"),
          row.names = FALSE)

cat("\nDone. Outputs written to:\n  ", output_dir, "\n")
cat("  - field_variable_selection.csv   (main result)\n")
cat("  - cluster_number_selection.csv\n")
cat("  - plot_complexity_groups.csv\n")
