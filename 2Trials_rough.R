####Rough####
#define the regression method to be used at each neighborhood 
#Provide the optimum no of components with lowest RMSD in above NNV statistic in the dissimilarity matrix above (taken as min of RMSD)
no_components
#define the dissimilarity method from the list (pca, pls, cor, euclid, cosine, sid)
#my_diss = "pls"
#define how to use the dissimilarity information (ignore it)
#ignore_diss = "none"
#Building regression object with instructions to build PLS models
#my_model = local_fit_pls(pls_c = no_components)
#define the neighborhood sizes to test
#my_ks = seq(50, 500, by = 50)
#for the moment use only "NNv" validation (it will be faster)
#nnv_val_control = mbl_control(validation_type = "NNv")
#Use Yu argument to validate the predictions directly
Y2_hat = mbl(Xr = foo,Yr = foo1,Xu = X2,Yu = Y2,
             k = seq(50, 500, by = 50), method = local_fit_pls(pls_c = no_components),
             diss_method = "pls", diss_usage = "none", pc_selection = list(method = "opc", value = 40),
             control = mbl_control(validation_type = "NNv"), scale = TRUE)
#Y2_hat = mbl(Xr = foo,Yr = foo1,Xu = X2,Yu = Y2,k = my_ks,method = my_model,diss_method = my_diss,diss_usage = ignore_diss,pc_selection = list(method = "opc", value = 40),control = nnv_val_control,scale = TRUE)
Y2_hat
#Plot predictions in testing data
plot(Y2_hat, main = "Global2Local")
Y2_hat$validation_results$nearest_neighbor_validation
#collect predictions and get the indices of the best results according to nearest neighbor validation statistics
c_val_name = "validation_results"
c_nn_val_name = "nearest_neighbor_validation"
index_min = which.min(Y2_hat[[c_val_name]][[c_nn_val_name]]$rmse)
Y2_hat_preds = get_predictions(Y2_hat)[, ..index_min]
Y2_hat_preds = as.matrix(Y2_hat_preds)[,1]
# Calculate statistics
statistics_result <- calculate_statistics(Y2, Y2_hat_preds)

#For parallel implementation
#Running the mbl function using multiple cores
n_cores = detectCores()/2
clust = makeCluster(n_cores)
registerDoParallel(clust)
Y2_hat = mbl(Xr = X,Yr = Y,Xu = X2,Yu = Y2,k = my_ks,method = my_model,diss_method = my_diss,diss_usage = ignore_diss,control = nnv_val_control,scale = TRUE)
stopCluster(clust)
Y2_hat

# Calculate statistics
statistics_result <- calculate_statistics(Y, Y2_hat)


####MBL source code####
function (Xr, Yr, Xu = NULL, Yu = NULL, k, k_diss, k_range, spike = NULL, method = local_fit_wapls(min_pls_c = 3, max_pls_c = min(dim(Xr), 15)), 
          diss_method = "pca", diss_usage = "predictors", 
          gh = TRUE, pc_selection = list(method = "opc", value = min(dim(Xr), 40)), control = mbl_control(), group = NULL, center = TRUE, 
          scale = FALSE, verbose = TRUE, documentation = character(), 
          seed = NULL, ...) 
{
  f_call <- match.call()
  "%mydo%" <- get("%do%")
  if (control$allow_parallel & getDoParRegistered()) {
    "%mydo%" <- get("%dopar%")
  }
  if (!is.logical(verbose)) {
    stop("'verbose' must be logical")
  }
  if (missing(k)) {
    k <- NULL
  }
  if (missing(k_diss)) {
    k_diss <- NULL
  }
  if (missing(k_range)) {
    k_range <- NULL
  }
  input_dots <- list(...)
  ini_cntrl <- control
  ortho_diss_methods <- c("pca", "pca.nipals", "pls")
  if (".local" %in% names(input_dots)) {
    if (isTRUE(input_dots$.local)) {
      if (!"pre_k" %in% names(input_dots)) {
        stop("When '.local = TRUE', argument 'pre_k' needs to be provided. See ortho_diss documentation")
      }
      if (!is.null(k)) {
        if (input_dots$pre_k < max(k)) {
          stop("'k' cannot be larger than 'pre_k'")
        }
      }
    }
  }
  if (!is.logical(center)) {
    stop("'center' argument must be logical")
  }
  if (!is.logical(scale)) {
    stop("'scale' argument must be logical")
  }
  if (missing(Xu)) {
    stop("Xu is missing")
  }
  if (!is.null(Xu)) {
    if (ncol(Xr) != ncol(Xu)) {
      stop("The number of predictor variables in Xr must be equal to the number of variables in Xu")
    }
  }
  if (ncol(Xr) < 4) {
    stop("This function works only with matrices containing more than 3 predictor variables")
  }
  if (length(Yr) != nrow(Xr)) {
    stop("length(Yr) must be equal to nrow(Xr)")
  }
  if (!is.null(Xu)) {
    if (any(is.na(Yr))) {
      stop("The current version of the mbl function does not handle NAs in the response variable of the reference observations (Yr)")
    }
  }
  Xr <- as.matrix(Xr)
  Yr <- as.matrix(Yr)
  rownames(Xr) <- 1:nrow(Xr)
  if (is.null(colnames(Xr))) {
    colnames(Xr) <- 1:ncol(Xr)
  }
  validation_type <- control$validation_type
  is_local_cv <- "local_cv" %in% validation_type
  is_nnv_val <- "NNv" %in% validation_type
  if (all(c("local_cv", "NNv") %in% control$validation_type)) {
    validation_type <- "both"
  }
  if (!is.null(Xu)) {
    pre_nms_ng <- "Xu_"
    n_xu <- ln <- nrow(Xu)
    Xu <- as.matrix(Xu)
    rownames(Xu) <- 1:nrow(Xu)
    if (is.null(colnames(Xu))) {
      colnames(Xu) <- 1:ncol(Xu)
    }
    if (sum(!colnames(Xu) == colnames(Xr)) != 0) {
      stop("Variable names in Xr do not match those in Xu")
    }
    if (validation_type %in% c("NNv", "both") & nrow(Xu) < 
        3) {
      stop(paste0("For nearest neighbor validation (control$validation_type == 'NNv')", 
                  " Xu must contain at least 3 observations"))
    }
    if (!is.null(Yu)) {
      Yu <- as.matrix(Yu)
      if (length(Yu) != nrow(Xu)) {
        stop("Number of observations in Yu and Xu differ")
      }
    }
    constellation <- FALSE
    first_nn <- 1
    observed <- Yu
    y_output_name <- "yu_obs"
    y_hat_output_name <- "pred"
    val_summary_name <- "Yu_prediction_statistics"
  }
  else {
    pre_nms_ng <- "Xr_"
    ln <- nrow(Xr)
    n_xu <- 0
    constellation <- TRUE
    first_nn <- 2
    observed <- Yr
    y_output_name <- "yr_obs"
    y_hat_output_name <- "fitted"
    val_summary_name <- "Yr_fitted_statistics"
  }
  n_xr <- nrow(Xr)
  n_total <- n_xr + n_xu
  diss_methods <- c("pca", "pca.nipals", "pls", "cor", "euclid", 
                    "cosine", "sid")
  if (!is.character(diss_method) & !is.matrix(diss_method)) {
    mtds <- paste(diss_methods, collapse = ", ")
    stop(paste0("'diss_method' must be one of: ", mtds, " or a matrix"))
  }
  if (!is.null(group)) {
    if (length(group) != nrow(Xr)) {
      stop(paste0("The length of 'group' must be equal to the number of ", 
                  "observations in 'Xr'"))
    }
  }
  if (length(pc_selection) != 2 | !is.list(pc_selection)) {
    stop("'pc_selection' must be a list of length 2")
  }
  if (!all(names(pc_selection) %in% c("method", "value")) | 
      is.null(names(pc_selection))) {
    names(pc_selection)[sapply(pc_selection, FUN = is.character)] <- "method"
    names(pc_selection)[sapply(pc_selection, FUN = is.numeric)] <- "value"
  }
  pc_sel_method <- match.arg(pc_selection$method, c("opc", 
                                                    "var", "cumvar", "manual"))
  pc_threshold <- pc_selection$value
  if (pc_sel_method %in% c("opc", "manual") & pc_selection$value > 
      min(n_total, ncol(Xr))) {
    pc_threshold <- min(n_total, ncol(Xr), n_total)
    if (!is.null(Xu)) {
      message_pc <- paste0("When pc_selection$method is 'opc' or 'manual', the value ", 
                           "specified in \npc_selection$value cannot be larger than ", 
                           "min(nrow(Xr) + nrow(Xu), ncol(Xr)) \n(i.e ", 
                           pc_threshold, "). Therefore the value was reset to ", 
                           pc_threshold)
    }
    else {
      message_pc <- paste0("When pc_selection$method is 'opc' or 'manual', the value ", 
                           "specified in \npc_selection$value cannot be larger than ", 
                           "min(dim(Xr)) \n(i.e ", pc_threshold, "). Therefore the value was reset to ", 
                           pc_threshold)
    }
    warning(message_pc)
  }
  match.arg(diss_usage, c("predictors", "weights", "none"))
  if (is.null(k) & is.null(k_diss)) {
    stop("Either k or k_diss must be specified")
  }
  k_max <- NULL
  if (!is.null(k)) {
    if (!is.null(k_diss)) {
      stop("Only one of k or k_diss can be specified")
    }
    if (!is.numeric(k)) {
      stop("k must be a vector of integers")
    }
    else {
      k <- unique(sort(ceiling(k)))
    }
    k <- sort(k)
    k_max <- max(k)
  }
  k_diss_max <- NULL
  if (!is.null(k_diss)) {
    k_diss <- unique(sort(k_diss))
    if (is.null(k_range)) {
      stop("If the k_diss argument is used, k_range must be specified")
    }
    if (length(k_range) != 2 | !is.numeric(k_range) | any(k_range < 
                                                          1)) {
      stop("k_range must be a vector of length 2 which specifies the minimum (first value, larger than 0) and the maximum (second value) number of neighbors")
    }
    k_range <- sort(k_range)
    k_min_range <- as.integer(k_range[1])
    k_max_range <- as.integer(k_range[2])
    if (k_min_range < 4) {
      stop("Minimum number of nearest neighbors allowed is 4")
    }
    if (k_max_range > nrow(Xr)) {
      stop("Maximum number of nearest neighbors cannot exceed the number of reference observations")
    }
    k_diss_max <- max(k_diss)
  }
  if (".local" %in% names(input_dots)) {
    if (isTRUE(input_dots$local)) {
      if (!"pre_k" %in% names(input_dots)) {
        stop(paste0("When .local = TRUE (passed to the ortho_diss method), the ", 
                    "'pre_k' argument must be specified"))
      }
      if (input_dots$pre_k < k_max) {
        stop(paste0("pre_k must be larger than ", ifelse(is.null(k), 
                                                         "max(k_range)", "max(k)")))
      }
    }
  }
  if (!"local_fit" %in% class(method)) {
    stop("Object passed to method must be of class local_fit")
  }
  if (!is.null(k)) {
    k <- as.integer(k)
    if (min(k) < 4) {
      stop("Minimum number of nearest neighbors allowed is 3")
    }
    if (max(k) > nrow(Xr)) {
      stop(paste0("The number of nearest neighbors cannot exceed the number ", 
                  "of observations in Xr"))
    }
  }
  has_projection <- FALSE
  if (!is.matrix(diss_method)) {
    if (is.null(Xu)) {
      rdiss <- FALSE
    }
    else {
      rdiss <- control$return_dissimilarity
    }
    spike <- c(spike, -which(is.na(Yr)))
    if (length(spike) == 0) {
      spike <- NULL
    }
    neighborhoods <- get_neighbor_info(Xr = Xr, Xu = Xu, 
                                       diss_method = diss_method, Yr = Yr, k = k_max, k_diss = k_diss_max, 
                                       k_range = k_range, spike = spike, pc_selection = pc_selection, 
                                       return_dissimilarity = rdiss, center = center, scale = scale, 
                                       gh = gh, diss_usage = diss_usage, allow_parallel = control$allow_parallel, 
                                       ...)
    if (is.null(Xu)) {
      neighborhoods$neighbors <- rbind(k_0 = 1:ncol(neighborhoods$neighbors), 
                                       neighborhoods$neighbors)
      neighborhoods$neighbors_diss <- rbind(k_0 = 0, neighborhoods$neighbors_diss)
      k <- k + 1
      diss_xr_xtarget <- neighborhoods$diss_xr_xr
    }
    else {
      diss_xr_xtarget <- neighborhoods$dissimilarity
    }
    if (!is.null(neighborhoods$projection)) {
      diss_xr_xtarget_projection <- neighborhoods$projection
      has_projection <- TRUE
    }
  }
  else {
    diss_xr_xr <- NULL
    dim_diss <- dim(diss_method)
    if (diss_usage == "predictors") {
      if (diff(dim_diss) != 0 | dim_diss[1] != n_total | 
          any(diag(diss_method) != 0)) {
        stop(paste0("If a matrix is passed to 'diss_method' ", 
                    "and diss_usage = 'predictors', this matrix must be ", 
                    "squared symmetric zeroes in its diagonal"))
      }
      diss_xr_xr <- diss_method[1:nrow(Xr), 1:nrow(Xr)]
      if (!is.null(Xu)) {
        diss_method <- diss_method[1:nrow(Xr), (1 + nrow(Xr)):ncol(diss_method)]
      }
      else {
        diss_method <- diss_xr_xr
      }
    }
    if (!is.null(Xu)) {
      if (diss_usage %in% c("weights", "none")) {
        if (dim_diss[1] != n_xr & dim_diss[2] != n_xu) {
          stop(paste0("If a matrix is passed to 'diss_method' ", 
                      "and 'diss_usage' argument is set to either 'weights' or  ", 
                      "'none', the number of rows and columns of this matrix ", 
                      "must be equal to the number of rows of 'Xr' and the ", 
                      "number of rows of 'Xu' respectively"))
        }
      }
    }
    diss_xr_xtarget <- diss_method
    diss_method <- "external_matrix"
    neighborhoods <- diss_to_neighbors(diss_xr_xtarget, k = k_max, 
                                       k_diss = k_diss_max, k_range = k_range, spike = spike, 
                                       return_dissimilarity = control$return_dissimilarity, 
                                       skip_first = ifelse(is.null(Xu), TRUE, FALSE))
    if (is.null(Xu)) {
      neighborhoods$neighbors <- rbind(k_0 = 0, neighborhoods$neighbors)
      neighborhoods$neighbors_diss <- rbind(k_0 = 0, neighborhoods$neighbors_diss)
    }
    if (gh) {
      neighborhoods$gh$projection <- pls_projection(Xr = Xr, 
                                                    Xu = Xu, Yr = Yr, pc_selection = pc_selection, 
                                                    scale = scale, ...)
      neighborhoods$gh$gh_Xr <- f_diss(neighborhoods$gh$projection$scores, 
                                       Xu = t(colMeans(neighborhoods$gh$projection$scores[1:nrow(Xr), 
                                       ])), diss_method = "mahalanobis", center = FALSE, 
                                       scale = FALSE)
      if (!is.null(Xu)) {
        neighborhoods$gh$gh_Xu <- neighborhoods$gh$gh_Xr[-c(1:nrow(Xr))]
      }
      else {
        neighborhoods$gh$gh_Xu <- NULL
      }
      neighborhoods$gh$gh_Xr <- neighborhoods$gh$gh_Xr[c(1:nrow(Xr))]
      neighborhoods$gh <- neighborhoods$gh[c("gh_Xr", "gh_Xu", 
                                             "projection")]
    }
    neighborhoods$diss_xr_xr <- diss_xr_xr
    rm(diss_xr_xr)
    gc()
  }
  if (!is.null(k)) {
    smallest_neighborhood <- neighborhoods$neighbors[1:min(k), 
                                                     , drop = FALSE]
    smallest_n_neighbors <- colSums(!is.na(smallest_neighborhood))
  }
  if (!is.null(k_diss)) {
    min_diss <- neighborhoods$neighbors_diss <= min(k_diss)
    if (!is.null(spike)) {
      min_diss[1:length(spike), ] <- TRUE
    }
    smallest_neighborhood <- neighborhoods$neighbors
    smallest_neighborhood[!min_diss] <- NA
    smallest_n_neighbors <- colSums(!is.na(smallest_neighborhood))
    smallest_n_neighbors[smallest_n_neighbors < min(k_range)] <- min(k_range)
    smallest_n_neighbors[smallest_n_neighbors > max(k_range)] <- max(k_range)
  }
  if (is_local_cv) {
    min_n_samples <- floor(min(smallest_n_neighbors) * control$p) - 
      1
    min_cv_samples <- floor(min(k, k_range) * (1 - control$p))
    if (min_cv_samples < 3) {
      stop(paste0("Local cross-validation requires at least 3 observations in ", 
                  "the hold-out set, the current cross-validation parameters ", 
                  "leave less than 3 observations in some neighborhoods."))
    }
  }
  else {
    min_n_samples <- smallest_n_neighbors - 1
  }
  if (method$method %in% c("pls", "wapls")) {
    max_pls <- max(method$pls_c)
    if (any(min_n_samples < max_pls)) {
      stop(paste0("More pls components than observations in some neighborhoods.\n", 
                  "If 'local_cv' is being used, consider that some ", 
                  "observations \nin the neighborhoods are hold-out for local ", 
                  "validation"))
    }
  }
  if (!".local" %in% names(input_dots)) {
    iter_neighborhoods <- ith_mbl_neighbor(Xr = Xr, Xu = Xu, 
                                           Yr = Yr, Yu = Yu, diss_usage = diss_usage, neighbor_indices = neighborhoods$neighbors, 
                                           neighbor_diss = neighborhoods$neighbors_diss, diss_xr_xr = neighborhoods$diss_xr_xr, 
                                           group = group)
  }
  else {
    iter_neighborhoods <- ith_mbl_neighbor(Xr = Xr, Xu = Xu, 
                                           Yr = Yr, Yu = Yu, diss_usage = "none", neighbor_indices = neighborhoods$neighbors, 
                                           neighbor_diss = neighborhoods$neighbors_diss, group = group)
  }
  r_fields <- c("o_index", "k_diss", "k_original", "k", "npls", 
                "min_pls", "max_pls", y_output_name, y_hat_output_name, 
                "yr_min_obs", "yr_max_obs", "index_nearest_in_Xr", "index_farthest_in_Xr", 
                "y_nearest", "y_nearest_pred", "y_farthest", "diss_nearest", 
                "diss_farthest", "loc_rmse_cv", "loc_st_rmse_cv", "loc_n_components", 
                "rep")
  n_ith_result <- ifelse(is.null(k_diss), length(k), length(k_diss))
  template_pred_results <- data.table(matrix(NA, n_ith_result, 
                                             length(r_fields), dimnames = list(NULL, r_fields)))
  template_pred_results$rep[1] <- 0
  if (!is.null(k_diss)) {
    template_pred_results$k_diss <- k_diss
  }
  else {
    template_pred_results$k <- k
  }
  pg_bar_width <- 10
  if (!is.null(Xu)) {
    n_characters <- nchar(n_xu)
    n_iter <- n_xu
  }
  else {
    n_characters <- nchar(n_xr)
    n_iter <- n_xr
  }
  to_erase <- pg_bar_width + (2 * n_characters) + 8
  to_erase <- paste(rep(" ", to_erase), collapse = "")
  if (verbose) {
    cat("\033[32m\033[3mPredicting...\n\033[23m\033[39m")
  }
  pred_obs <- foreach(i = 1:n_iter, ith_observation = iter_neighborhoods, 
                      .inorder = FALSE, .export = c("ortho_diss", "fit_and_predict", 
                                                    "pls_cv", "get_col_sds", "get_wapls_weights"), .noexport = c("Xr", 
                                                                                                                 "Xu")) %mydo% {
                                                                                                                   ith_pred_results <- template_pred_results
                                                                                                                   additional_results <- NULL
                                                                                                                   ith_pred_results$o_index[] <- i
                                                                                                                   if (".local" %in% names(input_dots) & diss_method %in% 
                                                                                                                       ortho_diss_methods) {
                                                                                                                     ith_observation <- get_ith_local_neighbors(ith_xr = ith_observation$ith_xr, 
                                                                                                                                                                ith_xu = ith_observation$ith_xu, ith_yr = ith_observation$ith_yr, 
                                                                                                                                                                ith_yu = ith_observation$ith_yu, diss_usage = diss_usage, 
                                                                                                                                                                ith_neig_indices = ith_observation$ith_neig_indices, 
                                                                                                                                                                k = k_max, k_diss = k_diss_max, k_range = k_range, 
                                                                                                                                                                spike = spike, diss_method = diss_method, pc_selection = pc_selection, 
                                                                                                                                                                center = center, scale = scale, ith_group = ith_observation$ith_group, 
                                                                                                                                                                ...)
                                                                                                                     ith_pred_results$loc_n_components[] <- ith_observation$ith_components
                                                                                                                     additional_results$ith_neig_indices <- ith_observation$ith_neig_indices
                                                                                                                     additional_results$ith_neigh_diss <- ith_observation$ith_neigh_diss
                                                                                                                   }
                                                                                                                   if (verbose) {
                                                                                                                     cat(paste0("\033[34m\033[3m", i, "/", n_iter, "\033[23m\033[39m"))
                                                                                                                     pb <- txtProgressBar(width = pg_bar_width, char = "\033[34m_\033[39m")
                                                                                                                   }
                                                                                                                   if (!is.null(k_diss)) {
                                                                                                                     ith_diss <- ith_observation$ith_neigh_diss
                                                                                                                     if (!is.null(spike)) {
                                                                                                                       ith_diss[1:length(spike)] <- 0
                                                                                                                     }
                                                                                                                     ith_pred_results$k_original <- sapply(k_diss, FUN = function(x, 
                                                                                                                                                                                  d) sum(d < x), d = ith_diss)
                                                                                                                     ith_pred_results$k <- ith_pred_results$k_original
                                                                                                                     ith_pred_results$k[ith_pred_results$k_original < 
                                                                                                                                          min(k_range)] <- min(k_range)
                                                                                                                     ith_pred_results$k[ith_pred_results$k_original > 
                                                                                                                                          max(k_range)] <- max(k_range)
                                                                                                                   }
                                                                                                                   else {
                                                                                                                     ith_pred_results$k <- k
                                                                                                                   }
                                                                                                                   for (kk in 1:nrow(ith_pred_results)) {
                                                                                                                     if (verbose) {
                                                                                                                       setTxtProgressBar(pb, kk/nrow(ith_pred_results))
                                                                                                                     }
                                                                                                                     current_k <- ith_pred_results$k[kk]
                                                                                                                     if (current_k != ifelse(kk == 1, 0, ith_pred_results$k[kk - 
                                                                                                                                                                            1])) {
                                                                                                                       if (diss_usage == "predictors") {
                                                                                                                         keep_cols <- c(1:current_k, (1 + ith_observation$n_k):ncol(ith_observation$ith_xr))
                                                                                                                         i_k_xr <- ith_observation$ith_xr[1:current_k, 
                                                                                                                                                          keep_cols]
                                                                                                                         i_k_xu <- ith_observation$ith_xu[, keep_cols, 
                                                                                                                                                          drop = FALSE]
                                                                                                                       }
                                                                                                                       else {
                                                                                                                         i_k_xr <- ith_observation$ith_xr[1:current_k, 
                                                                                                                         ]
                                                                                                                         i_k_xu <- ith_observation$ith_xu
                                                                                                                       }
                                                                                                                       i_k_yr <- ith_observation$ith_yr[first_nn:current_k, 
                                                                                                                                                        , drop = FALSE]
                                                                                                                       i_k_yu <- ith_observation$ith_yu
                                                                                                                       kth_diss <- ith_observation$ith_neigh_diss[first_nn:current_k]
                                                                                                                       i_idx <- ith_observation$ith_neig_indices[first_nn:current_k]
                                                                                                                       ith_pred_results$rep[kk] <- 0
                                                                                                                       ith_yr_range <- range(i_k_yr)
                                                                                                                       ith_pred_results$yr_min_obs[kk] <- ith_yr_range[1]
                                                                                                                       ith_pred_results$yr_max_obs[kk] <- ith_yr_range[2]
                                                                                                                       ith_pred_results$diss_farthest[kk] <- max(kth_diss)
                                                                                                                       ith_pred_results$diss_nearest[kk] <- min(kth_diss)
                                                                                                                       ith_pred_results$y_farthest[kk] <- i_k_yr[which.max(kth_diss)]
                                                                                                                       ith_pred_results$y_nearest[kk] <- i_k_yr[which.min(kth_diss)]
                                                                                                                       ith_pred_results$index_nearest_in_Xr[kk] <- i_idx[which.min(kth_diss)]
                                                                                                                       ith_pred_results$index_farthest_in_Xr[kk] <- i_idx[which.max(kth_diss)]
                                                                                                                       i_k_yr <- ith_observation$ith_yr[1:current_k, 
                                                                                                                                                        , drop = FALSE]
                                                                                                                       if (!is.null(group)) {
                                                                                                                         i_k_group <- factor(ith_observation$ith_group[1:current_k])
                                                                                                                       }
                                                                                                                       else {
                                                                                                                         i_k_group <- NULL
                                                                                                                       }
                                                                                                                       if (diss_usage == "weights") {
                                                                                                                         if (is.null(Xu)) {
                                                                                                                           stop("'weights' are not yet supported for diss_usage")
                                                                                                                         }
                                                                                                                         std_kth_diss <- kth_diss/max(kth_diss)
                                                                                                                         kth_weights <- (1 - (std_kth_diss^3))^3
                                                                                                                         kth_weights[which(kth_weights == 0)] <- 1e-04
                                                                                                                       }
                                                                                                                       else {
                                                                                                                         kth_weights <- rep(1, current_k)
                                                                                                                       }
                                                                                                                       ith_observation$ith_neig_indices
                                                                                                                       i_k_pred <- fit_and_predict(x = i_k_xr, y = i_k_yr, 
                                                                                                                                                   pred_method = method$method, scale = scale, 
                                                                                                                                                   pls_c = method$pls_c, weights = kth_weights, 
                                                                                                                                                   newdata = i_k_xu, CV = is_local_cv, tune = control$tune_locally, 
                                                                                                                                                   group = i_k_group, p = control$p, number = control$number, 
                                                                                                                                                   noise_variance = method$noise_variance, range_prediction_limits = control$range_prediction_limits, 
                                                                                                                                                   pls_max_iter = 1, pls_tol = 1e-06, seed = seed, 
                                                                                                                                                   modified = ifelse(is.null(method$modified), 
                                                                                                                                                                     FALSE, method$modified))
                                                                                                                       i_k_pred$validation$models$coefficients
                                                                                                                       ith_pred_results[[y_hat_output_name]][kk] <- i_k_pred$prediction
                                                                                                                       selected_pls <- NULL
                                                                                                                       if (is_local_cv) {
                                                                                                                         if (control$tune_locally) {
                                                                                                                           best_row <- which.min(i_k_pred$validation$cv_results$rmse_cv)
                                                                                                                         }
                                                                                                                         else {
                                                                                                                           best_row <- ifelse(method$method == "pls", 
                                                                                                                                              method$pls_c, 1)
                                                                                                                         }
                                                                                                                         if (method$method == "pls") {
                                                                                                                           ith_pred_results$npls[kk] <- i_k_pred$validation$cv_results$npls[best_row]
                                                                                                                           selected_pls <- ith_pred_results$npls[kk]
                                                                                                                         }
                                                                                                                         if (method$method == "wapls") {
                                                                                                                           ith_pred_results$min_pls[kk] <- i_k_pred$validation$cv_results$min_component[best_row]
                                                                                                                           ith_pred_results$max_pls[kk] <- i_k_pred$validation$cv_results$max_component[best_row]
                                                                                                                           selected_pls <- i_k_pred$validation$cv_results[best_row, 
                                                                                                                                                                          1:2]
                                                                                                                         }
                                                                                                                         ith_pred_results$loc_rmse_cv[kk] <- i_k_pred$validation$cv_results$rmse_cv[best_row]
                                                                                                                         ith_pred_results$loc_st_rmse_cv[kk] <- i_k_pred$validation$cv_results$st_rmse_cv[best_row]
                                                                                                                       }
                                                                                                                       else {
                                                                                                                         if (method$method == "pls") {
                                                                                                                           ith_pred_results$npls[kk] <- method$pls_c
                                                                                                                           selected_pls <- ith_pred_results$npls[kk]
                                                                                                                         }
                                                                                                                         if (method$method == "wapls") {
                                                                                                                           ith_pred_results$min_pls[kk] <- method$pls_c[[1]]
                                                                                                                           ith_pred_results$max_pls[kk] <- method$pls_c[[2]]
                                                                                                                           selected_pls <- method$pls_c
                                                                                                                         }
                                                                                                                       }
                                                                                                                       if (is_nnv_val) {
                                                                                                                         if (!is.null(group)) {
                                                                                                                           out_group <- which(i_k_group == i_k_group[[ith_observation$local_index_nearest]])
                                                                                                                         }
                                                                                                                         else {
                                                                                                                           out_group <- ith_observation$local_index_nearest
                                                                                                                         }
                                                                                                                         nearest_pred <- fit_and_predict(x = i_k_xr[-out_group, 
                                                                                                                         ], y = i_k_yr[-out_group, , drop = FALSE], 
                                                                                                                         pred_method = method$method, scale = scale, 
                                                                                                                         pls_c = selected_pls, noise_variance = method$noise_variance, 
                                                                                                                         newdata = i_k_xr[ith_observation$local_index_nearest, 
                                                                                                                                          , drop = FALSE], CV = FALSE, tune = FALSE, 
                                                                                                                         range_prediction_limits = control$range_prediction_limits, 
                                                                                                                         pls_max_iter = 1, pls_tol = 1e-06, seed = seed, 
                                                                                                                         modified = ifelse(is.null(method$modified), 
                                                                                                                                           FALSE, method$modified))$prediction
                                                                                                                         ith_pred_results$y_nearest_pred[kk] <- nearest_pred/kth_weights[1]
                                                                                                                       }
                                                                                                                     }
                                                                                                                     else {
                                                                                                                       ith_k_diss <- ith_pred_results$k_diss[kk]
                                                                                                                       ith_pred_results[kk, ] <- ith_pred_results[kk - 
                                                                                                                                                                    1, ]
                                                                                                                       ith_pred_results$rep[kk] <- 1
                                                                                                                       ith_pred_results$k_diss[kk] <- ith_k_diss
                                                                                                                     }
                                                                                                                   }
                                                                                                                   if (verbose) {
                                                                                                                     if (kk == nrow(ith_pred_results) & i != n_iter) {
                                                                                                                       cat("\r", to_erase, "\r")
                                                                                                                     }
                                                                                                                     if (i == n_iter) {
                                                                                                                       cat("\n")
                                                                                                                     }
                                                                                                                   }
                                                                                                                   list(results = ith_pred_results, additional_results = additional_results)
                                                                                                                 }
  iteration_order <- sapply(pred_obs, FUN = function(x) x$results$o_index[1])
  pred_obs <- pred_obs[order(iteration_order, decreasing = FALSE)]
  results_table <- do.call("rbind", lapply(pred_obs, FUN = function(x) x$results))
  if (is.null(Xu) & !is.null(k)) {
    results_table$k <- results_table$k - 1
    fix_k <- 1
  }
  else {
    fix_k <- 0
  }
  if (".local" %in% names(input_dots) & diss_method %in% ortho_diss_methods) {
    diss_xr_xtarget <- do.call("cbind", lapply(iteration_order, 
                                               FUN = function(x, m, ii) {
                                                 idc <- x[[ii]]$additional_results$ith_neig_indices
                                                 d <- x[[ii]]$additional_results$ith_neigh_diss
                                                 m[idc] <- d
                                                 m
                                               }, x = pred_obs, m = matrix(NA, nrow(Xr), 1)))
    class(diss_xr_xtarget) <- c("local_ortho_diss", "matrix")
    dimnames(diss_xr_xtarget) <- list(paste0("Xr_", 1:nrow(diss_xr_xtarget)), 
                                      paste0(pre_nms_ng, 1:ncol(diss_xr_xtarget)))
    neighborhoods$neighbors <- do.call("cbind", lapply(iteration_order, 
                                                       FUN = function(x, m, ii) {
                                                         idc <- x[[ii]]$additional_results$ith_neig_indices
                                                         m[1:length(idc)] <- idc
                                                         m
                                                       }, x = pred_obs, m = matrix(NA, max(results_table$k), 
                                                                                   1)))
  }
  out <- c(if (is.null(Yu) & !is.null(Xu)) {
    "yu_obs"
  }, if (all(is.na(results_table$k_original))) {
    "k_original"
  }, if (!(validation_type %in% c("NNv", "both"))) {
    "y_nearest_pred"
  }, if (method$method != "wapls") {
    c("min_pls", "max_pls")
  }, if (method$method != "pls") {
    "npls"
  }, if (!(validation_type %in% c("local_cv", "both"))) {
    c("loc_rmse_cv", "loc_st_rmse_cv")
  }, "rep")
  results_table[, `:=`((out), NULL)]
  if (is.null(Xu)) {
    names(results_table)[names(results_table) %in% "yu_obs"] <- "yr_obs"
    names(results_table)[names(results_table) %in% "pred"] <- "fitted"
  }
  if (!is.null(k_diss)) {
    param <- "k_diss"
    results_table <- lapply(get(param), FUN = function(x, 
                                                       sel, i) x[x[[sel]] == i, ], x = results_table, sel = param)
    names(results_table) <- paste0("k_diss_", k_diss)
    p_bounded <- sapply(results_table, FUN = function(x, 
                                                      k_range) {
      sum(x$k_original <= k_range[1] | x$k_original >= 
            k_range[2])
    }, k_range = k_range)
    col_ks <- data.table(k_diss = k_diss, p_bounded = paste0(round(100 * 
                                                                     p_bounded/nrow(Xu), 3), "%"))
  }
  else {
    param <- "k"
    results_table <- lapply(get(param) - fix_k, FUN = function(x, 
                                                               sel, i) x[x[[sel]] == i, ], x = results_table, sel = param)
    names(results_table) <- paste0("k_", k - fix_k)
    col_ks <- data.table(k = k - fix_k)
  }
  if (validation_type %in% c("NNv", "both")) {
    nn_stats <- function(x) {
      nn_rmse <- (mean((x$y_nearest - x$y_nearest_pred)^2))^0.5
      nn_st_rmse <- nn_rmse/diff(range(x$y_nearest))
      nn_rsq <- (cor(x$y_nearest, x$y_nearest_pred))^2
      c(nn_rmse = nn_rmse, nn_st_rmse = nn_st_rmse, nn_rsq = nn_rsq)
    }
    loc_nn_res <- do.call("rbind", lapply(results_table, 
                                          FUN = nn_stats))
    loc_nn_res <- cbind(col_ks, rmse = loc_nn_res[, "nn_rmse"], 
                        st_rmse = loc_nn_res[, "nn_st_rmse"], r2 = loc_nn_res[, 
                                                                              "nn_rsq"])
  }
  else {
    loc_nn_res <- NULL
  }
  if (validation_type %in% c("local_cv", "both")) {
    mean_loc_res <- function(x) {
      mean_loc_rmse <- mean(x$loc_rmse_cv)
      mean_loc_st_rmse <- mean(x$loc_st_rmse_cv)
      c(loc_rmse = mean_loc_rmse, loc_st_rmse = mean_loc_st_rmse)
    }
    loc_res <- do.call("rbind", lapply(results_table, mean_loc_res))
    loc_res <- cbind(col_ks, rmse = loc_res[, "loc_rmse"], 
                     st_rmse = loc_res[, "loc_st_rmse"])
  }
  else {
    loc_res <- NULL
  }
  if (!is.null(observed)) {
    for (i in 1:length(results_table)) {
      results_table[[i]][[y_output_name]] <- observed
    }
    yu_stats <- function(x, y_hat, y) {
      y_rmse <- mean((x[[y_hat]] - x[[y]])^2, na.rm = TRUE)^0.5
      y_st_rmse <- y_rmse/diff(range(x[[y_hat]]), na.rm = TRUE)
      y_rsq <- cor(x[[y_hat]], x[[y]], use = "complete.obs")^2
      c(rmse = y_rmse, st_rmse = y_st_rmse, rsq = y_rsq)
    }
    pred_res <- do.call("rbind", lapply(results_table, yu_stats, 
                                        y_hat = y_hat_output_name, y = y_output_name))
    pred_res <- cbind(col_ks, rmse = pred_res[, "rmse"], 
                      st_rmse = pred_res[, "st_rmse"], r2 = pred_res[, 
                                                                     "rsq"])
  }
  else {
    pred_res <- NULL
  }
  if ("local_ortho_diss" %in% class(diss_xr_xtarget)) {
    diss_method <- paste0(diss_method, " (locally computed)")
  }
  if (control$return_dissimilarity) {
    diss_list <- list(diss_method = diss_method, diss_xr_xu = diss_xr_xtarget)
    if (has_projection) {
      diss_list$global_projection <- diss_xr_xtarget_projection
    }
  }
  else {
    diss_list <- NULL
  }
  colnames(neighborhoods$neighbors) <- paste0(pre_nms_ng, 1:ln)
  rownames(neighborhoods$neighbors) <- paste0("k_", 1:nrow(neighborhoods$neighbors))
  val_list <- structure(list(loc_res, loc_nn_res, pred_res), 
                        names = c("local_cross_validation", "nearest_neighbor_validation", 
                                  val_summary_name))
  results_list <- list(call = f_call, cntrl_param = control, 
                       dissimilarities = diss_list, Xu_neighbors = list(neighbors = neighborhoods$neighbors, 
                                                                        neighbors_diss = neighborhoods$neighbors_diss), n_predictions = nrow(Xu), 
                       gh = neighborhoods$gh, validation_results = val_list, 
                       results = results_table, documentation = documentation, 
                       seed = seed)
  if (is.null(Xu)) {
    names(results_list)[names(results_list) %in% "Xu_neighbors"] <- "Xr_neighbors"
  }
  attr(results_list, "call") <- f_call
  class(results_list) <- c("mbl", "list")
  results_list
}







