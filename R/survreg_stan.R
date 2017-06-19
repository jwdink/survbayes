
#' Map the name of a PDF in Stan to its name in R
#'
#' @param stan_pdf Chraracter-string with stan name
#'
#' @return Character-string with R function name
#' @export
stan_to_r_pdf <- function(stan_pdf) {
  mapper <- c(`normal` = 'dnorm', `double_exponential` = 'dlaplace')
  mapper[stan_pdf]
}

#' #' Prediction method for \code{survreg_stan}
#' #'
#' #' @param object Of type \code{survreg_stan}
#' #' @param newdata Data.frame
#' #' @param type Either 'terms' (parameter estimates), or 'survival'/'hazard' (at 'times')
#' #' @param times A vector with length = nrow(newdata), specifying the times at which to estimate
#' #'   survival/hazard
#' #' @param time_grid Instead of specifying a vector of times with length = nrow(newdata), you can
#' #'   specify a grid of times which, for each row of the dataset, you'd like to get estimates. So
#' #'   \code{time_grid = 1:10} means that, for each row of the dataset, you'll get survival/hazard
#' #'   estimates for times 1 through 10.
#' #' @param start Start/truncation times, must be used with \code{times}, not \code{time_grid}
#' #' @param num_samples How many samples from the posterior should be used in getting
#' #'   parameter-estimates? To avoid running out of memory, you might want this to be less than the
#' #'   total number of samples. Default is 500.
#' #'
#' #' @return A vector or matrix, depending on the options above.
#' #' @export
#' predict.survreg_stan <- function(object, newdata = NULL, type = 'terms',
#'                                  times = NULL, start = NULL,
#'                                  time_grid = NULL,
#'                                  num_samples = 500) {
#'   type <- match.arg(type, c('terms','survival','hazard'))
#'   #summary <- match.arg(summary, c('mean_terms','mean_ci','none','mean'))
#'
#'   if (is.null(newdata)) newdata <- object$data
#'
#'   dist_info <- object$dist_info
#'
#'   ## form model matrices:
#'   #cat("\nMaking model mats...")
#'   formula_concat <- formula(object$concat_terms)
#'   forms <- purrr::map(object$form_terms, formula)
#'   lazyeval::f_lhs(formula_concat) <- NULL
#'   model_frame <- model.frame(formula = formula_concat, data = newdata, na.action = na.pass, xlev = object$xlevels)
#'   if (!is.null(object$scale_params)) {
#'     model_frame[,names(object$scale_params)] <-
#'       purrr::map2(.x = model_frame[,names(object$scale_params),drop=FALSE],
#'                   .y = object$scale_params[names(object$scale_params)],
#'                   .f = ~ (.x-.y$mean)/.y$sd )
#'   }
#'   model_mats <- purrr::map(forms, ~model.matrix(.x, model_frame) )
#'
#'   ## reduce number of samples if requested:
#'   num_samples_all <- nrow(object$posterior_mats[[1]])
#'   if (is.null(num_samples)) {
#'     num_samples <- num_samples_all
#'   } else {
#'     if (num_samples_all>num_samples) {
#'       object$posterior_mats <- purrr::map(object$posterior_mats, ~ .x[seq_len(num_samples),,drop=FALSE])
#'     } else{
#'       num_samples <- num_samples_all
#'     }
#'   }
#'
#'   ##
#'   #cat("\nPreallocating estimate-mats...")
#'   list_of_estimate_mats <- vector(mode = 'list', length = length(dist_info$pars))
#'   for (i in seq_len(length(dist_info$pars)))
#'     list_of_estimate_mats[[i]] <- bigmemory::big.matrix(nrow = nrow(newdata),
#'                                                         ncol = nrow(object$posterior_mats[[i]]))
#'   names(list_of_estimate_mats) <- dist_info$pars
#'
#'   #
#'   #cat("\nDoing matrix multiplciation...")
#'   object$posterior_mats <- purrr::map(.x = object$posterior_mats, as.matrix)
#'   transform_funcs <- purrr::map(dist_info$transforms, get, mode='function')
#'   inv_transform_funcs <- purrr::map(dist_info$inverse_transforms, get, mode='function')
#'   for (param_name in dist_info$pars) {
#'     this_beta_mat <- transform_funcs[[param_name]]( object$posterior_mats[[param_name]] )
#'     list_of_estimate_mats[[param_name]] <-
#'       inv_transform_funcs[[param_name]](model_mats[[param_name]] %*% t(this_beta_mat[,colnames(model_mats[[param_name]]),drop=FALSE]))
#'   }
#'
#'   HPD <- function(x) {
#'     if (length(na.omit(x))<10)
#'       tibble::data_frame(mean=mean(x,na.rm=TRUE), conf.lower=NA, conf.upper=NA)
#'     else
#'       structure(
#'         cbind(mean(x),coda::HPDinterval(coda::as.mcmc(x))) %>% tibble::as_data_frame()
#'         , names =  c('mean','conf.lower','conf.upper')
#'       )
#'   }
#'
#'   if (type == 'terms') {
#'
#'     stop(call. = FALSE, "Please report this error to the package maintainer.")
#'
#'   } else if (type %in% c('survival','hazard')) {
#'
#'     if (!is.null(time_grid)) {
#'       stop("to do")
#'       param_args <- data.frame(q = rep(time_grid, times = nrow(newdata)) )
#'       for (param in names(list_of_estimate_mats)) {
#'         temp <- rowMeans(list_of_estimate_mats[[param]], na.rm = TRUE)
#'         param_args[[param]] <- rep(temp, each = length(time_grid))
#'       }
#'     } else if (!is.null(times)) {
#'       param_args <- data.frame(q = times)
#'       for (param in names(list_of_estimate_mats))
#'         param_args[[param]] <- rowMeans(list_of_estimate_mats[[param]], na.rm = TRUE)
#'     } else {
#'       stop(call. = FALSE, "For `type != 'terms'`, please specify `times` or `time_grid`.")
#'     }
#'     if (!is.null(start))
#'       stop(call. = FALSE, "Please report this error to the package maintainer.")#newdata$..start <- s
#'
#'     df_tidy_posterior <- param_args
#'     cc <- which(complete.cases(param_args))
#'     if (is.null(object$all_knots)) {
#'       param_args <- as.list(param_args[cc,])
#'     } else {
#'       param_args <- list(q = param_args$q[cc], gamma = as.matrix(param_args[cc,-1,drop=FALSE]))
#'     }
#'
#'     #cat("\ndo.call...")
#'     df_tidy_posterior$result <- NA
#'     df_tidy_posterior$result[cc] <- do.call(dist_info$cdf_function, c(list(lower.tail=FALSE), param_args) )
#'
#'     if (type == 'hazard') {
#'       #cat("\ndo.call again...")
#'       param_args_lagged <- param_args
#'       param_args_lagged$q <- param_args$q - 1
#'       warning(call. = FALSE, "Using approximation for hazard-- please report this warning to the package maintainer.")
#'       df_tidy_posterior$result_lagged <- NA
#'       df_tidy_posterior$result_lagged[cc] <- do.call(dist_info$cdf_function, c(list(lower.tail=FALSE), param_args_lagged) )
#'       df_tidy_posterior$result <- with(df_tidy_posterior, 1 - result/result_lagged)
#'       df_tidy_posterior$result_lagged <- NULL
#'     }
#'
#'     if (!is.null(time_grid))
#'       stop("to do")
#'     else
#'       df_tidy_posterior$result
#'
#'   }
#'
#' }



# make_stancode <- function(survreg_data, save_model = NULL) {
#   string <- survreg_data$dist_info$model_code_function(model_mats = survreg_data$model_mats, priors = survreg_data$priors)
#   if (is.null(save_model))
#     string
#   else
#     cat(string, file = save_model)
# }
#
# convert_to_survstan <- function(stanfit, survreg_data) {
#   browser()
# }





