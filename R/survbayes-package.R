#' survbayes package
#'
#' @import lazyeval
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import survival
#'
#' @importFrom tibble rownames_to_column lst
#' @importFrom purrr map map2 pmap map_dbl map_chr map_lgl
#'             map2_chr map2_dbl map2_lgl transpose flatten_chr flatten_dbl flatten_lgl
#'             walk walk2 map_df map2_df
#' @importFrom tidysurv get_terms_mapper
#' @importFrom graphics plot
#' @importFrom stats sd terms update integrate delete.response
#'             as.formula coef fitted glm lm median na.omit
#'             predict var quantile model.response model.frame
#'             na.pass na.exclude na.fail model.matrix model.weights
#'             .getXlevels .checkMFClasses reformulate logLik model.extract optim
#'             contr.treatment formula setNames
#' @importFrom grDevices dev.off pdf
#' @importFrom glue glue
#'
#' @docType package
#' @name survbayes
NULL
#> NULL


# SurvReg Data --------------------------------------------------------------------------------

#' Contrasts where each factor-level is explicitly represented
#'
#' This contrast-type should only be used in \code{survreg_map}, where (a) the function knows to
#' center the contrasts so that each coefficient estimate can be interpreted as "deviation from
#' average," and (b) the fact that this contrast-type is linearly dependent isn't a problem, thanks
#' to shrinkage from priors.
#'
#' @param n a vector of levels for a factor, or the number of levels.
#' @param contrasts Just set this to TRUE
#' @param sparse logical indicating if the result should be sparse
#'
#' @return Contrast-matrix
#' @export
contr.full <- function(n, contrasts=TRUE, sparse=FALSE) {
  stop(call. = FALSE,
       "`contr.full` is currently broken, please report this error to the package maintainer.")
  if (!contrasts)
    stop(call. = FALSE, 'Please report this error to the package maintainer.')
  out <- contr.treatment(n = n, contrasts = contrasts, sparse = sparse, base = 1)
  to_add <- matrix(nrow = nrow(out), ncol = 1, data = c(1,rep(0,nrow(out)-1)),
                   dimnames = list(row.names(out), rownames(out)[1]))
  cbind(out,to_add)
}

#' Prepare components of a model (model.frame, response, etc) for bayesian survival-regression
#'
#' @param formula_concat A formula containing all the terms for all the parameters
#' @param forms A list of formula, one for each parameter
#' @param data The data.frame
#' @param na.action Function for NAs
#' @param xlev Passed to model.frame
#' @param standardize_y If TRUE, then response variable is standardized. This means taking the log
#'   of the 'end' times, standardizing (centering/scaling) them, then taking the exponent of that.
#'   If this is a list, it's assumed to have components 'center' and 'scale', for performing this
#'   standardization.
#'
#' @return A list with needed components, such as model.frame/mats, response object, and info for
#'   subsequent calls (params for standardizing).
prep_model_components <- function(formula_concat, forms, data, na.action, xlev, standardize_y) {

  model_frame <- model.frame(formula = formula_concat, data = data, na.action = na.action, drop.unused.levels = FALSE, xlev = xlev)
  xlevels <- .getXlevels(attr(model_frame, "terms"), model_frame)
  if (!is.null(xlevels)) xlevels <- xlevels[map_lgl(names(xlevels), ~is.character(data[[.x]]))]
  if (length(xlevels)==0) xlevels <- NULL
  form_terms <- purrr::map(forms, ~delete.response(terms(.x, data=data)) )
  model_mats <- purrr::map(form_terms, model.matrix,  data=model_frame)
  list_of_term_mappers <- map(forms, tidysurv::get_terms_mapper, data = data)

  terms_mapper <- tidysurv::get_terms_mapper(formula = formula_concat, data = data)

  scale_params <- map(seq_along(model_mats), function(i) {
    mm <- model_mats[[i]]
    the_terms <- form_terms[[i]]
    if (ncol(mm)==1 && colnames(mm)=="(Intercept)") {
      out <- data_frame(term = "(Intercept)",
                        center = 0,
                        scale = 1,
                        undo_for_coefs =FALSE)
      return(out)
    } else {
      # for each model-mat term, get all the cols in the original data that went into it:
      fact_mat <- attr(the_terms,'factors')
      mapping_from_mm_col_to_data_col <-
        map(attr(mm,'assign'), function(assign_idx) row.names(fact_mat)[1==fact_mat[,assign_idx,drop=TRUE]])
      # only consider model-mat term to be numeric if all the cols that went into it were numeric:
      is_numeric_term <- map_lgl(mapping_from_mm_col_to_data_col, function(col)
        if (length(col)>0 && col %in% colnames(model_frame)) all(map_lgl(as.list(model_frame[,col,drop=FALSE]),is.numeric)) else FALSE)
      names(is_numeric_term) <- colnames(mm)

      cols_to_center <- names(is_numeric_term)[is_numeric_term]
      cols_to_scale <- names(is_numeric_term)[is_numeric_term]

      # additionally, we want to center any factors with the 'contr.full' contrast:
      contrasts_attr <- attr(mm, 'contrasts')
      if (is.null(contrasts_attr)) contrasts_attr <- list()
      is_full <- vector(mode = 'logical', length = length(contrasts_attr))
      is_mat_lgl <- map_lgl(contrasts_attr, is.matrix)
      is_full[is_mat_lgl] <- map_lgl(contrasts_attr[is_mat_lgl], ~ncol(.x)==nrow(.x))
      is_full[!is_mat_lgl] <- map_lgl(contrasts_attr[!is_mat_lgl], ~.x=='contr.full')
      full_contrast_mf_cols <- names(contrasts_attr)[is_full]
      if (!is.null(full_contrast_mf_cols) && length(full_contrast_mf_cols)>0) {
        full_contrast_od_cols <- flatten_chr(list_of_term_mappers[[i]](model_frame_cols=full_contrast_mf_cols))
        cols_to_center <- c(cols_to_center,
                            flatten_chr(list_of_term_mappers[[i]](original_cols=full_contrast_od_cols)))
      }

      centers <- setNames(rep(NA, ncol(mm)), colnames(mm))
      centers[cols_to_center] <- map_dbl(as.data.frame(mm[,cols_to_center]), mean, na.rm=TRUE)
      scales <- setNames(rep(NA, ncol(mm)), colnames(mm))
      scales[cols_to_scale] <- map_dbl(as.data.frame(mm[,cols_to_scale]), sd, na.rm=TRUE)
      undo_for_coefs <- names(centers) %in% names(is_numeric_term)[is_numeric_term]
      data_frame(term = colnames(mm),
                 center = centers,
                 scale = scales,
                 undo_for_coefs =undo_for_coefs)
    }
  })

  y <- model.extract(model_frame, "response")
  if (!is.null(y)) {
    if (class(y)[1] == 'Surv') {
      # need to add support for non-interval censoring. will need dbasis
      stop(call. = FALSE, "Please report this error to the package maintainer.")
    } else if (class(y)[1] == 'Survint') {
      if (is.logical(standardize_y) && standardize_y)
        standardize_y <- list(center = mean(log(y[,3,drop=TRUE]), na.rm=TRUE),
                              scale = sd(log(y[,3,drop=TRUE]), na.rm = TRUE) )
      if (is.list(standardize_y)) {
        y[,-4] <- exp(scale(x = log(y[,-4,drop=FALSE]), center = rep(standardize_y$center,3), scale = rep(standardize_y$scale,3)))
      } else {
        standardize_y <- list(center=0, scale=1)
      }
    } else {
      stop(call. = FALSE, "Type of survival-object on formula's left-hand-side not recognized.")
    }
  } else {
    standardize_y <- NULL
  }

  list(model_frame = model_frame,
       xlevels = xlevels,
       model_mats = model_mats,
       scale_params = scale_params,
       y = y, standardize_y = standardize_y)
}

#' Prepare data for a bayesian survival-regression model
#'
#' @param formula Formula with both rhs and lhs. Rhs will be applied to 'location' parameter of the distribution, or for spline-based models, the first parameter (gamma0).
#' @param anc A list of lhs-only formulae. Must be named for the parameters of the distribution (e.g., \code{list(gamma1 = ~ pred1 + pred2)}).
#' @param data A data.frame
#' @param distribution Character naming the distribution
#' @param dist_config A list customizing the distribution. For example, if \code{distribution = 'roy_parm_splines'}, then you can pass \code{list(knots=)} to specify the location of the knots.
#' @param standardize_y If TRUE (the default), then response variable is standardized. This means taking the log
#'   of the 'end' times, standardizing (centering/scaling) them, then taking the exponent of that.
#'   If this is a list, it's assumed to have components 'center' and 'scale', for performing this
#'   standardization. Standardization is helpful because it makes the default priors more meaningful
#'   across datasets.
#' @param na.action Function for NAs
#'
#' @return A list with components that can be passed to \code{rstan::stan}, or the list can be passed to \code{survreg_map}.
#' @export
prep_survreg_data <- function(formula,
                              anc = NULL,
                              data,
                              distribution = 'roy_parm_splines',
                              dist_config = list(knots=NULL),
                              standardize_y = TRUE,
                              na.action = na.exclude) {
  stopifnot(is.data.frame(data))

  ## set ancilliary formula:
  if (length(formula)!=3) stop(call. = FALSE, "`formula` must have rhs and lhs.")
  if (!is.null(anc)) {
    anc_arg <- anc
    if (!is.list(anc_arg) || is.null(names(anc_arg)) || any(names(anc_arg)=="") )
      stop("`anc` must be a named list")
    if (!all(purrr::map_lgl(anc_arg, ~inherits(.x, "formula"))))
      stop("`anc` must be a list of formulae")

    dist_info <- get_dist_info(distribution, k = length(anc) )

    non_loc_pars <- setdiff(dist_info$pars, dist_info$location)
    anc <- structure(names = non_loc_pars, purrr::rerun(length(non_loc_pars), ~1))
    anc[names(anc_arg)] <- anc_arg
    unmatched_params <- setdiff( names(anc_arg), dist_info$pars)
    if (length(unmatched_params)>0)
      stop("The following names in your `anc` argument are not parameters of '", distribution, "':", paste0(", ", unmatched_params))
  } else {
    dist_info <- get_dist_info(distribution, k = 1)
    non_loc_pars <- setdiff(dist_info$pars, dist_info$location)
    anc <- structure(names = non_loc_pars, purrr::rerun(length(non_loc_pars), ~1))
  }

  forms <- c(location=formula, anc)
  names(forms)[[1]] <- dist_info$location

  if (distribution=='roy_parm_splines') {
    if (is.null(dist_config$knots)) {
      if (length(anc)>1) stop("For this distribution, please specify `knots`).")
      else dist_config$knots <- numeric(0)
    } else {
      stopifnot( length(dist_config$knots)==(length(anc)-1) )
    }
  }

  ##
  forms <- purrr::map(forms, function(.x) {
    lazyeval::f_lhs(.x) <- lazyeval::f_lhs(formula)
    return(.x)# add lhs to each formula, so that . will be interpreted correctly
  })
  formula_concat <- concatenate_formula(formula,forms,data)

  model_components <- prep_model_components(formula_concat, forms, data, na.action, xlev = NULL, standardize_y=standardize_y)
  names(model_components$model_mats) <- dist_info$pars_real
  names(model_components$scale_params) <- dist_info$pars_real

  ## boundary knots:
  if (!is.null(dist_config$knots)) {
    inner_knots <- (log(dist_config$knots)-model_components$standardize_y$center)/model_components$standardize_y$scale
    # TODO: user-specified boundary knots
    all_knots <- c(log(min(model_components$y[,3], na.rm=TRUE)), inner_knots, log(max(model_components$y[,3], na.rm=TRUE)))
    if (all_knots[1]==-Inf) {
      all_knots[1] <- (-30) # a very small number
      warning(call. = FALSE, "At least one `event` occurs at `time=0`.")
    }
    if ( any(all_knots == -Inf, na.rm = TRUE) | any(is.na(all_knots)) )
      stop(call. = FALSE,
           "When taking log of the knots, NaNs or -Infs resulted. Did you place the knots on the unit of log-time, instead of on the unit of time?")
  }

  ## priors
  priors <- dist_info$default_priors(model_components$model_mats, model_components$y)
  df_scale_params <- bind_rows(model_components$scale_params, .id = 'parameter')
  priors <- left_join(x = priors, y = df_scale_params, by = c('parameter','term'))
  inits <- purrr::map(split(priors, priors$parameter), ~with(.x, structure(mu,names=term)))
  inits <- purrr::map(inits, ~array(.x, dim = length(.x)))

  ## to be passed as the 'data' arg to `stan`.
  stan_data <- list(Y = model_components$y,
                    N = nrow(model_components$y),
                    obs_idx = which(model_components$y[,4]==1),
                    N_obs = sum(model_components$y[,4]==1),
                    cens_idx = which(model_components$y[,4]==0),
                    N_cens = sum(model_components$y[,4]==0),
                    scaling_factor = 50 ) # <---- better system for this
  if (!is.null(dist_config$knots)) {
    stan_data$knots <- all_knots
    stan_data$k <- length(all_knots)
  }
  for (i in seq_along(model_components$model_mats)) {
    stan_data[[paste0('num_covs_gamma',i-1)]] <- ncol(model_components$model_mats[[i]])
    stan_data[[paste0('X_gamma',i-1)]] <- scale(model_components$model_mats[[i]],
                                                center = coalesce(model_components$scale_params[[i]]$center,0),
                                                scale = coalesce(model_components$scale_params[[i]]$scale,1))
  }

  ## out:
  out <- list(data = data,
              stan_data = stan_data,
              formula_concat = formula_concat,
              scale_params = model_components$scale_params,
              log_time_scaling = model_components$standardize_y,
              forms = forms,
              priors = priors,
              inits = inits,
              na.action = na.action,
              xlevels = model_components$xlevels,
              dist_info = dist_info)
  if (!is.null(dist_config$knots))
    out$all_knots <- all_knots
  class(out) <- 'survreg_data'
  return(out)

}

#' @export
print.survreg_data <- function(x, title=TRUE, ...) {
  if (title)
  cat(sep = "",
      "A Dataset ready for parametric survival-regression with the '",x$dist_info$name, "' family.\n\n")

  cat("Formulae: ====\n")
  cat("\n$DV\n")
  print(x$formula_concat[[2]])
  if ( near(x$log_time_scaling$center,0) & near(x$log_time_scaling$scale,1) )
    cat("DV not standardized")
  else
    cat(glue::glue("Standardized so that `standardized(time) = exp((log(time)-{x$log_time_scaling$center})/{x$log_time_scaling$scale}))`"))
  cat("\n\n")
  print(map(map(map(x$forms, terms),delete.response),formula))

  cat("Priors: ====\n")
  cat("(Priors applied to standardized (centered/scaled) variables)\n")
  print(select(x$priors,-undo_for_coefs))
}


# General Helpers -------------------------------------------------------------------------------------

#' Concatenate a list of formula
#'
#' @param formula Main formula
#' @param forms List of rhs-only formula, from `anc` argument.
#' @param data Data.frame
#'
#' @return Formula
concatenate_formula <- function(formula, forms, data) {
  # TO DO: do i really need/want this function?
  term_objs <- purrr::map(forms, ~terms(.x, data=data))
  covnames <- unique(purrr::flatten_chr( purrr::map(term_objs, ~attr(.x,'term.labels')) ))
  cov_formula_char <- if (length(covnames)==0) "1" else paste0(collapse = " + ", covnames)
  concat_formula_char <- paste0(paste0(deparse(lazyeval::f_lhs(formula)),collapse=""), " ~ ", cov_formula_char)
  concat_formula <- as.formula(concat_formula_char)
  environment(concat_formula) <- environment(formula)
  covnames_bare <- unique(purrr::flatten_chr( purrr::map(term_objs, ~all.vars(delete.response(.x)) ) ))
  attr(concat_formula, 'covnames') <- covnames_bare
  attr(concat_formula, 'covnames.orig') <- covnames
  concat_formula
}

#' Plot model coefficients
#'
#' @param object A model object
#' @param ... Other arguments to be passed to methods
#'
#' @return A ggplot object
#' @export
plot_coefs <- function(object, ...) {
  UseMethod('plot_coefs')
}

#' Get information about a survival distribution for survreg_stan
#'
#' @param dist Character naming a distribution
#' @param k For spline-based distributions, the number of knots
#'
#' @return A list with information, functions, etc.
#' @export
get_dist_info <- function(dist, k = NULL) {
  dist_infos <- list()

  dist <- match.arg(arg = dist, choices = c('roy_parm_splines'))
  if (dist=='roy_parm_splines' & is.null(k))
    stop(call. = FALSE, "Please specify 'k' (number of knots) for this distribution.")

  if (!is.null(k)) {
    dist_infos[['roy_parm_splines']] <- tibble::lst(
      name = 'Royston-Parmesian Cumulative-Splines',
      location = 'gamma0',
      pars = paste0('gamma', 0:k),
      transforms = as.list(structure(names=pars, c('identity', 'log', rep('identity',k-1)))),
      inverse_transforms = as.list(structure(names=pars, c('identity', 'exp', rep('identity',k-1)))),
      pars_real = paste0(ifelse(transforms=="identity","",paste0(transforms,"_")),pars),

      default_priors = function(model_mats, y = NULL) {
        param_names <- names(model_mats)

        if (is.null(y)) gamma0_init <- 0
        else gamma0_init <- -log(mean(y[,3]))

        priors <- purrr::map_df(
          model_mats, .id = 'parameter',
          function(mat) {
            data_frame(term = colnames(mat),
                       mu = rep(0, ncol(mat)),
                       sigma = rep(1, ncol(mat)) )
          })
        priors$mu <- ifelse(priors$parameter == 'gamma0' & priors$term == "(Intercept)",
                                  yes = gamma0_init,
                                  no = priors$mu)
        priors
      },

      distribution_function_factory = function(all_knots, scaling_factor) {

        cdf_function <- function(q, gamma, lower.tail = TRUE, log.p =FALSE, b = NULL) {
          invalid_idx <- which(q<=0)
          q[invalid_idx] <- NA_real_
          log_cumu_haz <- roy_parm_log_cumu_haz(log_time = log(q), gamma = gamma, all_knots = all_knots, scaling_factor = scaling_factor, b = b)
          log_cumu_haz[invalid_idx] <- -Inf
          cumu_haz <- exp(log_cumu_haz)
          surv <- exp(-cumu_haz)
          if (lower.tail) out <- 1-surv
          else out <- surv
          if (log.p) return(log(out))
          else return(out)
        }

        list(pdf_function = function() stop(call. = FALSE, "Please report this error to the package maintainer."), # TO DO!
             cdf_function = cdf_function
        )
      }
      ,
      model_code_function = make_roy_parm_spline_stan_model_code

    )
  }



  return(dist_infos[[dist]])
}

#' Survival object allowing for data that is both interval-censored *and* truncated
#'
#' @param start   Truncation time. If unspecified
#' @param end     Last observed time.
#' @param event   Event indicator. When `event==0`, `end_lb` is ignored, and `end` is considered
#'   right-censored. When `event==1` and `end > end_lb`, the event is considered interval-censored
#'   (occurring sometime between these two times); when `event==1` and `end == end_lb`, the `end` is
#'   considered an exact time of event.
#' @param end_lb  Optional. Lower bound of last-observed time. If included, then an `event`
#'   indicator signals inteval-censoring (i.e., an event some time between `end_lb` and `end`). If
#'   some observations are interval-censored, while for some the exact time is known, the latter can
#'   be indicated by rows where `end == end_lb`. If the `end_lb` argument is omitted entirely, then
#'   `end == end_lb` is assumed for all data.
#'
#' @return An object of class 'Survint', which can be used as the response-object in
#'   \code{survreg_stan} or \code{survreg_map}.
#' @export
Survint <- function(end, event, start = NULL, end_lb = NULL) {
  if (is.null(start)) start <- numeric(length(end))
  stopifnot(is.numeric(start))
  stopifnot(is.numeric(end))
  if (is.null(end_lb)) end_lb <- end
  stopifnot(is.numeric(end_lb))
  stopifnot(is.numeric(event) & all(dplyr::between(event, 0, 1), na.rm = TRUE) )
  out <- tibble::data_frame(start = start, end_lb = end_lb, end = end, event = event)
  structure(as.matrix(out), class = c('Survint','matrix'))
}



# SurvReg MAP ---------------------------------------------------------------------------------

#' @export
terms.survreg_map <- function(x, ...) {
  terms(x$formula_concat)
}

#' Predict method for 'survreg_map'
#'
#' @param object 'survreg_map' object
#' @param newdata A data.frame, optional
#' @param times Only required if type != 'parameters. A vector of times with length 1 or nrow(newdata)
#' @param starts Optional, for type = 'survival' only. A vector of start/truncation times with length 1 or nrow(newdata). Survival times will be conditional up to survival at this point.
#' @param type Type of prediction. Can be predicted 'parameters' of the distribution for each row, or 'survival' probabilities (at 'times', given 'starts').
#' @param na.action Function for dealing with NAs
#' @param ... Ignored
#'
#' @return A matrix if type = 'parameters', or a vector if type = 'survival'.
#' @export
predict.survreg_map <- function(object, newdata = NULL, times, starts = NULL, type = 'survival', na.action = na.pass, ...) {
  if (is.null(newdata)) newdata <- object$data

  type <- match.arg(arg = type, choices = c('survival','parameters'))

  object$formula_concat[[2]]<- NULL
  model_components <- prep_model_components(formula_concat = object$formula_concat,
                                            forms = object$forms,
                                            data = newdata,
                                            na.action = na.action,
                                            xlev = object$xlevels,
                                            standardize_y = object$log_time_scaling)
  names(model_components$model_mats) <- object$dist_info$pars_real

  model_mats_std <- map(seq_along(model_components$model_mats),
                        function(i)
                          scale(model_components$model_mats[[i]],
                                center = coalesce(object$scale_params[[i]]$center,0),
                                scale = coalesce(object$scale_params[[i]]$scale,1) ))
  names(model_mats_std) <- object$dist_info$pars_real

  unrolled_par <- with(object$res_std, structure( estimate, names = paste0(parameter, "__sep__", term)))
  dist_params <- object$helper_funs$get_distribution_params_from_unrolled_coefs(unrolled_par,
                                                                                list_of_standardized_model_mats = model_mats_std)
  dist_params <- as.matrix(dist_params)

  if (length(times)==1)
    times <- rep(x=times, times=nrow(newdata))
  if (is.null(starts))
    starts <- rep(0, length(times))
  times <- with(object$log_time_scaling, exp( (log(times) - center)/scale ) )
  starts <- with(object$log_time_scaling, exp( (log(starts) - center)/scale ) )

  na_dropped_idx <- attr(model_components$model_frame,'na.action')
  if (length(na_dropped_idx)>0) {
    times <- times[-na_dropped_idx]
    starts <- starts[-na_dropped_idx]
  }

  if (type == 'parameters')
    return(dist_params)

  if (!is.null(object$all_knots)) {
    if (type == 'survival')
      surv <- object$helper_funs$distribution_functions$cdf_function(q = times, gamma = dist_params, lower.tail =FALSE)
    if (any( !near(starts, 0) ))
      surv_at_start <- object$helper_funs$distribution_functions$cdf_function(q = starts, gamma = dist_params, lower.tail =FALSE)
    else
      surv_at_start <- rep(1, length(surv))
    return(surv/surv_at_start)
  } else {
    stop(call. = FALSE, "Please report this error to the package maintainer.")
  }

}

#' Fit a bayesian survival model using MAP
#'
#' This function usually isn't called directly. Instead, you'd call \code{survreg_map}, which will
#' call \code{make_survreg_data}, then call this function.
#'
#' @param survreg_data An object from \code{make_survreg_data}
#' @param newdata New data to use in fitting
#' @param optim_args A list of arguments to be passed to \code{optim}.
#'
#' @return An object of class \code{survreg_map}, with print/plot/predict/logLik methods.
#' @export
fit_survreg_map <- function(survreg_data, newdata=NULL,
                            optim_args = list(method = "BFGS",
                                              control = list(trace=as.integer(interactive()),
                                                             maxit=250) )
                            ) {

  if (is.null(newdata))
    newdata <- survreg_data$data

  ## get model-components, standardize
  model_components <- prep_model_components(formula_concat = survreg_data$formula_concat,
                                            forms = survreg_data$forms,
                                            data = newdata,
                                            na.action = survreg_data$na.action,
                                            xlev = survreg_data$xlevels,
                                            standardize_y = survreg_data$log_time_scaling)
  names(model_components$model_mats) <- survreg_data$dist_info$pars_real

  model_mats_std <- map(seq_along(model_components$model_mats),
                        function(i)
                          scale(model_components$model_mats[[i]],
                                center = coalesce(model_components$scale_params[[i]]$center,0),
                                scale = coalesce(model_components$scale_params[[i]]$scale,1))
  )
  names(model_mats_std) <- survreg_data$dist_info$pars_real

  ## unroll inits:
  survreg_data$priors$inits <- purrr::flatten_dbl(survreg_data$inits)
  unrolled_par_init <- with(survreg_data$priors, structure( inits, names = paste0(parameter, "__sep__", term)))

  ## helper functions:
  helper_funs <- list()

  helper_funs$inv_trans <- purrr::map(survreg_data$dist_info$inverse_transforms, get, mode='function')

  helper_funs$get_distribution_params_from_unrolled_coefs <- function(unrolled_par, list_of_standardized_model_mats = NULL) {

    if (is.null(list_of_standardized_model_mats))
      list_of_standardized_model_mats <- model_mats_std

    df_par <- tidyr::separate(enframe(unrolled_par), col = name, into = c('parameter', 'term'), sep = "__sep__")
    par <- purrr::map(split(df_par, df_par$parameter), ~with(.x, structure(value, names=term)))
    par <- par[survreg_data$dist_info$pars_real]

    list_of_pars_per_obs <- list()

    for (i in seq_along(survreg_data$dist_info$pars)) {
      param_name <- survreg_data$dist_info$pars[[i]]
      param_name_real <- survreg_data$dist_info$pars_real[[i]]
      list_of_pars_per_obs[[param_name]] <-
        helper_funs$inv_trans[[param_name]](as.matrix(list_of_standardized_model_mats[[param_name_real]]) %*% matrix( par[[param_name_real]] ) )
    }
    as_data_frame(purrr::map(list_of_pars_per_obs, as.numeric))
  }

  helper_funs$distribution_functions <-
    survreg_data$dist_info$distribution_function_factory(all_knots = survreg_data$all_knots,
                                                         scaling_factor = survreg_data$stan_data$scaling_factor)

  helper_funs$lik_fun <- function(dist_params, Y = NULL) {
    if (is.null(Y))
      Y <- model_components$y
    if (class(Y)[1] == "Survint") {
      if (!is.null(survreg_data$all_knots)) {

        event_lgl <- (Y[,4,drop=TRUE] == 1)
        exact_event_lgl <- event_lgl & dplyr::near(Y[,3,drop=TRUE],Y[,2,drop=TRUE])

        if ( any(exact_event_lgl) ) # need dbasis
          stop("Please report this error to the package maintainer.", call. = FALSE)

        lik_uncond <- helper_funs$distribution_functions$cdf_function(q = Y[,3,drop=TRUE],
                                                                      gamma = as.matrix(dist_params), lower.tail = FALSE)
        lik_at_lb <- helper_funs$distribution_functions$cdf_function(q = Y[,2,drop=TRUE],
                                                                     gamma = as.matrix(dist_params), lower.tail = FALSE)
        lik_uncond[event_lgl] <- (lik_at_lb - lik_uncond)[event_lgl]
        lik_at_start <-  helper_funs$distribution_functions$cdf_function(q = Y[,1,drop=TRUE],
                                                                         gamma = as.matrix(dist_params), lower.tail = FALSE)
      } else {
        # add when you add other distributions
        stop("Please report this error to the package maintainer.", call. = FALSE)
      }
      lik <- lik_uncond / lik_at_start
      lik[is.na(lik)] <- 10*.Machine$double.eps
      lik[near(lik,0)] <- 10*.Machine$double.eps
      lik[lik<0] <- 10*.Machine$double.eps
      return(lik)
    } else {
      # add support for `Surv` response object
      stop(call. = FALSE, "Please report this error to the package maintainer.")
    }
  }

  helper_funs$prior_fun <- function(unrolled_par) {
    df_par <- tidyr::separate(enframe(unrolled_par), col = name, into = c('parameter', 'term'), sep = "__sep__")
    with(left_join(survreg_data$priors, df_par, by = c("parameter", "term")),
         dnorm(x = value, mean = mu, sd = sigma))
  }

  neg_loglik_fun <- function(unrolled_par) {
    dist_params <- helper_funs$get_distribution_params_from_unrolled_coefs(unrolled_par)
    lik <- helper_funs$lik_fun(dist_params)
    - (sum(log(lik)) + sum(log(helper_funs$prior_fun(unrolled_par))))
  }

  ## optimize:
  out <- list()
  out$helper_funs <- helper_funs
  out$optim <- do.call(optim,
                       c(
                         list( par = unrolled_par_init,
                               fn = neg_loglik_fun,
                               hessian = TRUE ),
                         optim_args
                       ))

  ## collect results, add prior
  if (!is.null(out$optim$hessian) &&
      all(!is.na(out$optim$hessian)) &&
      all(!is.nan(out$optim$hessian)) &&
      all(is.finite(out$optim$hessian)) &&
      all(eigen(out$optim$hessian)$values > 0)) {
    out$cov <- solve(out$optim$hessian)
    se <- sqrt(diag(out$cov))
  } else {
    warning(immediate. = TRUE, call. = FALSE,
            "Optimisation has probably not converged - Hessian is not positive definite. ")
    out$cov <- NA
    se <- setNames(rep(NA, length(unrolled_par_init)), nm=names(unrolled_par_init))
  }

  out$res_std <- survreg_data$priors[,c('parameter','term'),drop=FALSE]
  out$res_std$to_join <- with(out$res_std, paste(parameter, term, sep="__sep__"))

  out$res_std <- out$res_std %>%
    left_join(.,enframe(out$optim$par,'to_join','estimate'),by='to_join') %>%
    left_join(.,enframe(se,'to_join','std.err'),by='to_join') %>%
    mutate(ci.low = estimate-std.err*1.96,
           ci.hi = estimate+std.err*1.96,
           std.err=NULL) %>%
    left_join(x = ., y = survreg_data$priors, by=c('parameter','term')) %>%
    select(parameter, term, estimate, ci.low, ci.hi,
           scaled.center = center, scaled.scale = scale,
           prior.mu = mu, prior.sigma = sigma, undo_for_coefs)

  ## get coefficients on raw (non-standardized) scale:
  df_est_real <- out$res_std %>%
    split(., .$parameter) %>%
    map_df(function(data) {
      intercept_lgl <- (data$term == "(Intercept)")
      data[c('estimate','ci.low','ci.hi')] <-
        map(data[c('estimate','ci.low','ci.hi')],
            function(ests) {
              scaling_vec <- ifelse(data$undo_for_coefs, coalesce(data$scaled.scale, 1), 1)
              centering_vec <- ifelse(data$undo_for_coefs, coalesce(data$scaled.center, 0), 0)
              estimate_raw <- ests*scaling_vec
              intercept_raw <- ests[intercept_lgl] - sum( (estimate_raw*centering_vec)[!intercept_lgl] )
              out <- numeric(length(ests))
              out[intercept_lgl] <- intercept_raw
              out[!intercept_lgl] <- estimate_raw[!intercept_lgl]
              out
            })
      data
    }) %>%
    select(-scaled.center, -scaled.scale, -prior.mu,-prior.sigma, -undo_for_coefs)
  out$res_std$undo_for_coefs <- NULL

  ## get coefficients in terms of non-transformed parameters:
  out$res <- map2_df(.x = split(df_est_real, df_est_real$parameter)[survreg_data$dist_info$pars_real],
                     .y = survreg_data$dist_info$inverse_transforms,
                     function(data,func_char) {
                       data[c('estimate','ci.low','ci.hi')] <-
                         map(data[c('estimate','ci.low','ci.hi')],
                             ~get(func_char, mode = 'function')(.x))
                       data
                     }) %>%
    rename(parameter_real = parameter) %>%
    left_join(data_frame(parameter_real = survreg_data$dist_info$pars_real,
                         parameter = survreg_data$dist_info$pars),
              by = 'parameter_real') %>%
    select(-parameter_real) %>%
    select(parameter, everything())

  to_write <- setdiff(names(survreg_data), names(out))
  out[to_write] <- survreg_data[to_write]

  class(out) <- 'survreg_map'
  out

}

#' Set prior
#'
#' @param parameter The parameter
#' @param terms The terms. Can be a character vector, or the result of \code{dplyr::vars}, which can be convenient for setting the prior on all but some terms.
#' @param mu A value for mu, to be appled to all `terms` for the `parameter`
#' @param sigma A value for sigma, to be appled to all `terms` for the `parameter`
#'
#' @return An object of class 'prior', to be (optionally inserted into a list with other priors and)
#'   passed to `priors` arg in `survreg_map`
#' @export
set_prior <- function(parameter = NULL, terms, mu = NULL, sigma = NULL) {
  stopifnot(is.character(terms) || inherits(terms, "col_list"))

  if (is.character(terms)) terms <- paste0("`",terms,"`")
  out <- function(prior_df) {
    if (!is.null(parameter)) {
      if (!parameter %in% prior_df$parameter) {
        warning(call. = FALSE,immediate. = TRUE,
                "Parameter ", parameter,
                " not found; available parameters are: ", paste0(unique(prior_df$parameter), collapse=", "))
        return(prior_df)
      }
      param_lgl <- (prior_df$parameter == parameter)
    } else {
      param_lgl <- rep(TRUE, nrow(prior_df))
    }
    terms_to_modify <- unname(dplyr::select_vars_(vars = unique(prior_df$term[param_lgl]), args = terms))
    if (!is.null(mu))
      prior_df$mu[prior_df$term%in%terms_to_modify & param_lgl] <- mu
    if (!is.null(sigma))
      prior_df$sigma[prior_df$term%in%terms_to_modify & param_lgl] <- sigma
    prior_df
  }
  class(out) <- c('prior',class(out))
  out
}

#' Build a bayesian survival model using MAP
#'
#' @param formula Formula with both rhs and lhs. Rhs will be applied to 'location' parameter of the
#'   distribution, or for spline-based models, the first parameter (gamma0).
#' @param anc A list of lhs-only formulae. Must be named for the parameters of the distribution
#'   (e.g., \code{list(gamma1 = ~ pred1 + pred2)}).
#' @param data A data.frame
#' @param distribution Character naming the distribution
#' @param dist_config A list customizing the distribution. For example, if \code{distribution =
#'   'roy_parm_splines'}, then you can pass \code{list(knots=)} to specify the location of the
#'   knots.
#' @param standardize_y If TRUE (the default), then response variable is standardized. This means
#'   taking the log of the 'end' times, standardizing (subtract mean, divide by sd) them, then
#'   taking the exponent of that. If a list is passed for this argument, it's assumed to have
#'   components 'center' and 'scale', for performing this standardization. Standardization is
#'   helpful because it makes the default priors more meaningful across datasets.
#' @param na.action Function for NAs
#' @param priors An object resulting from \code{set_prior} (or a list of these).
#' @param optim_args A list of arguments to be passed to \code{optim}.
#'
#' @return An object of class \code{survreg_map}, with print/plot/predict/logLik methods.
#' @export
survreg_map <- function(formula,
                        anc = NULL,
                        data,
                        distribution = 'roy_parm_splines',
                        dist_config = list(knots=NULL),
                        standardize_y = TRUE,
                        na.action = na.exclude,
                        priors = NULL,
                        optim_args = list(method = "BFGS", control = list(trace=as.integer(interactive()),
                                                                          maxit=250) )
                        ) {
  survreg_data <- prep_survreg_data(formula = formula, anc = anc, data = data,
                                    distribution = distribution,
                                    dist_config = dist_config,
                                    standardize_y = standardize_y,
                                    na.action = na.action)

  if (!is.null(priors)) {
    if (class(priors)[1]=='prior')
      priors <- list(priors)
    if (is.list(priors)) {
      walk(priors, ~stopifnot(class(.x)[1]=='prior'))
      for (prior in priors)
        survreg_data$priors <- prior(survreg_data$priors)
    } else {
      stop("`priors` should be a list of outputs from `set_prior`.")
    }

  }
  fit_survreg_map(survreg_data, optim_args=optim_args)
}

#' Get Cross-Validated Log-Likelihood
#'
#' @param object An object to be cross-validated
#'
#' @return A vector of log-liklihoods, one for each fold
#' @export
crossv_loglik <- function(object, ...)  UseMethod('crossv_loglik')

#' Get Cross-Validated Log-Likelihood for a survreg_map model
#' @describeIn crossv_loglik
#'
#' @param folds Either (a) the number of folds, or (b) a list of indices for the *test* group.
#' @param seed Allows you to set the seed for reproducible folds. This is essential if you want to
#'   compare cross-validation estimates for different calls to this function.
#' @param mc.cores Passed to \code{parallel::mclapply} for running folds in parallel.
#' @param ... Ignored
#'
#' @export
crossv_loglik.survreg_map <- function(object, folds = 5, seed = NULL, mc.cores = NULL, ...) {
  if (is.numeric(folds)) {
    if (is.null(seed)) stop("Please set the seed.")
    else set.seed(seed)
    num_folds <- folds
    folds <- list()
    n <- nrow(object$data)
    test_n <- floor(n/num_folds)
    remaining <- seq_len(n)
    for (i in seq_len(num_folds-1)) {
      folds[[i]] <- sample(remaining, size = test_n, replace = FALSE)
      remaining <- setdiff(remaining, folds[[i]])
    }
    folds[[num_folds]] <- remaining
  }
  if (is.list(folds)) {
    if (is.null(mc.cores))
      mc.cores <- min(length(folds), parallel::detectCores()-1)
    fits <- parallel::mclapply(X = folds, mc.cores = mc.cores,
             FUN = function(.x)
               update.survreg_map(object, data = object$data[-.x,]))
    conv_fail_lgl <- !map_lgl(map(fits,'cov'), is.matrix)
    out <- map2_dbl(fits, folds, ~logLik(.x, newdata = object$data[.y,]))
    out[conv_fail_lgl] <- NA
    out
  } else {
    stop(call. = FALSE, "`folds` should either be an integer or a list of row-indices for 'test'.")
  }
}

#' Get logliklihood from an object of class \code{survreg_map}.
#'
#' @param object An object of class \code{survreg_map}.
#' @param newdata Optional, newdata on which to calculate the likelihood.
#' @param ... Ignored
#'
#' @return The log-likelihood.
#' @export
logLik.survreg_map <- function(object, newdata = NULL, ...) {
  if (is.null(newdata))
    newdata <- object$data

  model_components <- prep_model_components(formula_concat = object$formula_concat,
                                            forms = object$forms,
                                            data = newdata,
                                            na.action = object$na.action,
                                            xlev = object$xlevels,
                                            standardize_y = object$log_time_scaling)
  names(model_components$model_mats) <- object$dist_info$pars_real

  model_mats_std <- map(seq_along(model_components$model_mats),
                        function(i)
                          scale(model_components$model_mats[[i]],
                                center = coalesce(object$scale_params[[i]]$center,0),
                                scale = coalesce(object$scale_params[[i]]$scale,1)))
  names(model_mats_std) <- object$dist_info$pars_real

  unrolled_par <- with(object$res_std, structure( estimate, names = paste0(parameter, "__sep__", term)))
  dist_params <- object$helper_funs$get_distribution_params_from_unrolled_coefs(unrolled_par,
                                                                                list_of_standardized_model_mats = model_mats_std)
  lik <- object$helper_funs$lik_fun(dist_params, Y = model_components$y)
  return( sum(log(lik)) )
}

#' @export
coef.survreg_map <- function(object, ...) {
  purrr::map(split(object$res, object$res$parameter), ~setNames(.x$estimate, nm = .x$term))
}

#' @export
print.survreg_map <- function(x, standarized = TRUE, ...) {
  if (standarized) {
    cat("\nResults (standardized): ====\n")
    print(x$res_std)
  } else {
    cat("\nResults (unstandardized): ====\n")
    print(x$res)
  }
}

#' Plot coefficients for 'survreg_map'
#' @describeIn plot_coefs
#'
#' @param standardized Should the coefficient-estimates be standardized? If so, the prior will be
#'   plotted in the background.
#'
#' @export
plot_coefs.survreg_map <- function(object, standardized=FALSE, ...) {

  if (standardized) {
    object$res_std$parameter <- factor(object$res_std$parameter, levels = object$dist_info$pars_real)
    ggplot(object$res_std, aes(x=term, y = estimate, color = parameter, shape = ci.low>0|ci.hi<0)) +
      scale_shape_discrete(guide=FALSE)+
      geom_point()+
      geom_linerange(aes(ymin = ci.low, ymax = ci.hi)) +
      facet_wrap(~parameter, scales = 'free') +
      geom_hline(yintercept = 0) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_pointrange(alpha=.25,
                      mapping = aes(y=prior.mu, ymin = prior.mu-prior.sigma*1.96, ymax = prior.mu+prior.sigma*1.96)) +
      guides(colour=FALSE)
  } else {
    object$res$parameter <- factor(object$res$parameter, levels = object$dist_info$pars)
    ggplot(object$res, aes(x=term, y = estimate, color = parameter)) +
      geom_point()+
      geom_linerange(aes(ymin = ci.low, ymax = ci.hi)) +
      facet_wrap(~parameter, scales = 'free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(colour=FALSE)
  }

}

#' Update method for \code{survreg_map}
#'
#' This is a rudimentary update method for \code{survreg_map}. Unlike most update methods in R, this
#' doesn't work by modifying and replaying the call, but instead uses the slots in the object. It's
#' limited to updating the formula(s), the data, and/or the prior(s).
#'
#' @param object An object of class \code{survreg_map}
#' @param ... Args to be updated. Currently just `data` and `priors`. Even if `data` is updated, the
#'   centering/scaling of predictors will *not* be updated, so that the updated model can be
#'   compared to the original apples-to-apples.
#'
#' @return An updated model.
#' @export
update.survreg_map <- function(object, ...) {
  new_args <- list(...)
  if (!all(names(new_args)%in%c('data','priors')))
    stop(call. = FALSE,
         "Currently, the `update` method for `survreg_map` only supports updating the formula(e), data, and/or priors.")
  if ('data' %in% names(new_args))
    object$data <- new_args$data

  if ('priors' %in% names(new_args)) {
    priors <- new_args$priors
      if (class(priors)[1]=='prior')
        priors <- list(priors)
      if (is.list(priors)) {
        walk(priors, ~stopifnot(class(.x)[1]=='prior'))
        for (prior in priors)
          object$priors <- prior(object$priors)
      } else {
        stop("`priors` should be a list of outputs from `set_prior`.")
      }
  }
  fit_survreg_map(object, newdata = new_args$data)

}

