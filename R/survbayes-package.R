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
#'             walk walk2
#' @importFrom graphics plot
#' @importFrom stats sd terms update integrate delete.response
#'             as.formula coef fitted glm lm median na.omit
#'             predict var quantile model.response model.frame
#'             na.pass na.exclude na.fail model.matrix model.weights
#'             .getXlevels .checkMFClasses reformulate logLik
#' @importFrom grDevices dev.off pdf
#' @importFrom glue glue
#'
#' @docType package
#' @name survbayes
NULL
#> NULL


# SurvReg Data --------------------------------------------------------------------------------

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
  model_frame <- model.frame(formula = formula_concat, data = data, na.action = na.action, drop.unused.levels = TRUE, xlev = xlev)
  xlevels <- .getXlevels(attr(model_frame, "terms"), model_frame)
  form_terms <- purrr::map(forms, ~delete.response(terms(.x, data=data)) )
  model_mats <- purrr::map(form_terms, ~model.matrix(.x, model_frame) )

  classes <- map_chr(map(model_frame, class), 1)
  if (any(classes %in% c('character','factor','logical')))
    warning(call. = FALSE,
            "Some factor-variables found. You should set a prior for these manually; ",
            "or use a contrast-scheme where the default prior makes sense.")

  scale_params <- map(seq_along(model_mats), function(i) {
    mm <- model_mats[[i]]
    the_terms <- form_terms[[i]]
    # for each model-mat term, get all the cols in the original data that went into it:
    mapping_from_mm_col_to_data_col <-
      map(attr(mm,'assign'), function(assign_idx) names(which(1==attr(the_terms,'factors')[,assign_idx])))
    # only consider model-mat term to be numeric if all the cols that went into it were numeric:
    is_numeric_term <- map_lgl(mapping_from_mm_col_to_data_col, function(col)
      if (length(col)>0 && col %in% colnames(model_frame)) all(map_lgl(as.list(model_frame[,col,drop=FALSE]),is.numeric)) else FALSE)
    names(is_numeric_term) <- colnames(mm)

    # scale/center the numeric ones, set 0,1 for others (meaning no effect)
    scaled_mm_subset <- scale(mm[,names(which(is_numeric_term)),drop=FALSE])
    centers <- setNames(rep_along(colnames(mm),0), nm = colnames(mm))
    centers[is_numeric_term] <- attr(scaled_mm_subset,"scaled:center")
    scales <- setNames(rep_along(colnames(mm),1), nm = colnames(mm))
    scales[is_numeric_term] <- attr(scaled_mm_subset,"scaled:scale")
    list(center=centers,scale=scales)
  })

  y <- model.extract(model_frame, "response")
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
  df_scale_params <-
    map_df(map(model_components$scale_params,~map_df(.x,.id = 'type',enframe,'term')), ~spread(.x,type,value),.id = 'parameter')
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
                                                center = model_components$scale_params[[i]]$center,
                                                scale = model_components$scale_params[[i]]$scale)
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
  print(x$priors)
}

#' Concatenate a list of formula
#'
#' @param formula Main formula
#' @param forms List of rhs-only formula, from `anc` argument.
#'
#' @return Formula
concatenate_formula <- function(formula, forms, data) {
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


# General Helpers -------------------------------------------------------------------------------------

#' Plot model coefficients
#'
#' @param x A model object
#' @param ... Other arguments to be passed to methods
#'
#' @return A ggplot object
#' @export
plot_coefs <- function(x, ...) {
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

  if (!is.null(k)) {
    dist_infos[['roy_parm_splines']] <- tibble::lst(
      name = 'Royston-Parmesian Cumulative-Splines',
      location = 'gamma0',
      pars = paste0('gamma', 0:k),
      transforms = as.list(structure(names=pars, c('identity', 'log', rep('identity',k-1)))),
      inverse_transforms = as.list(structure(names=pars, c('identity', 'exp', rep('identity',k-1)))),
      pars_real = paste0(ifelse(transforms=="identity","",paste0(transforms,"_")),pars),

      default_priors = function(model_mats, y) {
        param_names <- names(model_mats)

        gamma0_init <- -log(mean(y[,3]))

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

      distribution_function_factory = function(all_knots, y, scaling_factor) {

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

  dist <- match.arg(arg = dist, choices = names(dist_infos))

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

#' Fit a bayesian survival model using MAP
#'
#' @param survreg_data An object from \code{make_survreg_data}
#'
#' @return An object of class \code{survreg_map}, with print/plot/predict/logLik methods.
#' @export
fit_survreg_map <- function(survreg_data, newdata=NULL) {

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
                                center = survreg_data$scale_params[[i]]$center,
                                scale = survreg_data$scale_params[[i]]$scale))
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
  out$optim <- optim(par = unrolled_par_init, fn = neg_loglik_fun, method = "BFGS", hessian = TRUE, control = list(trace=1))

  ## collect results, add prior:
  fisher_info <- solve(out$optim$hessian)
  out$res_std <- survreg_data$priors[,c('parameter','term'),drop=FALSE]
  out$res_std$to_join <- with(out$res_std, paste(parameter, term, sep="__sep__"))

  out$res_std <- out$res_std %>%
    left_join(.,enframe(out$optim$par,'to_join','estimate'),by='to_join') %>%
    left_join(.,enframe(sqrt(diag(fisher_info)),'to_join','std.err'),by='to_join') %>%
    mutate(ci.low = estimate-std.err*1.96,
           ci.hi = estimate+std.err*1.96,
           std.err=NULL) %>%
    left_join(x = ., y = survreg_data$priors, by=c('parameter','term')) %>%
    select(parameter, term, estimate, ci.low, ci.hi,
           scaled.center = center, scaled.scale = scale,
           prior.mu = mu, prior.sigma = sigma)

  ## get coefficients on raw (non-standardized) scale:
  df_est_real <- out$res_std %>%
    split(., .$parameter) %>%
    map_df(function(data) {
      intercept_lgl <- (data$term == "(Intercept)")
      data[c('estimate','ci.low','ci.hi')] <-
        map(data[c('estimate','ci.low','ci.hi')],
            function(ests) {
              estimate_raw <- ests*data$scaled.scale
              intercept_raw <- ests[intercept_lgl] - sum( (estimate_raw*data$scaled.center)[!intercept_lgl] )
              out <- numeric(length(ests))
              out[intercept_lgl] <- intercept_raw
              out[!intercept_lgl] <- estimate_raw[!intercept_lgl]
              out
            })
      data
    }) %>%
    select(-scaled.center, -scaled.scale, -prior.mu,-prior.sigma)

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

  out[names(survreg_data)] <- survreg_data

  class(out) <- 'survreg_map'
  out

}

#' Get logliklihood from an object of class \code{survreg_map}.
#'
#' @param object An object of class \code{survreg_map}.
#' @param newdata Optional, newdata on which to calculate the likelihood.
#' @param ...
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
                                center = object$scale_params[[i]]$center,
                                scale = object$scale_params[[i]]$scale))
  names(model_mats_std) <- object$dist_info$pars_real

  unrolled_par <- with(object$res_std, structure( estimate, names = paste0(parameter, "__sep__", term)))
  dist_params <- object$helper_funs$get_distribution_params_from_unrolled_coefs(unrolled_par,
                                                                                list_of_standardized_model_mats = model_mats_std)
  lik <- object$helper_funs$lik_fun(dist_params, Y = model_components$y)
  return( sum(log(lik)) )
}

coef.survreg_map <- function(object, ...)
  purrr::map(split(object$res, object$res$parameter), ~setNames(.x$estimate, nm = .x$term))

print.survreg_map <- function(x, ...) {
  cat("\nResults (unstandardized): ====\n")
  print(x$res)

  cat("\nResults (standardized): ====\n")
  print(x$res_std)
  cat("\n")
  print.survreg_data(x, title=FALSE)
}

plot_coefs.survreg_map <- function(object, standardized=FALSE) {

  if (standardized) {
    object$res_std$parameter <- factor(object$res_std$parameter, levels = object$dist_info$pars_real)
    ggplot(object$res_std, aes(x=term, y = estimate, color = parameter)) +
      geom_pointrange(aes(ymin = ci.low, ymax = ci.hi)) +
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
      geom_pointrange(aes(ymin = ci.low, ymax = ci.hi)) +
      facet_wrap(~parameter, scales = 'free') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(colour=FALSE)
  }

}

