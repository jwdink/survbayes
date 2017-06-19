
# Main Fxn & Internal Helpers ---------------------------------------------------------------------------------

standardize_knots <- function(formula, anc, data, dist_config) {
  get_knot_boundaries <- function(formula, data) {
    response_object <- model.frame(formula = update(formula, .~1), data=data, na.action=na.pass)[[1]]
    if (class(response_object)[1]=="Survint")
      range(response_object[,3,drop=TRUE])
    else if (class(response_object)[1]=="Surv")
      stop(call. = FALSE,
           "Please use `Survint` for response in formula, not `Surv`; this will be fixed in a future release.")
    else
      stop(call. = FALSE,
           "Response in formula not a recognizable survival-object.")
  }
  if ( 'all_knots' %in% names(dist_config)) {
    stopifnot( length(dist_config$all_knots) == (length(anc)+1) )
  } else {
    if (! 'all_knots' %in% names(dist_config)) {
      if (! 'knots' %in% names(dist_config))
        stop(call. = FALSE,
             "Please specify `dist_config$knots` on the original timescale.")
      stopifnot( length(dist_config$knots) == (length(anc)-1) )
      bound_knots <- get_knot_boundaries(formula = formula, data = data)
      dist_config$all_knots <- c(bound_knots[1], dist_config$knots, bound_knots[2])
      dist_config$knots <- NULL
    }
  }
  dist_config
}

get_prior_df <- function(priors, dist_info, model_components) {
  if (is.data.frame(priors)) {
    df_prior <- priors
  } else {
    df_prior <- dist_info$default_priors(model_components$list_of_model_mats_std, model_components$response_object)
    if (!is.null(priors)) {
      if (class(priors)[1]=='prior')
        priors <- list(priors)
      if (is.list(priors)) {
        walk(priors, ~stopifnot(class(.x)[1]=='prior'))
        for (prior in priors)
          df_prior <- prior(df_prior)
      } else {
        stop("`priors` should be a list of outputs from `set_prior`.")
      }
    }
  }
  df_prior
}

survreg_map.fit <- function(model_components, df_prior, dist_info, optim_args) {

  ## get inits, unroll:
  unrolled_par_all <- with(df_prior, structure(location, names=paste0(parameter, "__sep__", term)))
  fixed_lgl <- setNames(df_prior$fixed, nm=names(unrolled_par_all))
  unrolled_par_init <- unrolled_par_all[!fixed_lgl]

  ##
  roll_pars <- function(unrolled_par) {
    df_par <- tidyr::separate(tibble::enframe(unrolled_par), col = name, into = c('parameter', 'term'), sep = "__sep__")
    par <- purrr::map(split(df_par, df_par$parameter), ~with(.x, structure(value, names=term)))
    par <- par[dist_info$pars_real]
    par
  }

  ##
  itrans <- map(dist_info$inverse_transforms, get, mode='function')
  get_distribution_params <- function(par, list_of_mm) {
    par <- par[dist_info$pars_real]
    list_of_mm <- list_of_mm[dist_info$pars]

    list_of_pars_per_obs <- list()
    for (i in seq_along(dist_info$pars)) {
      param_name <- dist_info$pars[[i]]
      param_name_real <- dist_info$pars_real[[i]]
      list_of_pars_per_obs[[param_name]] <-
        itrans[[param_name]](as.matrix(list_of_mm[[param_name]]) %*% matrix( par[[param_name_real]] ) )
    }
    as_data_frame(purrr::map(list_of_pars_per_obs, as.numeric))
  }

  ## distribution functions:
  if (dist_info$spline_dist)
    dfuns <- dist_info$distribution_function_factory(all_log_knots = log(dist_info$config$all_knots))
  else
    dfuns <- dist_info$distribution_function_factory()

  ## likelihood function:
  lik_fun <- function(dist_params, Y) {
    if (class(Y)[1] == "Survint") {
      if (dist_info$spline_dist) {
        event_lgl <- (Y[,4,drop=TRUE] == 1)
        exact_event_lgl <- event_lgl & dplyr::near(Y[,3,drop=TRUE],Y[,2,drop=TRUE])
        lik_uncond <- dfuns$cdf_function(q = Y[,3,drop=TRUE],gamma = as.matrix(dist_params), lower.tail = FALSE)
        lik_at_lb <- dfuns$cdf_function(q = Y[,2,drop=TRUE],gamma = as.matrix(dist_params), lower.tail = FALSE)
        lik_uncond[event_lgl] <- (lik_at_lb - lik_uncond)[event_lgl]
        if ( any(exact_event_lgl) )
          lik_uncond[exact_event_lgl] <- dfuns$pdf_function(q=Y[exact_event_lgl,3,drop=TRUE],
                                                            gamma=as.matrix(dist_params[exact_event_lgl,,drop=FALSE]))
        lik_at_start <- dfuns$cdf_function(q = Y[,1,drop=TRUE],gamma = as.matrix(dist_params), lower.tail = FALSE)
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

  ## prior function(s):
  prior_funs <- split(df_prior, df_prior$parameter) %>%
    purrr::map(function(df) {
      if (n_distinct(df$family)==1) {
        pfun <- get(stan_to_r_pdf(df$family[1]), mode = 'function')
        function(x) pfun(x, df$location, df$spread)
      } else {
        pfuns <- purrr::map(stan_to_r_pdf(df$family), get, mode='function')
        function(x) purrr::map_dbl(seq_along(x), ~pfuns[[i]](x[i], df$location[i], df$spread[i]))
      }
    })

  ##
  neg_loglik_fun <- function(unrolled_par_variable, model_mats, Y) {
    unrolled_par_all[!fixed_lgl] <- unrolled_par_variable
    # dist params, liklihood:
    par <- roll_pars(unrolled_par_all)
    dist_params <- get_distribution_params(par, list_of_mm=model_mats)
    lik <- lik_fun(dist_params, Y = Y)

    # prior:
    prior_p <- purrr::flatten_dbl(purrr::map2(prior_funs, par, ~.x(.y)))

    # neg loglik:
    - (sum(log(lik)) + sum(log(prior_p)))
  }

  ## optimize:
  out <- list()
  optim_args$par <- unrolled_par_init
  optim_args$fn <- neg_loglik_fun
  optim_args$hessian <- TRUE
  optim_args$model_mats <- model_components$list_of_model_mats_std
  optim_args$Y <- model_components$response_object
  out$optim <- do.call(optim,optim_args)

  ## collect results, add prior
  out$estimate <- unrolled_par_all
  out$estimate[!fixed_lgl] <- out$optim$par
  out$cov <- matrix(nrow = length(unrolled_par_all), ncol = length(unrolled_par_all),
                    dimnames = list(names(unrolled_par_all),names(unrolled_par_all)))
  if (!is.null(out$optim$hessian) &&
      all(!is.na(out$optim$hessian)) &&
      all(!is.nan(out$optim$hessian)) &&
      all(is.finite(out$optim$hessian)) &&
      all(eigen(out$optim$hessian)$values > 0)) {
    replace_idx <- matrix(!fixed_lgl,nrow = length(fixed_lgl), ncol = length(fixed_lgl)) &
      matrix(!fixed_lgl,nrow = length(fixed_lgl), ncol = length(fixed_lgl), byrow = TRUE)
    out$cov[replace_idx] <- solve(out$optim$hessian)

  } else {
    warning(immediate. = TRUE, call. = FALSE,
            "Optimisation has probably not converged - Hessian is not positive definite. ")
  }
  out$se <- sqrt(diag(out$cov))

  out$predict  <- function(unrolled_par, model_mats) {
    par <- roll_pars(unrolled_par)
    get_distribution_params(par, list_of_mm=model_mats)
  }

  out$loglik <- function(unrolled_par, model_mats, Y) {
    par <- roll_pars(unrolled_par)
    dist_params <- get_distribution_params(par, list_of_mm=model_mats)
    log(lik_fun(dist_params, Y = Y))
  }

  return(out)

}

prepare_model_components <- function(forms, data, predvars, dist_info,
                                     standardize_x, na.action, contrasts, xlev, drop.unused.levels) {

  ## make model-frame, get mapping: --
  formula_merged <- merge_formulae(forms, data)
  terms_merged <- terms(formula_merged)
  if (!is.null(predvars)) {
    if (is.character(predvars)) predvars <- parse(text = predvars)[[1]]
    if (length(formula_merged)==2)
      predvars <- predvars[-2] # if no response in formula, removed it from predvars
    attr(terms_merged,'predvars') <- predvars
  }
  model_frame_merged <-
    model.frame(terms_merged, data = data, na.action = na.action, drop.unused.levels = drop.unused.levels)
  xlevels <- .getXlevels(attr(model_frame_merged, "terms"), model_frame_merged)
  if (length(xlevels)==0) xlevels<-NULL

  # separate out response object: --
  response_idx <- attr(terms_merged,'response')
  if (response_idx!=0) {
    terms_full <- terms(model_frame_merged)
    response_object <- model_frame_merged[[response_idx]]
    model_frame_merged <- model_frame_merged[,-response_idx,drop=FALSE]
    attr(model_frame_merged,'terms') <- delete.response(terms_full)
  } else {
    response_object <- NULL
  }
  if (is.null(predvars)) predvars <- attr(terms_full,'predvars')

  ## standardize model-frame: --
  mf_is_numeric_lgl <- purrr::map_lgl(model_frame_merged, is.numeric) & !purrr::map_lgl(model_frame_merged, is.matrix)
  standardize_x_arg <- standardize_x
  standardize_x <- list(
    center= purrr::map_dbl(model_frame_merged[,mf_is_numeric_lgl,drop=FALSE], mean, na.rm=TRUE),
    scale= purrr::map_dbl(model_frame_merged[,mf_is_numeric_lgl,drop=FALSE], sd, na.rm=TRUE)
  )
  if (!is.logical(standardize_x_arg)) {
    stopifnot(is.list(standardize_x_arg))
    stopifnot(c('center','scale')%in%names(standardize_x_arg))
    for (nm in c('center','scale')) {
      # check for factor names:
      if ( any( names(standardize_x_arg[[nm]]) %in% names(which(!mf_is_numeric_lgl)) ) )
        stop(call. = FALSE, "The following variables cannot be standardized because they aren't numeric:",
             paste0(deparse(names(which(!mf_is_numeric_lgl & does_something_lgl))),collapse="\n"),
             "\n(Hint: for factors, use contrasts for centering instead.)")

      # add to defaults:
      standardize_x[[nm]][names(standardize_x_arg[[nm]])] <- standardize_x_arg[[nm]]
    }
  } else {
    if (!standardize_x_arg) {
      standardize_x <- FALSE
    }
  }
  model_frame_std <- model_frame_merged
  model_frame_std[names(standardize_x$center)] <- purrr::map(.x = names(standardize_x$center),
                                                             .f= ~model_frame_std[[.x]]-standardize_x$center[[.x]])
  model_frame_std[names(standardize_x$scale)] <- purrr::map(.x = names(standardize_x$scale),
                                                            .f= ~model_frame_std[[.x]]/standardize_x$scale[[.x]])

  ## smarter contrasts: --
  if (is.null(contrasts)) {
    contr_full_lgl <-
      purrr::map_lgl(purrr::map(model_frame_std, attr,'contrasts'), identical, 'contr.full')
    for (col in names(which(contr_full_lgl)))
      attr(model_frame_std[[col]],'contrasts') <- contr.full(n = levels(model_frame_std[[col]]), f = model_frame_std[[col]])
    matrix_contrasts_lgl <- purrr::map_lgl(purrr::map(model_frame_std, attr,'contrasts'),is.matrix)
    if (any(matrix_contrasts_lgl)&!is.null(xlevels))  # any variables with matrix-contrasts doesn't need an xlev (i think it'd override it)
      xlevels <- xlevels[intersect(names(xlevels),names(which(!matrix_contrasts_lgl)))]

    contrasts <- purrr::compact(purrr::map(model_frame_std, attr,'contrasts'))
  }

  ## get model-mat(s): --
  model_matrix_merged_std <-
    model.matrix(terms(model_frame_std), data = model_frame_std, contrasts.arg = contrasts)
  term_mapping <- get_terms_mapping(formula_merged, data, contrasts.arg = contrasts)
  list_of_model_mats_std <- purrr::map(forms, function(this_form) {
    this_terms <- terms(this_form)
    needed_cols <- unique(flatten_chr(term_mapping$term_labels$model_matrix[attr(this_terms,'term.labels')]))
    if (attr(this_terms,'intercept')==1) needed_cols <- c("(Intercept)", needed_cols)
    model_matrix_merged_std[,needed_cols,drop=FALSE]
  })

  # done--
  out <- list(model_frame_merged= model_frame_merged,
              list_of_model_mats_std= list_of_model_mats_std,
              response_object = response_object,
              standardize_x= standardize_x,
              predvars = predvars,
              xlevels = xlevels,
              contrasts = attr(model_matrix_merged_std,'contrasts'))

  return(out)
}


#' Survival Regression with MAP
#'
#' Fits a survival-regression model with MAP (Maximum a posteriori estimation).
#'
#' This function centers and scales all numeric predictors before the fitting process, and all
#' priors are placed on this standardized scaled. There are some advantages to this. First,
#' default-priors are meaningful and applicable for all numeric predictors: they can be interpreted
#' as mean/expectation at no effect of each predictor, with a standard-deviation on this prior equal
#' to the standard-deviation for that predictor. (For factors, you can set the contrasts to
#' \code{contr.full} to acheive a similar effect.) Second, it's easy to set a single value for the
#' 'spread' of the prior on all (non-intercept) predictors, which means it's easy to use priors for
#' their regularization properties -- e.g., trying different values for the 'spread' and picking the
#' one that maximizes cross-validation performance. See \code{crossv_loglik}.
#'
#' The disadvantage to standardizing the predictors is that a little more care is needed in
#' preserving the standardization-parameters across model-calls. Default behavior for R's
#' \code{update} function would recompute these parameters on each call: in this case it would mean
#' refitting with new data (updating nothing else), which would have the side-effect of updating the prior
#' (since updating the data would change the mean and standard-deviation of your predictors). This
#' function avoids this unexpected behavior with the \code{standardize_x}, \code{predvars}, and
#' \code{contrasts} arguments.
#'
#' The \code{standardize_x} argument, if set to TRUE, will center and scale all numeric predictors.
#' If a list specifying standardization-parameters is passed, then these will be used and not
#' recomputed. The \code{update} method for \code{survreg_map} is smart enough to take advantage of
#' this functionality: if only the data is being updated, it will make sure to replace the
#' \code{standardize_x=TRUE} from the original call with \code{standardize_x=[the parameters from
#' the first call]}.
#'
#' Some R transformation-functions perform standardization for you: for example,
#' \code{stats::scale} or \code{stats::poly} (use \code{methods(makepredictcall)} to see them all).
#' This function handles these by saving the 'predvars' attribute after the parameters are first
#' computed. Just like \code{standardize_x}, the update method can then avoid re-computing these.
#'
#' Finally, \code{contr.full} centers your contrast-codes based on the data it sees when first
#' calling this function; if you call \code{update}, these contrasts will be preserved with the
#' \code{contrasts} argument.
#'
#' @param formula Formula, with both rhs and lhs, for the main parameter of the distribution.
#' @param anc A list of formulae, rhs only, named for other parameters of the distribution.
#' @param data A data.frame
#' @param distribution Character-string for distribution.
#' @param dist_config A list with options controlling the distribution.
#' @param na.action Function for NAs
#' @param priors Either object resulting from \code{set_prior} (or a list of these), or a data.frame
#'   with columns 'parameter', 'term', 'family', 'location', 'spread'. An optional logical column
#'   'fixed' can be used to specify parameters that should be fixed at their initial value and not
#'   updated during optimization. Note that priors are placed on the predictors *after* these
#'   predictors are centered and scaled, according to the \code{standardize_x} argument.
#' @param standardize_x Either a logical specifying whether to center/scale numeric predictors, or a
#'   list with names 'center' and 'scale', each of which in turn are named numeric vectors
#'   specifying centering/scaling. Because priors are placed on the standardized scale, this is an
#'   important argument; see Details.
#' @param contrasts Contrasts that will determine how factors are formatted. Often the user doesn't
#'   want to use this argument, but instead it's useful for `update`. See 'Details' and
#'   \code{?update.survreg_map}.
#' @param xlevels The levels for each factor (for when contrasts are not explicit). See
#'   \code{?update.survreg_map}.
#' @param optim_args Arguments to pass to \code{stats::optim}
#' @param predvars The 'predvars' attribute of a terms object. See \code{?update.survreg_map}.
#'
#' @return An object of type \code{survreg_map}
#' @export
survreg_map <-
  function(
    formula,
    anc = NULL,
    data,
    distribution = 'roy_parm_splines',
    dist_config = list(knots=NULL),
    na.action = na.exclude,
    priors = NULL,
    standardize_x = TRUE,
    contrasts = NULL,
    xlevels = NULL,
    optim_args = list(method = "BFGS", control = list(trace=as.integer(interactive()),maxit=250)),
    predvars = NULL
  ) {

    the_call <- match.call()

    stopifnot(is.data.frame(data))

    ## formula, distribution:
    if (is.null(lazyeval::f_lhs(formula))) stop(call. = FALSE, "`formula` has no left-hand side.")
    dist_info <- get_dist_info(distribution, k = max(length(anc),1) )
    forms <- standardize_formulae(formula, anc, dist_info)
    anc <- map(forms[-1], ~.x[-2])
    the_call$anc <- parse(text = deparse(anc))[[1]]

    # deal with knots:
    if (dist_info$spline_dist)
      dist_config <- standardize_knots(formula, anc, data, dist_config)
    the_call$dist_config <- dist_config
    dist_info$config <- dist_config

    ## model-componenets:
    model_components <-
      prepare_model_components(forms, data, predvars=predvars, dist_info, standardize_x, na.action,
                               contrasts=contrasts, xlev = xlevels, drop.unused.levels = FALSE)

    ## prior:
    df_prior <- get_prior_df(priors, dist_info, model_components)

    ## fit:
    fit_res <- survreg_map.fit(model_components, df_prior, dist_info, optim_args)

    ## organize results:
    out <- list()
    out$res_std <- df_prior[,c('parameter','term','family','location','spread'),drop=FALSE]
    out$res_std$to_join <- with(out$res_std, paste(parameter, term, sep="__sep__"))

    out$res_std <- out$res_std %>%
      left_join(.,tibble::enframe(fit_res$estimate,'to_join','estimate'),by='to_join') %>%
      left_join(.,tibble::enframe(fit_res$se,'to_join','std.err'),by='to_join') %>%
      mutate(ci.low = estimate-std.err*1.96,
             ci.hi = estimate+std.err*1.96,
             std.err=NULL) %>%
      left_join(x = .,
                y = data_frame(term = names(model_components$standardize_x[[1]]),
                               center = model_components$standardize_x$center,
                               scale = model_components$standardize_x$scale),
                by='term') %>%
      select(parameter, term, estimate, ci.low, ci.hi,
             scaled.center = center, scaled.scale = scale,
             prior.family = family, prior.location = location, prior.spread = spread)

    ## get coefficients on raw (non-standardized) scale:
    df_est_real <- out$res_std %>%
      split(., .$parameter) %>%
      map_df(function(data) {
        intercept_lgl <- (data$term == "(Intercept)")
        data[c('estimate','ci.low','ci.hi')] <-
          map(data[c('estimate','ci.low','ci.hi')],
              function(ests) {
                scaling_vec <- coalesce(data$scaled.scale, 1)
                centering_vec <- coalesce(data$scaled.center, 0)
                estimate_raw <- ests*scaling_vec
                intercept_raw <- ests[intercept_lgl] - sum( (estimate_raw*centering_vec)[!intercept_lgl] )
                out <- numeric(length(ests))
                out[intercept_lgl] <- intercept_raw
                out[!intercept_lgl] <- estimate_raw[!intercept_lgl]
                out
              })
        data
      }) %>%
      select(-scaled.center, -scaled.scale, -matches('prior\\.'))

    ## get coefficients in terms of non-transformed parameters:
    out$res <- map2_df(.x = split(df_est_real, df_est_real$parameter)[dist_info$pars_real],
                       .y = dist_info$inverse_transforms,
                       function(data,func_char) {
                         data[c('estimate','ci.low','ci.hi')] <-
                           map(data[c('estimate','ci.low','ci.hi')],
                               ~get(func_char, mode = 'function')(.x))
                         data
                       }) %>%
      rename(parameter_real = parameter) %>%
      left_join(data_frame(parameter_real = dist_info$pars_real,
                           parameter = dist_info$pars),
                by = 'parameter_real') %>%
      select(-parameter_real) %>%
      select(parameter, everything())

    #
    out[c('loglik','optim','predict','cov')] <- fit_res[c('loglik','optim','predict','cov')]
    out[c('standardize_x','predvars','contrasts','xlevels')] <-
      model_components[c('standardize_x','predvars','contrasts','xlevels')]
    out$data <- data
    out$df_prior <- df_prior
    out$forms <- forms
    out$formula_merged <- merge_formulae(forms, data)
    out$na.action <- na.action
    out$dist_info <- dist_info
    out$call <- the_call

    class(out) <- 'survreg_map'

    return(out)

  }


# Survreg Map Methods -------------------------------------------------------------------------

#' @export
formula.survreg_map <- function(x, all=FALSE, ...) {
  if (all)
    x$forms
  else
    eval(x$call$formula, envir = environment(x$forms[[1]]))
}

#' @export
terms.survreg_map <- function(x, ...) {
  terms(x$formula_merged)
}

#' Update method for \code{survreg_map}
#'
#' A critical difference between this update method and the typical method for other R objects
#' concerns how centering/scaling of the data is (not) updated. The standardization of variables is
#' important in \code{survreg_map}, because priors are placed on the standardized scale. So if
#' updating the data meant re-computing this standardization, the updated vs. original models would
#' effectively have different priors. Instead, for this update method, whether the standardization
#' is recomputed is controlled by the `reeval_scaling_and_terms` argument. If only data are being
#' updated, this defaults to FALSE, meaning that standardization will *not* be recomputed. This
#' includes both scaling accomplished by code within \code{survreg_map} (controlled by the
#' \code{standardize_x} argument), as well as scaling accomplished by transformations in the
#' formula, such as \code{stats::scale}, \code{stats::poly}, \code{splines::ns}, and other functions
#' with \code{makepredictcall} methods (this is controlled by the \code{predvars} argument).
#'
#' @param object Object of class \code{survreg_map}
#' @param formula. Passed to \code{update.formula} for the main formula argument.
#' @param anc. A named list of formulas, each of which are passed to \code{update.formula} for the
#'   corresponding anc formula.
#' @param reeval_scaling_and_terms Should scaling and terms be re-computed with new data? For
#'   details on scaling see 'Description'. This also controls whether the 'contrasts' and 'xlevels'
#'   arguments are updated. Broadly, the idea is to be cross-validation friendly: for example, if
#'   new factor-levels are present in a validation-fold that weren't in the training fold, that's OK
#'   because all levels of the factor were remembered from the original model-object.
#' @param ... Slots to update.
#' @param evaluate If true evaluate the new call else return the call.
#'
#' @return If evaluate = TRUE the fitted object, otherwise the updated call.
#' @export
update.survreg_map <- function(object, formula. = NULL, anc. = NULL,
                               reeval_scaling_and_terms = NULL,
                               ...,
                               evaluate = TRUE) {
  call <- getCall(object)

  the_dots <- pryr::dots(...)
  update_data <- "data" %in% names(the_dots)
  update_others <- (length(the_dots)>as.numeric(update_data)) | (!is.null(formula.)) | (!is.null(anc.))

  if (is.null(reeval_scaling_and_terms)) {
    if (update_others) {
      reeval_scaling_and_terms <- TRUE
      if (update_data)
        warning(call. = FALSE,
                "Other components are being updated aside from the data, so will re-evaluate scaling/terms. ",
                "To override this, set `reeval_scaling_and_terms` to FALSE")
    } else {
      reeval_scaling_and_terms <- FALSE
    }
  } else {
    if (!reeval_scaling_and_terms&update_others&!update_data)
      warning(call. = FALSE,
              "You're trying to use the original scaling/terms, even though you're not updated the data component.")
  }

  if (!reeval_scaling_and_terms & update_others & update_data) {
    object <- do.call(update, c(list(object=quote(object), formula. = formula., anc.=anc., reeval_scaling_and_terms=TRUE),
                                the_dots[names(the_dots)!='data']) )
    call <- getCall(object)
    browser()
    formula. <- NULL
    anc. <- NULL
    the_dots <- the_dots[names(the_dots)=='data']
  }

  if (!is.null(formula.))
    call$formula <- update.formula(formula(object), formula.)

  if (!is.null(anc.)) {
    original_anc <- eval(call$anc, envir = environment(formula(object)))
    if (!all(names(anc.) %in% names(original_anc) ))
      stop(call. = FALSE, "The `anc.` argument should be a list with names:\n",
           paste0(collapse="\n",deparse(names(original_anc))))
    anc. <- purrr::map(anc., ~.x[-2])
    for (p in names(anc.))
      original_anc[[p]] <- update.formula(original_anc[[p]], anc.[[p]])
    call$anc <- parse(text = deparse(original_anc))[[1]]
  }

  if (!reeval_scaling_and_terms) {
    call$standardize_x <- parse(text = deparse(object$standardize_x))[[1]]
    call$contrasts <- parse(text = deparse(object$contrasts))[[1]]
    call$xlevels <- parse(text = deparse(object$xlevels))[[1]]
    call$predvars <- deparse(object$predvars) # not sure how to 'double-escape', so prep-model-components just handles it
  }

  new_call <- pryr::modify_call(call, the_dots)
  if (evaluate)
    eval(new_call, envir = parent.frame(), enclos = environment(formula(object)))
  else
    new_call
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

  ## get estimates, unroll:
  unrolled_par_all <- with(object$res_std, structure(estimate, names = paste0(parameter, "__sep__", term)))

  ## model-componenets:
  forms_rhs <- purrr::map(object$forms, ~.x[c(1,3)])
  model_components <-
    prepare_model_components(forms_rhs, newdata, predvars=object$predvars, object$dist_info, object$standardize_x, object$na.action,
                             contrasts=object$contrasts, xlev = object$xlevels, drop.unused.levels = FALSE)

  df_params<-
    object$predict(unrolled_par = unrolled_par_all, model_mats = model_components$list_of_model_mats_std)

  if (type=='parameters')
    return(as.matrix(df_params))

  if (object$dist_info$spline_dist)
    dfuns <- object$dist_info$distribution_function_factory(all_log_knots = log(object$dist_info$config$all_knots))
  else
    dfuns <- object$dist_info$distribution_function_factory()

  if (length(times)==1)
    times <- rep(x=times, times=nrow(newdata))
  if (is.null(starts))
    starts <- rep(0, length(times))

  na_dropped_idx <- attr(model_components$model_frame_merged,'na.action')
  if (length(na_dropped_idx)>0) {
    times <- times[-na_dropped_idx]
    starts <- starts[-na_dropped_idx]
  }

  if (object$dist_info$spline_dist) {
    if (type == 'survival') {
      surv <- dfuns$cdf_function(q = times, gamma = as.matrix(df_params), lower.tail =FALSE)
      if (any( !near(starts, 0) ))
        surv_at_start <- dfuns$cdf_function(q = starts, gamma = as.matrix(df_params), lower.tail =FALSE)
      else
        surv_at_start <- rep(1, length(surv))
      return(surv/surv_at_start)
    }
  } else {
    stop(call. = FALSE, "Please report this error to the package maintainer.")
  }

}

#' Method for extracting logLiklihood from \code{survreg_map}.
#'
#' @param object Object of class \code{survreg_map}
#' @param newdata New data to compute the log-liklihood on.
#' @param ... Ignored
#'
#' @return Log-liklihood
#' @export
logLik.survreg_map <- function(object, newdata = NULL, ...) {
  if (is.null(newdata))
    newdata <- object$data

  ## get estimates, unroll:
  unrolled_par_all <- with(object$res_std, structure(estimate, names = paste0(parameter, "__sep__", term)))

  ## model-componenets:
  model_components <-
    prepare_model_components(object$forms, newdata, predvars=object$predvars, object$dist_info, object$standardize_x, object$na.action,
                             contrasts=object$contrasts, xlev = object$xlevels, drop.unused.levels = FALSE)

  ll <- object$loglik(unrolled_par = unrolled_par_all,
                      model_mats = model_components$list_of_model_mats_std,
                      Y = model_components$response_object)
  sum(ll)
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
                                 update.survreg_map(object, data = object$data[-.x,], reeval_scaling_and_terms = FALSE))
    conv_fail_lgl <- purrr::map_lgl(.x = purrr::map(fits,'cov'), .f = ~all(is.na(.x)))
    out <- purrr::map2_dbl(fits, folds, ~logLik(.x, newdata = object$data[.y,]))
    out[conv_fail_lgl] <- NA
    out
  } else {
    stop(call. = FALSE, "`folds` should either be an integer or a list of row-indices for 'test'.")
  }
}

#' Plot coefficients for 'survreg_map'
#' @describeIn plot_coefs
#'
#' @param standardized Should the coefficient-estimates be standardized? If so, the prior will be
#'   plotted in the background.
#'
#' @export
plot_coefs.survreg_map <- function(object, standardized=TRUE, ...) {

  if (standardized) {
    object$res_std$parameter <- factor(object$res_std$parameter, levels = object$dist_info$pars_real)
    ggplot(object$res_std, aes(x=term, y = estimate, color = parameter, shape = coalesce(ci.low>0|ci.hi<0,FALSE)) ) +
      scale_shape_discrete(guide=FALSE)+
      geom_point()+
      geom_linerange(aes(ymin = ci.low, ymax = ci.hi)) +
      facet_wrap(~parameter, scales = 'free') +
      geom_hline(yintercept = 0) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_pointrange(alpha=.25,
                      mapping = aes(y=prior.location,
                                    ymin = prior.location-prior.spread*1.96,
                                    ymax = prior.location+prior.spread*1.96)) +
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
