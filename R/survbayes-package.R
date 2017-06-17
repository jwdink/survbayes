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


#' Survival Regression with MAP
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
#'   updated during optimization.
#' @param standardize_x If this is a logical, it controls whether numeric variables in the
#'   model-frame should be centered and scaled before being passed for optimization, which is
#'   recommended to improve optimization performance, and because it makes default (standard-normal)
#'   priors meaningful/applicable. Note that matrix variables--like those generated from
#'   \code{stats::scale} or \code{states::poly}--will not be (re)centered/scaled. This argument can
#'   also be a list with elements 'center' and 'scale', which are each named vectors for each
#'   numeric element of the model-frame. When the model is fit, the call will be modified to replace
#'   `TRUE` with such a list. This means that in calls to \code{update.survreg_map}, the same
#'   standardization will apply to new data (again, with the exception of matrix-variables; to save
#'   the parameters for these you can pass a `terms` object to `update`; see
#'   \code{?update.survreg_map}).
#' @param contrasts Contrasts that will determine how factors are formatted. Often the user doesn't
#'   want to use this argument, but instead it's useful for `update`. See
#'   \code{?update.survreg_map}.
#' @param xlevels The levels for each factor (for when contrasts are not explicit). See
#'   \code{?update.survreg_map}.
#' @param optim_args Arguments to pass to \code{stats::optim}
#' @param terms A terms object. See \code{?update.survreg_map}.
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
    terms = NULL
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
      prepare_model_components(forms, data, terms=terms, dist_info, standardize_x, na.action,
                               contrasts=contrasts, xlev = xlevels, drop.unused.levels = FALSE)
    the_call$standardize_x <- parse(text = deparse(model_components$standardize_x))[[1]]
    the_call$contrasts <- parse(text = deparse(model_components$contrasts))[[1]]
    the_call$xlevels <- parse(text = deparse(model_components$xlevels))[[1]]

    ## prior:
    df_prior <- get_prior_df(priors, dist_info, model_components)

    ## fit:
    fit_res <- survreg_map.fit(model_components, df_prior, dist_info, optim_args)

    ## organize results:
    out <- list()
    out$res_std <- df_prior[,c('parameter','term','family','location','spread'),drop=FALSE]
    out$res_std$to_join <- with(out$res_std, paste(parameter, term, sep="__sep__"))

    out$res_std <- out$res_std %>%
      left_join(.,enframe(fit_res$optim$par,'to_join','estimate'),by='to_join') %>%
      left_join(.,enframe(fit_res$se,'to_join','std.err'),by='to_join') %>%
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

    out[c('loglik','optim')] <- fit_res[c('loglik','optim')]
    out$call <- the_call
    out$terms <- model_components$terms
    out$dist_info <- dist_info

    class(out) <- 'survreg_map'

    return(out)

  }

formula.survreg_map <- function(x, ...) {
  eval(x$call$formula)
}

update.survreg_map <- function(object, formula. = NULL, anc. = NULL, ...) {
  call <- getCall(object)

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

  new_call <- pryr::modify_call(call, pryr::dots(...))
  eval(new_call, environment(formula(object)))
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
    ggplot(object$res_std, aes(x=term, y = estimate, color = parameter, shape = ci.low>0|ci.hi<0)) +
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
  inits <- purrr::map(split(df_prior, df_prior$parameter), ~with(.x, structure(location,names=term)))
  inits <- purrr::map(inits, ~array(.x, dim = length(.x))) # really just needed for stan
  unrolled_par_init <- purrr::flatten_dbl(inits)
  unrolled_par_init <- with(df_prior, structure( unrolled_par_init, names = paste0(parameter, "__sep__", term)))

  ##
  roll_pars <- function(unrolled_par) {
    df_par <- tidyr::separate(enframe(unrolled_par), col = name, into = c('parameter', 'term'), sep = "__sep__")
    par <- purrr::map(split(df_par, df_par$parameter), ~with(.x, structure(value, names=term)))
    par <- par[dist_info$pars_real]
  }

  ##
  itrans <- map(dist_info$inverse_transforms, get, mode='function')
  get_distribution_params <- function(par, list_of_mm) {
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
    dfuns <- dist_info$distribution_function_factory(all_log_knots = dist_info$config$all_knots)
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
  neg_loglik_fun <- function(unrolled_par, model_mats, Y) {
    # dist params, liklihood:
    par <- roll_pars(unrolled_par)
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
  if (!is.null(out$optim$hessian) &&
      all(!is.na(out$optim$hessian)) &&
      all(!is.nan(out$optim$hessian)) &&
      all(is.finite(out$optim$hessian)) &&
      all(eigen(out$optim$hessian)$values > 0)) {
    out$cov <- solve(out$optim$hessian)
    out$se <- sqrt(diag(out$cov))
  } else {
    warning(immediate. = TRUE, call. = FALSE,
            "Optimisation has probably not converged - Hessian is not positive definite. ")
    out$cov <- NA
    out$se <- setNames(rep(NA, length(unrolled_par_init)), nm=names(unrolled_par_init))
  }

  out$loglik <- function(unrolled_par, model_mats, Y) {
    par <- roll_pars(unrolled_par)
    dist_params <- get_distribution_params(par, list_of_mm=model_mats)
    lik_fun(dist_params, Y = Y)
  }

  return(out)

}

stan_to_r_pdf <- function(stan_pdf) {
  mapper <- c(`normal` = 'dnorm', `double_exponential` = 'dlaplace')
  mapper[stan_pdf]
}

prepare_model_components <- function(forms, data, terms, dist_info,
                                     standardize_x, na.action, contrasts, xlev, drop.unused.levels) {

  ## make model-frame, get mapping: --
  formula_merged <- merge_formulae(forms, data)
  model_frame_merged <-
    model.frame(formula_merged, data = data, na.action = na.action, drop.unused.levels = drop.unused.levels)
  xlevels <- .getXlevels(attr(model_frame_merged, "terms"), model_frame_merged)

  if (is.null(terms))
    terms <- terms(model_frame_merged)

  # separate out response object: --
  response_idx <- attr(terms(formula_merged),'response')
  if (response_idx!=0) {
    response_object <- model_frame_merged[[response_idx]]
    model_frame_merged <- model_frame_merged[,-response_idx,drop=FALSE]
    attr(model_frame_merged,'terms') <- delete.response(terms)
  } else {
    response_object <- NULL
  }

  ## standardize model-frame: --
  mf_is_numeric_lgl <- purrr::map_lgl(model_frame_merged, is.numeric)
  standardize_x_arg <- standardize_x
  standardize_x <- list(center = map_dbl(as.data.frame(model_frame_merged[mf_is_numeric_lgl]), ~0),
                        scale = map_dbl(as.data.frame(model_frame_merged[mf_is_numeric_lgl]), ~1) )

  if (!is.logical(standardize_x_arg)) {
    stopifnot(is.list(standardize_x_arg))
    stopifnot(c('center','scale')%in%names(standardize_x_arg))
    for (nm in c('center','scale')) {
      # check for names not there:
      unexpected_names <- setdiff(names(standardize_x_arg[[nm]]), colnames(model_frame_merged))
      if (length(unexpected_names)>0)
        stop(call. = FALSE,
             "The following names are not in the model.frame:\n",
             paste0(deparse(unexpected_names),collapse="\n"),
             "\n\nNames in model.matrix:\n",
             paste0(deparse(colnames(model_frame_merged)),collapse="\n"),
             "\n(if you'd like to center a factor variable, use contrasts)")

      # check for factor names:
      if ( any( names(standardize_x_arg[[nm]]) %in% names(which(!mf_is_numeric_lgl)) ) )
        stop(call. = FALSE, "The following variables cannot be standardized because they aren't numeric:",
                paste0(deparse(names(which(!mf_is_numeric_lgl & does_something_lgl))),collapse="\n"),
                "\n(Hint: for factors, use contrasts for centering instead.)")

      # add to defaults:
      standardize_x[[nm]][names(standardize_x_arg[[nm]])] <- standardize_x_arg[[nm]]
    }
  } else {
    if (standardize_x_arg) {
      standardize_x$center <- map_dbl(model_frame_merged[,mf_is_numeric_lgl,drop=FALSE], mean, na.rm=TRUE)
      standardize_x$scale <- map_dbl(model_frame_merged[,mf_is_numeric_lgl,drop=FALSE], sd, na.rm=TRUE)
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
    if (any(matrix_contrasts_lgl))  # any variables with matrix-contrasts doesn't need an xlev (i think it'd override it)
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
              terms = terms,
              xlevels = xlevels,
              contrasts = attr(model_matrix_merged_std,'contrasts'))

  return(out)
}

#' Get mapping from column names to model-terms and vice versa
#'
#' @param formula The formula to be passed to model.frame
#' @param data The data.frame
#' @param ... Arguments to be passed to \code{stats::model.matrix}
#'
#' @return A nested list, where the first level is 'from', the second level is 'to'.
#' @export
get_terms_mapping <- function(formula, data, ...) {

  flip_mapping <- function(mapping_list)  {
    data_frame(map_from = names(mapping_list), map_to = mapping_list) %>%
      tidyr::unnest() %>%
      group_by(map_to) %>%
      dplyr::do(model_mat_cols = .$map_from) %>%
      ungroup() %>%
      tibble::deframe()
  }

  lazyeval::f_lhs(formula) <- NULL
  model_frame <- model.frame(formula = formula, data, na.action=na.pass)
  model_matrix <- do.call(model.matrix, c(list(terms(model_frame), data=model_frame), list(...)))

  if (ncol(model_frame)==0) {
    # check for intercept?
    from_mm_to_od <- from_od_to_mm <- from_mf_to_od <- list()
    stop(call. = FALSE, "Please report this error to the package maintainer.")
  } else {
    # for each column in the model.matrix, get corresponding model.frame col(s):
    cols_in_mm <- colnames(model_matrix)
    fact_mat <- attr(terms(model_frame),'factors')
    from_mm_to_mf <- map(attr(model_matrix,'assign'),
                         function(assign_idx) row.names(fact_mat)[1==fact_mat[,assign_idx,drop=TRUE]])
    names(from_mm_to_mf) <- cols_in_mm

    # for each term.label, get corresponding model.frame col(s):
    from_tl_to_mf <-
      purrr::map(.x = seq_len(ncol(fact_mat)), .f = ~row.names(fact_mat)[1==fact_mat[,.x,drop=TRUE]])
    names(from_tl_to_mf) <- colnames(fact_mat)

    # for each column in the model.matrix, get corresponding original-data col(s):
    is_null_lgl <- map_lgl(from_mm_to_mf, is.null)
    from_mm_to_od <- vector(mode = 'list', length = length(from_mm_to_mf))
    names(from_mm_to_od) <- names(from_mm_to_mf)
    if (any(!is_null_lgl)) {
      from_mm_to_od[!is_null_lgl] <- map(
        .x = from_mm_to_mf[!is_null_lgl],
        .f = function(vec_of_mm_cols) flatten_chr(map(vec_of_mm_cols, ~all.vars(parse(text = .x)[[1]]))))
      # we parse the syntax for each term, getting any corresponding to a variable. then we have to check
      # if that variable is a column in the original-data (e.g., in case a term is
      # `poly(column, degree=my_variable_specifying_degree_that_isnt_a_column)`
      from_mm_to_od[!is_null_lgl] <- map(from_mm_to_od[!is_null_lgl], ~.x[.x%in%colnames(data)])
    }

    # for each column in the model.frame, get corresponding original-data col(s):
    from_mf_to_od <- map(colnames(model_frame), ~all.vars(parse(text = .x)[[1]])) %>%
      map(~.x[.x%in%colnames(data)])
    names(from_mf_to_od) <- colnames(model_frame)
  }

  #
  list(
    original_data = list(
      term_labels = purrr::map(flip_mapping(from_mf_to_od), ~purrr::flatten_chr(flip_mapping(from_tl_to_mf)[.x])),
      model_frame = flip_mapping(from_mf_to_od),
      model_matrix = flip_mapping(from_mm_to_od)
    ),
    term_labels = list(
      original_data = purrr::map(from_tl_to_mf, ~purrr::flatten_chr(from_mf_to_od[.x])),
      model_frame = from_tl_to_mf,
      model_matrix = purrr::map(from_tl_to_mf, ~purrr::flatten_chr(flip_mapping(from_mm_to_mf)[.x]))
    ),
    model_frame = list(
      original_data = from_mf_to_od,
      term_labels = flip_mapping(from_tl_to_mf),
      model_matrix = flip_mapping(from_mm_to_mf)
    ),
    model_mat = list(
      original_data = from_mm_to_od,
      term_labels = purrr::map(from_mm_to_mf, ~purrr::flatten_chr(flip_mapping(from_tl_to_mf)[.x])),
      model_frame = from_mm_to_mf
    )
  )

}

standardize_formulae <- function(formula, anc, dist_info) {
  if (is.null(anc)) {
    anc <- structure(names = dist_info$pars[-1], purrr::rerun(length(dist_info$pars[-1]), ~1))
  } else {
    anc_arg <- anc
    if (!is.list(anc_arg) || is.null(names(anc_arg)) || any(names(anc_arg)=="") )
      stop("`anc` must be a named list")
    if (!all(purrr::map_lgl(anc_arg, ~inherits(.x, "formula"))))
      stop("`anc` must be a list of formulae")
    anc <- structure(names = dist_info$pars[-1], purrr::rerun(length(dist_info$pars[-1]), ~1))
    anc[names(anc_arg)] <- anc_arg
    unmatched_params <- setdiff( names(anc_arg), dist_info$pars)
    if (length(unmatched_params)>0)
      stop("The following names in your `anc` argument are not parameters of '", distribution, "':", paste0(", ", unmatched_params))
  }
  out <- c(list(formula), anc)
  names(out)[1] <- dist_info$pars[1]
  # this is necessary so `response ~ .` works as expected:
  for (i in seq_along(out)[-1]) lazyeval::f_lhs(out[[i]]) <- lazyeval::f_lhs(out[[1]])
  return(out)
}

# concatenate_formulae <- function(forms, data) {
#   list_of_term_labels <- purrr::map(purrr::map(forms, terms, data=data), attr, 'term.labels')
#   purrr::map2(list_of_term_labels, names(list_of_term_labels), .f = ~paste0(".param_",.y,"(",.x,")"))
# }

merge_formulae <- function(forms, data) {
  list_of_term_labels <- purrr::map(purrr::map(forms, terms, data=data), attr, 'term.labels')
  form_out <- reformulate(unique(purrr::flatten_chr(list_of_term_labels)))
  lazyeval::f_lhs(form_out) <- lazyeval::f_lhs(forms[[1]])
  environment(form_out) <- environment(forms[[1]])
  form_out
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
      spline_dist = TRUE,
      location = 'gamma0',
      pars = paste0('gamma', 0:k),
      transforms = as.list(structure(names=pars, c('identity', 'log', rep('identity',k-1)))),
      inverse_transforms = as.list(structure(names=pars, c('identity', 'exp', rep('identity',k-1)))),
      pars_real = paste0(ifelse(transforms=="identity","",paste0(transforms,"_")),pars)
      ,
      distribution_function_factory = function(all_log_knots, scaling_factor = 50) {

        cdf_function <- function(q, gamma, lower.tail = TRUE, log.p =FALSE, b = NULL) {
          invalid_idx <- which(q<=0)
          q[invalid_idx] <- NA_real_
          log_cumu_haz <- roy_parm_log_cumu_haz(log_time = log(q), gamma = gamma, all_log_knots = all_log_knots, scaling_factor = scaling_factor, b = b)
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
    dist_infos[['roy_parm_splines']]$default_priors <-
      with(dist_infos[['roy_parm_splines']],
           function(model_mats, y = NULL) {
             if (is.null(y)) gamma0_init <- 0
             else gamma0_init <- -log(mean(y[,3]))
             model_mats <- model_mats[pars]
             names(model_mats) <- pars_real
             priors <- purrr::map_df(
               model_mats, .id = 'parameter',
               function(mat) {
                 data_frame(term = colnames(mat),
                            family = "normal", #double_exponential
                            location = rep(0, ncol(mat)),
                            spread = rep(1, ncol(mat)) )
               })
             priors$location <- ifelse(priors$parameter == 'gamma0' & priors$term == "(Intercept)",
                                 yes = gamma0_init,
                                 no = priors$location)
             priors
           })
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
#' @param f The actual factor vector
#'
#' @return Contrast-matrix
#' @export
contr.full <- function(n, contrasts=TRUE, sparse=FALSE, f = NULL) {
  if (!contrasts)
    stop(call. = FALSE, 'Please report this error to the package maintainer.')
  out <- contr.treatment(n = n, contrasts = contrasts, sparse = sparse, base = 1)
  to_add <- matrix(nrow = nrow(out), ncol = 1, data = c(1,rep(0,nrow(out)-1)),
                   dimnames = list(row.names(out), rownames(out)[1]))
  out <- cbind(out,to_add)
  if (!is.null(f)) {
    centers <- setNames(purrr::map_dbl(levels(f),~mean(f==.x,na.rm=TRUE)), nm = levels(f))
    out <- out - matrix(centers, nrow = length(centers), ncol = length(centers))
  }
  return(out)
}

#' Set prior
#'
#' @param parameter The parameter
#' @param terms The terms. Can be a character vector, or the result of \code{dplyr::vars}, which can
#'   be convenient for setting the prior on all but some terms.
#' @param location A value for location on the distribution.
#' @param spread A value for the spread (scale) on the distribution
#'
#' @return An object of class 'prior', to be (optionally inserted into a list with other priors and)
#'   passed to `priors` arg in `survreg_map`
#' @export
set_prior <- function(parameter = NULL, terms, location = NULL, spread = NULL) {
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
    browser()
    if (!is.null(location))
      prior_df$location[prior_df$term%in%terms_to_modify & param_lgl] <- location
    if (!is.null(spread))
      prior_df$spread[prior_df$term%in%terms_to_modify & param_lgl] <- spread
    prior_df
  }
  class(out) <- c('prior',class(out))
  out
}
