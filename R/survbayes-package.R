#' survbayes package
#'
#' @import lazyeval
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import survival
#'
#' @importFrom tibble rownames_to_column lst enframe
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
#'             contr.treatment formula setNames family getCall update.formula
#' @importFrom grDevices dev.off pdf
#' @importFrom glue glue
#' @importFrom pryr dots
#'
#' @docType package
#' @name survbayes
NULL
#> NULL



#' Get Cross-Validated Log-Likelihood
#'
#' @param object An object to be cross-validated
#'
#' @return A vector of log-liklihoods, one for each fold
#' @export
crossv_loglik <- function(object, ...)  UseMethod('crossv_loglik')


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

#' Standardize formula arguments
#'
#' @param formula The main formula
#' @param anc List of others
#' @param dist_info Information about the distribution
#'
#' @return A list of formulas, now matched to the distribution
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

#' Merge a list of formulae
#'
#' @param forms List of formulae/formulas
#' @param data Data.frame
#'
#' @return Single formula, with environment(forms[[1]])
merge_formulae <- function(forms, data) {
  list_of_term_labels <- purrr::map(purrr::map(forms, terms, data=data), attr, 'term.labels')
  term_labels_unique <- unique(purrr::flatten_chr(list_of_term_labels))
  if (length(term_labels_unique)>0) form_out <- reformulate(term_labels_unique)
  else form_out <- ~1
  lazyeval::f_lhs(form_out) <- lazyeval::f_lhs(forms[[1]])
  environment(form_out) <- environment(forms[[1]])
  form_out
}

#' Get information about a survival distribution
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
                            spread = rep(1, ncol(mat)),
                            fixed = FALSE)
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
#' @param family The family for the prior distribution
#' @param location A value for location on the distribution.
#' @param spread A value for the spread (scale) on the distribution
#' @param fixed In the fitting process, should this parameter be fixed at its prior value, or updated?
#'
#' @return An object of class 'prior', to be (optionally inserted into a list with other priors and)
#'   passed to `priors` arg in `survreg_map`
#' @export
set_prior <- function(parameter = NULL, terms, family = NULL, location = NULL, spread = NULL, fixed=NULL) {
  stopifnot(is.character(terms) || inherits(terms, "col_list") || inherits(terms,"quosures") )

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

    if (!is.null(location)) {
      stopifnot(length(location)==1)
      prior_df$location[prior_df$term%in%terms_to_modify & param_lgl] <- location
    }
    if (!is.null(spread)) {
      stopifnot(length(spread)==1)
      prior_df$spread[prior_df$term%in%terms_to_modify & param_lgl] <- spread
    }
    if (!is.null(family)) {
      stopifnot(length(family)==1)
      prior_df$family[prior_df$term%in%terms_to_modify & param_lgl] <- family
    }
    if (!is.null(fixed)) {
      stopifnot(length(fixed)==1)
      prior_df$fixed[prior_df$term%in%terms_to_modify & param_lgl] <- fixed
    }
    prior_df
  }
  class(out) <- c('prior',class(out))
  out
}
