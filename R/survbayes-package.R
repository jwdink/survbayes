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
#'   with columns 'parameter', 'term', 'family', 'location', 'scale'. An optional logical column
#'   'fixed' can be used to specify parameters that should be fixed at their initial value and not
#'   updated during optimization.
#' @param standardize_x If this is a logical, it controls whether the model-matrix should be
#'   centered and scaled before being passed for optimization. If this is a data.frame with columns
#'   'parameter', 'term', 'center' and 'scale', allows the user to pass their own centering/scaling
#'   (or for c/s to be used on calls to `update`). Centering and scaling is strongly recommended,
#'   because the default inits (and priors) are designed to work with centered and scaled data.
#' @param optim_args Arguments to pass to \code{stats::optim}
#'
#' @return
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
    optim_args = list(method = "BFGS", control = list(trace=as.integer(interactive()),maxit=250))
  ) {

    the_call <- match.call()

    ## formula, distribution:
    if (is.null(lazyeval::f_lhs(formula))) stop(call. = FALSE, "`formula` has no left-hand side.")
    dist_info <- get_dist_info(distribution, k = max(length(anc),1) )
    forms <- standardize_formulae(formula, anc, dist_info)

    ## model-componenets:
    model_components <-
      prepare_model_components(forms, data, dist_info, standardize_x, na.action, xlev = NULL, drop.unused.levels = FALSE)
    the_call$standardize_x <- model_components$standardize_x

    ## prior:
    # [to do]
    the_call$priors <- priors

    ## fit:
    # [to do]

  }

prepare_model_components <- function(forms, data, dist_info, standardize_x, na.action, xlev, drop.unused.levels) {
  formula_merged <- merge_formulae(forms, data)
  model_frame_merged <-
    model.frame(formula_merged, data = data, na.action = na.action, drop.unused.levels = drop.unused.levels, xlev = xlev)
  model_matrix_merged <-
    model.matrix(terms(model_frame_merged), data = model_frame_merged, xlev = xlev)
  term_mapping <- get_terms_mapping(formula_merged, data)

  purrr::map(forms, function(this_form) {
    browser()
  })
}

#' Get mapping from column names to model-terms and vice versa
#'
#' @param formula The formula to be passed to model.frame
#' @param data The data.frame
#' @param ... Arguments to be passed to \code{stats::model.matrix}
#'
#' @return
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

    ## flip the mm -> original mapping to give original -> mm mapping:
    from_od_to_mm <- flip_mapping(from_mm_to_od)

    flip_mapping(from_tl_to_mf)
    browser()
  }

  #
  list(
    original_data = list(
      term_labels = ,
      model_frame = ,
      model_matrix = from_od_to_mm
    ),
    term_labels = list(
      original_data = ,
      model_frame = from_tl_to_mf,
      model_matrix =
    ),
    model_frame = list(
      original_data = from_mf_to_od,
      term_labels = ,
      model_matrix =
    ),
    model_mat = list(
      original_data = from_mm_to_od,
      term_labels = ,
      model_frame =
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
