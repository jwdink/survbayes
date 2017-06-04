#' Get the log-cumulative hazard for a set of log-times given a Royston/Parmar-spline-survival
#' distribution
#'
#' @param log_time log(time)
#' @param gamma Matrix of coefficients to be multiplied elementwise by basis
#' @param all_knots Knots, including boundary knots, on the scale of log(time)
#' @param scaling_factor A constant which gammas > 1 will be divided by, so that coefficients are on
#'   similar scales; helpful for optimization/MCMC
#' @param b computing the basis can be costly. if you have already have the basis, you can pass it
#'   to this fxn, instead of passing log-time.
#'
#' @return The log-cumulative hazard.
#' @export
roy_parm_log_cumu_haz <- function(log_time, gamma, all_knots, scaling_factor, b= NULL) {
  if (is.null(b))
    b <- roy_parm_basis(x = log_time, all_knots = all_knots, scaling_factor = scaling_factor)
  if (!is.matrix(gamma) || nrow(gamma)==1) {
    # check if same length as b. if so, replicate into matrix
    if ( length(gamma) != ncol(b) )
      stop(call. = FALSE, "If `gamma` is a vector, must have length equal to the number of num-knots + 2.")
    gamma <- matrix(data = rep(gamma, nrow(b)), nrow=nrow(b), byrow = TRUE)
  } else if (nrow(gamma)!=length(log_time) | ncol(gamma)!=length(all_knots) ) {
    stop(call. = FALSE,
         "If `gamma` is a matrix, must have num-rows equal to `length(log_time)`, and num-cols equal to num-knots + 2")
  }
  rowSums(b * gamma)
}

#' Get the basis matrix for the Royston/Parmar-spline parameterization
#'
#' @param x Log-time,
#' @param all_knots Knots, including boundary knots, on the scale of log(time)
#' @param scaling_factor A constant which gammas > 1 will be divided by, so that coefficients are on
#'   similar scales; helpful for optimization/MCMC
#'
#' @return The basis matrix
#' @export
#'
roy_parm_basis <- function(x, all_knots, scaling_factor) {
  nx <- length(x)
  nk <- length(all_knots)
  b <- matrix(nrow = nx, ncol = nk)
  knots <- matrix(rep(all_knots, nx), byrow = TRUE, ncol = nk)
  if (nk>0) {
    b[,1] <- rep(1,nx)
    b[,2] <- x
  }
  if (nk>2) {
    lam <- (knots[, nk] - knots)/(knots[, nk] - knots[, 1])
    for (j in 1:(nk - 2)) {
      part1 <- (x - knots[,j+1])^3
      part2 <- lam[,j+1] * (x - knots[,1])^3
      part3 <- (x - knots[,nk])^3
      b[,j+2] <- ( pmax(0, part1) - pmax(0, part2) - pmax(0, part3) ) / scaling_factor
    }
  }
  return(b)
}

#' Get Stan code for helper functions for Royston/Parmar-splines
#'
#' @return a character string of model code for the stan language
get_roy_parm_function_block_code <- function() {
  "
functions {

    vector row_sums(matrix X) {
      vector[rows(X)] s ;
      s = X * rep_vector(1, cols(X)) ;
      return s ;
    }

    // basis fxn for knots. makes a matrix where each col captures part of the curvature (like stats::poly)
    // x is a vector of (log) times; knots is a vector of knots (on the log scale), including boundary knots.
    // scaling factor is a number that ensures the columns 3 through num-knots are on roughly the same scale
    // as the first two columns. I'll be using median(time)^3
    matrix basis(vector x, vector knots, real scaling_factor) {
      int nx; // number of obs
      int nk; // number of knots
      matrix[num_elements(x),num_elements(knots)] b; // resulting basis matrix
      vector[num_elements(x)] lam;
      matrix[num_elements(x),num_elements(knots)-2] part1;
      matrix[num_elements(x),num_elements(knots)-2] part2;
      matrix[num_elements(x),num_elements(knots)-2] part3;

      nx = num_elements(x);
      nk = num_elements(knots);
      if (nk>0) {
        b[,1] = rep_vector(1,nx);
        b[,2] = x;
      }
      if (nk>2) {
      for (j in 1:nk) {
        lam[j] = (knots[nk] - knots[j]) / (knots[nk] - knots[1]);
      }
      for (j in 1:(nk-2)) { // for each of the remaining columns...
        for (r in 1:nx) { // for each row...
          part1[r,j] = pow( x[r] - knots[j+1] , 3);                   // pmax(x - knots[,j+1], 0)^3
          part2[r,j] = lam[j+1] * pow( x[r] - knots[1] , 3);        // lam[,j+1]*pmax(x - knots[,1], 0)^3
          part3[r,j] = (1 - lam[j+1]) * pow( x[r] - knots[nk] , 3); // (1 - lam[,j+1])*pmax(x - knots[,nk], 0)^3
          b[r,j+2] = ( ( (part1[r,j]<0) ? 0 : part1[r,j]) - ( (part2[r,j]<0) ? 0 : part2[r,j]) - ( (part3[r,j]<0) ? 0 : part3[r,j]) ) / scaling_factor;
          }
        }
      }
      return b;
    }

    // the spline log cumulative hazard
    vector spline_log_cumu_haz(vector t, vector knots, matrix gamma, real scaling_factor) {
      matrix[num_elements(t),num_elements(knots)] b;
      b = basis(log(t+.001), knots, scaling_factor);
      return row_sums(b .* gamma);
    }

    // the spline survival (ccdf) function. returns probability for each element of t
    vector spline_surv(vector t, vector knots, matrix gamma, real scaling_factor) {
      vector[num_elements(t)] log_cumu_haz;
      vector[num_elements(t)] cumu_haz;
      vector[num_elements(t)] surv;
      log_cumu_haz = spline_log_cumu_haz(t,knots,gamma, scaling_factor);
      cumu_haz = exp(log_cumu_haz);
      surv = exp(-cumu_haz);
      return surv;
    }

    real spline_interval_obs_lpdf(matrix x, vector knots, matrix gamma, real scaling_factor) {
      vector[rows(x)] prob;
      prob = (spline_surv(x[,2], knots, gamma, scaling_factor) - spline_surv(x[,3], knots, gamma, scaling_factor)) ./ spline_surv(x[,1], knots, gamma, scaling_factor);
      return sum(log(prob));
    }
    real spline_interval_cens_lpdf(matrix x, vector knots, matrix gamma, real scaling_factor) {
      vector[rows(x)] prob;
      prob = spline_surv(x[,3], knots, gamma, scaling_factor) ./ spline_surv(x[,1], knots, gamma, scaling_factor);
      return sum(log(prob));
    }
  }
  "
  }


#' Generate Stan model-code for fitting a model with the Royston/Parmar-spline-survival distribution
#'
#' @param model_mats List of model-matrices
#' @param priors A dataframe specifying priors
#'
#' @return a character string of model code for the stan language
make_roy_parm_spline_stan_model_code <- function(model_mats, priors) {

  template <-
           "
          data {{
             int<lower=1> N_obs;
             int obs_idx[N_obs];
             int<lower=1> N_cens;
             int cens_idx[N_cens];

             int<lower=1> k;
             vector[k] knots;
             real scaling_factor;

             int<lower=1> N;
             matrix[N, 4] Y;

             int<lower=1> num_covs_gamma0;
             matrix[N, num_covs_gamma0] X_gamma0;

             int<lower=1> num_covs_gamma1;
             matrix[N, num_covs_gamma1] X_gamma1;

             {OTHER_GAMMA_MODEL_MAT_INITS}
           }}
          transformed data {{
           }}
          parameters {{
            vector[num_covs_gamma0] gamma0;
            vector[num_covs_gamma1] log_gamma1;
            {OTHER_GAMMA_PARAM_INITS}
           }}
          transformed parameters {{
            vector[num_covs_gamma1] gamma1;

            gamma1 = exp(log_gamma1);
           }}
          model {{
            matrix[N,k] eta_gamma; // k gamma estimates for each row of the dataset

             {ALL_GAMMA_PRIORS}

             // matrix multiplication for value of parameter per-row:
             eta_gamma[,1] = X_gamma0 * gamma0;
             eta_gamma[,2] = exp(X_gamma1 * log_gamma1);
             {OTHER_GAMMA_MAT_MULTI}

             target += spline_interval_cens_lpdf(Y[cens_idx,] | knots, eta_gamma[cens_idx,], scaling_factor) ;
             target += spline_interval_obs_lpdf(Y[obs_idx,] | knots, eta_gamma[obs_idx,], scaling_factor) ;

             }}
          generated quantities {{
             // For model comparison, we'll want to keep the likelihood contribution of each point
             vector[N_obs+N_cens] log_lik;
             matrix[N,k] eta_gamma; // k gamma estimates for each row of the dataset

             // matrix multiplication for value of parameter per-row:
             eta_gamma[,1] = X_gamma0 * gamma0;
             eta_gamma[,2] = exp(X_gamma1 * log_gamma1);
             {OTHER_GAMMA_MAT_MULTI}

             log_lik[obs_idx] = log( (spline_surv(Y[obs_idx,2], knots, eta_gamma[obs_idx,], scaling_factor) - spline_surv(Y[obs_idx,3], knots, eta_gamma[obs_idx,], scaling_factor)) ./ spline_surv(Y[obs_idx,1], knots, eta_gamma[obs_idx,], scaling_factor) );
             log_lik[cens_idx] = log( spline_surv(Y[cens_idx,3], knots, eta_gamma[cens_idx,], scaling_factor) ./ spline_surv(Y[cens_idx,1], knots, eta_gamma[cens_idx,], scaling_factor) );

           }}
           "

  other_gamma_nums <- seq_along(model_mats)[-(1:2)]
  if (length(other_gamma_nums)>0) other_gamma_nums <- other_gamma_nums - 1

  prior_model_code <- purrr::map2(.x = split(priors, priors$parameter)[names(model_mats)],
              .y = names(model_mats),
              .f = function(this_param_prior_df, param_name) {
                prior_char_vec <- paste0("normal(", this_param_prior_df$mu, ",", this_param_prior_df$sigma, ")")
                if (n_distinct(prior_char_vec)==1) paste0(param_name, " ~ ", unique(prior_char_vec), ";\n")
                else paste(param_name, "[", seq_along(prior_char_vec), "] ~", prior_char_vec, ";\n")
              })
  prior_model_code <- prior_model_code[!purrr::map_lgl(prior_model_code, is.null)]

  ALL_GAMMA_PRIORS <- paste0(collapse="",purrr::flatten_chr(prior_model_code))

  if (length(other_gamma_nums)>0) {
    OTHER_GAMMA_MODEL_MAT_INITS <- glue::glue(
      "\tint<lower=1> num_covs_gamma{other_gamma_nums};\n",
        "\tmatrix[N, num_covs_gamma{other_gamma_nums}] X_gamma{other_gamma_nums};\n"
      )

    OTHER_GAMMA_PARAM_INITS <- glue::glue("vector[num_covs_gamma{other_gamma_nums}] gamma{other_gamma_nums};\n")
    OTHER_GAMMA_MAT_MULTI <- glue::glue("eta_gamma[,{other_gamma_nums+1}] = X_gamma{other_gamma_nums} * gamma{other_gamma_nums};\n")

  } else {
    OTHER_GAMMA_MODEL_MAT_INITS <- ""
    OTHER_GAMMA_PARAM_INITS <- ""
    OTHER_GAMMA_MAT_MULTI <- ""
  }

  return(paste0(get_roy_parm_function_block_code(), glue::glue(template)))

  }
