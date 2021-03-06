% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roy-parm-distribution.R
\name{roy_parm_log_cumu_haz}
\alias{roy_parm_log_cumu_haz}
\title{Get the log-cumulative hazard for a set of log-times given a Royston/Parmar-spline-survival
distribution}
\usage{
roy_parm_log_cumu_haz(log_time, gamma, all_log_knots, scaling_factor,
  b = NULL)
}
\arguments{
\item{log_time}{log(time)}

\item{gamma}{Matrix of coefficients to be multiplied elementwise by basis}

\item{all_log_knots}{Knots, including boundary knots, on the scale of log(time)}

\item{scaling_factor}{A constant which gammas > 1 will be divided by, so that coefficients are on
similar scales; helpful for optimization/MCMC}

\item{b}{computing the basis can be costly. if you have already have the basis, you can pass it
to this fxn, instead of passing log-time.}
}
\value{
The log-cumulative hazard.
}
\description{
Get the log-cumulative hazard for a set of log-times given a Royston/Parmar-spline-survival
distribution
}
