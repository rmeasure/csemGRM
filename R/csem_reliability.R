#' Compute Conditional Reliability for GRM
#'
#' Derives the conditional reliability coefficient as a function of the
#' latent trait level theta, following Samejima (1994). The conditional
#' reliability at theta is defined as 1 - CSEM(theta)^2 / Var(theta),
#' where Var(theta) is the variance of the latent trait distribution.
#' This provides a more nuanced picture of reliability than a single
#' marginal coefficient, showing how reliability varies across trait levels.
#'
#' @param fit An object of class \code{"csemgrm_fit"} from
#'   \code{csem_estimate()}.
#' @param analytical An object of class \code{"csemgrm_analytical"} from
#'   \code{csem_analytical()}. If \code{NULL}, computed internally.
#' @param var_theta Numeric. Variance of the latent trait distribution.
#'   Default: \code{1.0} (standard normal assumption).
#'
#' @return A list of class \code{"csemgrm_reliability"} containing:
#'   \item{theta_grid}{Numeric vector of theta grid points.}
#'   \item{csem}{Analytical CSEM at each theta point.}
#'   \item{crel}{Conditional reliability at each theta point.}
#'   \item{rel_df}{Data frame of theta, csem, and crel.}
#'   \item{marginal_reliability}{Marginal reliability from GRM fit.}
#'   \item{mean_crel}{Mean conditional reliability across theta.}
#'
#' @references
#'   Samejima, F. (1994). Estimation of reliability coefficients using the
#'   test information function and its modifications.
#'   \emph{Applied Psychological Measurement, 18}(3), 229-244.
#'   \doi{10.1177/014662169401800304}
#'
#' @examples
#' data(attitude_data)
#' items <- attitude_data[, 1:10]
#' fit   <- csem_estimate(items)
#' rel   <- csem_reliability(fit)
#' print(rel)
#'
#' @export
csem_reliability <- function(fit, analytical = NULL, var_theta = 1.0) {

  if (!inherits(fit, "csemgrm_fit"))
    stop("'fit' must be a 'csemgrm_fit' object.")
  if (var_theta <= 0)
    stop("'var_theta' must be positive.")

  # Compute analytical CSEM if not provided
  if (is.null(analytical))
    analytical <- csem_analytical(fit)

  theta <- analytical$theta_grid
  csem  <- analytical$csem

  # Conditional reliability: 1 - CSEM^2 / Var(theta)
  crel <- pmax(1 - csem^2 / var_theta, 0)

  rel_df <- data.frame(
    theta = round(theta, 4),
    csem  = round(csem,  4),
    crel  = round(crel,  4)
  )

  result <- list(
    theta_grid           = theta,
    csem                 = csem,
    crel                 = crel,
    rel_df               = rel_df,
    marginal_reliability = fit$reliability,
    mean_crel            = round(mean(crel), 4),
    var_theta            = var_theta
  )
  class(result) <- "csemgrm_reliability"
  return(result)
}


#' @export
print.csemgrm_reliability <- function(x, ...) {
  cat("=== Conditional Reliability — Graded Response Model ===\n\n")
  cat(sprintf("Marginal reliability:       %.4f\n", x$marginal_reliability))
  cat(sprintf("Mean conditional reliability: %.4f\n\n", x$mean_crel))
  cat("Conditional reliability at selected theta values:\n\n")
  idx <- c(which.min(abs(x$theta_grid - (-2))),
           which.min(abs(x$theta_grid - (-1))),
           which.min(abs(x$theta_grid -   0)),
           which.min(abs(x$theta_grid -   1)),
           which.min(abs(x$theta_grid -   2)))
  print(x$rel_df[idx, ], row.names = FALSE)
  invisible(x)
}
