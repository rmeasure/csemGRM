#' Compute Analytical Conditional Standard Error of Measurement for GRM
#'
#' Derives the analytical conditional standard error of measurement (CSEM)
#' at each point along the theta continuum using the GRM test information
#' function. The analytical CSEM is the reciprocal of the square root of
#' the test information function: CSEM(theta) = 1 / sqrt(I(theta)).
#'
#' @param fit An object of class \code{"csemgrm_fit"} from
#'   \code{csem_estimate()}.
#'
#' @details
#' The test information function I(theta) in GRM summarizes the total
#' amount of measurement precision at each level of the latent trait.
#' The analytical CSEM derived from I(theta) provides a continuous
#' profile of measurement precision across the theta continuum, showing
#' where the scale measures most precisely (low CSEM) and where
#' precision is poorest (high CSEM).
#'
#' For attitude and character measurement instruments, the analytical
#' CSEM profile reveals whether the instrument provides differential
#' precision across trait levels — a critical consideration for
#' instrument development and scale refinement.
#'
#' @return A list of class \code{"csemgrm_analytical"} containing:
#'   \item{theta_grid}{Numeric vector of theta grid points.}
#'   \item{test_info}{Test information at each theta point.}
#'   \item{csem}{Analytical CSEM at each theta point.}
#'   \item{csem_df}{Data frame of theta, test_info, and csem.}
#'   \item{theta_min_csem}{Theta value where CSEM is minimized.}
#'   \item{min_csem}{Minimum CSEM value.}
#'   \item{theta_max_info}{Theta value where test information is maximized.}
#'   \item{max_info}{Maximum test information value.}
#'
#' @references
#'   Samejima, F. (1994). Estimation of reliability coefficients using the
#'   test information function and its modifications.
#'   \emph{Applied Psychological Measurement, 18}(3), 229-244.
#'   \doi{10.1177/014662169401800304}
#'
#'   Lord, F. M. (1980). \emph{Applications of item response theory to
#'   practical testing problems.} Lawrence Erlbaum Associates.
#'
#' @examples
#' data(attitude_data)
#' items <- attitude_data[, 1:10]
#' fit   <- csem_estimate(items)
#' csem  <- csem_analytical(fit)
#' print(csem)
#'
#' @export
csem_analytical <- function(fit) {

  if (!inherits(fit, "csemgrm_fit"))
    stop("'fit' must be a 'csemgrm_fit' object from csem_estimate().")

  theta <- fit$theta_grid
  info  <- fit$test_info

  # Avoid division by zero
  info_safe <- pmax(info, 1e-6)

  # CSEM = 1 / sqrt(I(theta))
  csem <- 1 / sqrt(info_safe)

  csem_df <- data.frame(
    theta     = round(theta, 4),
    test_info = round(info,  4),
    csem      = round(csem,  4)
  )

  idx_min  <- which.min(csem)
  idx_max  <- which.max(info)

  result <- list(
    theta_grid    = theta,
    test_info     = info,
    csem          = csem,
    csem_df       = csem_df,
    theta_min_csem = round(theta[idx_min], 4),
    min_csem       = round(csem[idx_min],  4),
    theta_max_info = round(theta[idx_max], 4),
    max_info       = round(info[idx_max],  4),
    n_items        = fit$n_items,
    n_persons      = fit$n_persons
  )
  class(result) <- "csemgrm_analytical"
  return(result)
}


#' @export
print.csemgrm_analytical <- function(x, ...) {
  cat("=== Analytical CSEM — Graded Response Model ===\n\n")
  cat(sprintf("Peak information:  %.4f at theta = %.4f\n",
              x$max_info, x$theta_max_info))
  cat(sprintf("Minimum CSEM:      %.4f at theta = %.4f\n",
              x$min_csem, x$theta_min_csem))
  cat(sprintf("CSEM range:        %.4f to %.4f\n\n",
              min(x$csem), max(x$csem)))
  cat("CSEM at selected theta values:\n\n")
  idx <- c(which.min(abs(x$theta_grid - (-2))),
           which.min(abs(x$theta_grid - (-1))),
           which.min(abs(x$theta_grid -   0)),
           which.min(abs(x$theta_grid -   1)),
           which.min(abs(x$theta_grid -   2)))
  print(x$csem_df[idx, ], row.names = FALSE)
  invisible(x)
}
