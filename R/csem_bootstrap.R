#' Estimate Bootstrap Conditional Standard Error of Measurement for GRM
#'
#' Estimates the conditional standard error of measurement (CSEM) and its
#' confidence band via parametric bootstrap resampling. In each bootstrap
#' replication, item responses are simulated from the fitted GRM parameters,
#' the model is re-estimated, and the CSEM profile is computed. The resulting
#' distribution of CSEM profiles provides bootstrap standard errors and
#' confidence intervals for the CSEM at each theta point.
#'
#' @param fit An object of class \code{"csemgrm_fit"} from
#'   \code{csem_estimate()}.
#' @param n_boot Integer. Number of bootstrap replications. Default: 100.
#'   Use at least 500 for stable confidence intervals in practice.
#' @param conf_level Numeric. Confidence level for bootstrap CI.
#'   Default: 0.95.
#' @param seed Integer. Random seed for reproducibility. Default: 42.
#' @param verbose Logical. Print bootstrap progress. Default \code{FALSE}.
#'
#' @return A list of class \code{"csemgrm_bootstrap"} containing:
#'   \item{theta_grid}{Numeric vector of theta grid points.}
#'   \item{csem_analytical}{Analytical CSEM from the observed fit.}
#'   \item{csem_boot_mean}{Bootstrap mean CSEM at each theta point.}
#'   \item{csem_boot_se}{Bootstrap SE of CSEM at each theta point.}
#'   \item{csem_boot_lower}{Lower confidence bound.}
#'   \item{csem_boot_upper}{Upper confidence bound.}
#'   \item{csem_df}{Data frame combining all CSEM estimates.}
#'   \item{n_boot}{Number of successful bootstrap replications.}
#'   \item{conf_level}{Confidence level used.}
#'
#' @references
#'   Efron, B., & Tibshirani, R. J. (1993).
#'   \emph{An introduction to the bootstrap.} Chapman & Hall.
#'
#'   Samejima, F. (1994). Estimation of reliability coefficients using the
#'   test information function and its modifications.
#'   \emph{Applied Psychological Measurement, 18}(3), 229-244.
#'   \doi{10.1177/014662169401800304}
#'
#' @examples
#' data(attitude_data)
#' items <- attitude_data[, 1:10]
#' fit   <- csem_estimate(items)
#' boot  <- csem_bootstrap(fit, n_boot = 50, seed = 42)
#' print(boot)
#'
#' @export
csem_bootstrap <- function(fit,
                           n_boot     = 100,
                           conf_level = 0.95,
                           seed       = 42,
                           verbose    = FALSE) {

  if (!inherits(fit, "csemgrm_fit"))
    stop("'fit' must be a 'csemgrm_fit' object from csem_estimate().")
  if (n_boot < 10)
    stop("'n_boot' must be at least 10.")

  set.seed(seed)

  theta_grid <- fit$theta_grid
  theta_mat  <- matrix(theta_grid, ncol = 1)
  n_persons  <- fit$n_persons
  n_items    <- fit$n_items
  alpha      <- 1 - conf_level

  # Analytical CSEM from observed fit
  info_obs  <- fit$test_info
  csem_obs  <- 1 / sqrt(pmax(info_obs, 1e-6))

  # Bootstrap loop
  csem_boot_mat <- matrix(NA, nrow = n_boot, ncol = length(theta_grid))
  success       <- 0

  if (verbose) message("Running ", n_boot, " bootstrap replications...")

  for (b in seq_len(n_boot)) {

    # Simulate data from fitted GRM using mirt::simdata
    sim_data <- tryCatch({
      # Extract parameters directly from fitted model
      params <- mirt::coef(fit$fit, simplify = TRUE, IRTpars = TRUE)$items
      a_mat  <- matrix(params[, "a"], ncol = 1)
      # Get d parameters (mirt native format) for simulation
      params_d <- mirt::coef(fit$fit, simplify = TRUE)$items
      d_cols   <- grep("^d", colnames(params_d))
      d_mat    <- as.matrix(params_d[, d_cols, drop = FALSE])
      suppressMessages(
        mirt::simdata(
          a        = a_mat,
          d        = d_mat,
          N        = n_persons,
          itemtype = "graded",
          Theta    = matrix(rnorm(n_persons), ncol = 1)
        )
      )
    }, error = function(e) NULL)
    if (is.null(sim_data)) next
    # Convert to data.frame with responses starting from 1
    sim_data <- as.data.frame(sim_data + 1L)

    # Re-fit GRM
    fit_b <- tryCatch(
      suppressMessages(
        mirt::mirt(sim_data, 1, itemtype = "graded",
                   verbose = FALSE, SE = FALSE)
      ),
      error = function(e) NULL
    )
    if (is.null(fit_b)) next

    # Compute CSEM
    info_b <- tryCatch(
      as.numeric(mirt::testinfo(fit_b, theta_mat)),
      error = function(e) NULL
    )
    if (is.null(info_b)) next

    csem_b <- 1 / sqrt(pmax(info_b, 1e-6))
    csem_boot_mat[success + 1, ] <- csem_b
    success <- success + 1

    if (verbose && success %% 10 == 0)
      message("  Completed ", success, " replications...")
  }

  if (success < 5)
    stop("Fewer than 5 bootstrap replications succeeded. ",
         "Check that the model converges on simulated data.")

  # Trim failed rows
  csem_boot_mat <- csem_boot_mat[seq_len(success), ]

  # Summary statistics
  csem_mean  <- colMeans(csem_boot_mat,  na.rm = TRUE)
  csem_se    <- apply(csem_boot_mat, 2, sd, na.rm = TRUE)
  csem_lower <- apply(csem_boot_mat, 2,
                      quantile, probs = alpha/2,   na.rm = TRUE)
  csem_upper <- apply(csem_boot_mat, 2,
                      quantile, probs = 1 - alpha/2, na.rm = TRUE)

  csem_df <- data.frame(
    theta          = round(theta_grid, 4),
    csem_analytical = round(csem_obs,   4),
    csem_boot_mean = round(csem_mean,  4),
    csem_boot_se   = round(csem_se,    4),
    csem_lower     = round(csem_lower, 4),
    csem_upper     = round(csem_upper, 4)
  )

  result <- list(
    theta_grid      = theta_grid,
    csem_analytical = csem_obs,
    csem_boot_mean  = csem_mean,
    csem_boot_se    = csem_se,
    csem_boot_lower = csem_lower,
    csem_boot_upper = csem_upper,
    csem_df         = csem_df,
    n_boot          = success,
    conf_level      = conf_level
  )
  class(result) <- "csemgrm_bootstrap"
  return(result)
}


#' @export
print.csemgrm_bootstrap <- function(x, ...) {
  cat("=== Bootstrap CSEM — Graded Response Model ===\n\n")
  cat(sprintf("Successful replications: %d\n", x$n_boot))
  cat(sprintf("Confidence level: %.0f%%\n\n", x$conf_level * 100))
  cat("CSEM at selected theta values:\n\n")
  idx <- c(which.min(abs(x$theta_grid - (-2))),
           which.min(abs(x$theta_grid - (-1))),
           which.min(abs(x$theta_grid -   0)),
           which.min(abs(x$theta_grid -   1)),
           which.min(abs(x$theta_grid -   2)))
  print(x$csem_df[idx, ], row.names = FALSE)
  invisible(x)
}
