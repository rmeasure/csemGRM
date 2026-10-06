#' Estimate GRM Parameters for CSEM Analysis
#'
#' Fits a Graded Response Model (GRM) to polytomous item response data and
#' returns item parameters, person theta estimates, and test information values
#' needed for conditional standard error of measurement (CSEM) computation.
#'
#' @param data A \code{data.frame} or \code{matrix} of item responses.
#'   Rows = persons, columns = items. Responses must be consecutive integers
#'   starting from 1 (e.g., 1-5 for a 5-point Likert scale). Remove any
#'   non-item columns before passing.
#' @param theta_range Numeric vector of length 2. Range of theta values for
#'   CSEM computation. Default: \code{c(-4, 4)}.
#' @param theta_points Integer. Number of theta grid points. Default: 100.
#' @param method Character. Theta estimation method: \code{"EAP"} (default),
#'   \code{"MAP"}, or \code{"ML"}.
#' @param verbose Logical. Print estimation progress. Default \code{FALSE}.
#'
#' @return A list of class \code{"csemgrm_fit"} containing:
#'   \item{fit}{The fitted \code{mirt} model object.}
#'   \item{item_params}{Data frame of item parameters (a and b).}
#'   \item{theta_grid}{Numeric vector of theta grid points.}
#'   \item{test_info}{Numeric vector of test information at each theta point.}
#'   \item{item_info}{Matrix of item information (items x theta points).}
#'   \item{theta_persons}{Numeric vector of person theta estimates.}
#'   \item{se_persons}{Numeric vector of SE(theta) for each person.}
#'   \item{reliability}{Marginal reliability.}
#'   \item{n_persons}{Number of persons.}
#'   \item{n_items}{Number of items.}
#'   \item{n_cats}{Number of response categories.}
#'   \item{theta_range}{Theta range used.}
#'
#' @references
#'   Samejima, F. (1969). Estimation of latent ability using a response
#'   pattern of graded scores. \emph{Psychometrika Monograph Supplement, 17}.
#'   \doi{10.1007/BF03372160}
#'
#'   Chalmers, R. P. (2012). mirt: A multidimensional item response theory
#'   package for the R environment. \emph{Journal of Statistical Software,
#'   48}(6), 1-29. \doi{10.18637/jss.v048.i06}
#'
#' @examples
#' data(attitude_data)
#' items <- attitude_data[, 1:10]
#' fit   <- csem_estimate(items)
#' print(fit)
#'
#' @export
csem_estimate <- function(data,
                          theta_range  = c(-4, 4),
                          theta_points = 100,
                          method       = "EAP",
                          verbose      = FALSE) {

  # ── Input validation ─────────────────────────────────────────────────────────
  if (!is.data.frame(data) && !is.matrix(data))
    stop("'data' must be a data.frame or matrix.")

  data <- as.data.frame(data)

  not_num <- !sapply(data, function(x) is.numeric(x) || is.integer(x))
  if (any(not_num))
    stop("All columns must be numeric. Non-numeric columns found: ",
         paste(names(data)[not_num], collapse = ", "))

  if (any(is.na(data)))
    stop("Missing values detected. Remove or impute before estimating.")

  n_persons <- nrow(data)
  n_items   <- ncol(data)
  n_cats    <- max(unlist(data))

  if (n_persons < 200)
    warning("Sample size (N = ", n_persons, ") is below the recommended ",
            "minimum of 200 for stable GRM parameter estimation.")
  if (n_items < 5)
    warning("Number of items (k = ", n_items, ") is below the recommended ",
            "minimum of 5.")

  # ── Fit GRM ──────────────────────────────────────────────────────────────────
  message("Fitting GRM to ", n_persons, " persons x ", n_items, " items...")

  fit <- tryCatch(
    mirt::mirt(data     = data,
               model    = 1,
               itemtype = "graded",
               SE       = TRUE,
               verbose  = verbose),
    error = function(e)
      stop("GRM estimation failed: ", e$message)
  )

  # ── Extract item parameters ──────────────────────────────────────────────────
  coefs <- tryCatch(
    as.data.frame(round(
      mirt::coef(fit, simplify = TRUE, IRTpars = TRUE)$items, 4)),
    error = function(e)
      as.data.frame(round(
        mirt::coef(fit, simplify = TRUE)$items, 4))
  )
  coefs$item <- rownames(coefs)
  coefs <- coefs[, c("item", setdiff(names(coefs), "item"))]
  rownames(coefs) <- NULL

  # ── Theta grid ───────────────────────────────────────────────────────────────
  theta_grid <- seq(theta_range[1], theta_range[2],
                    length.out = theta_points)
  theta_mat  <- matrix(theta_grid, ncol = 1)

  # ── Test and item information ─────────────────────────────────────────────────
  test_info <- as.numeric(mirt::testinfo(fit, theta_mat))
  item_info <- sapply(seq_len(n_items), function(j) {
    as.numeric(mirt::iteminfo(
      mirt::extract.item(fit, j), theta_mat))
  })
  colnames(item_info) <- names(data)

  # ── Person estimates ──────────────────────────────────────────────────────────
  fs       <- mirt::fscores(fit, method = method, full.scores.SE = TRUE)
  theta_p  <- as.numeric(fs[, 1])
  se_p     <- as.numeric(fs[, 2])

  # ── Marginal reliability ──────────────────────────────────────────────────────
  rel <- tryCatch(
    as.numeric(mirt::empirical_rxx(fs)),
    error = function(e) {
      vt  <- var(theta_p)
      mse <- mean(se_p^2)
      (vt - mse) / vt
    }
  )

  result <- list(
    fit          = fit,
    item_params  = coefs,
    theta_grid   = theta_grid,
    test_info    = test_info,
    item_info    = item_info,
    theta_persons = theta_p,
    se_persons   = se_p,
    reliability  = round(rel, 4),
    n_persons    = n_persons,
    n_items      = n_items,
    n_cats       = n_cats,
    theta_range  = theta_range
  )
  class(result) <- "csemgrm_fit"
  return(result)
}


#' @export
print.csemgrm_fit <- function(x, ...) {
  cat("=== GRM Estimation for CSEM Analysis ===\n\n")
  cat(sprintf("N persons:     %d\n", x$n_persons))
  cat(sprintf("N items:       %d\n", x$n_items))
  cat(sprintf("N categories:  %d\n", x$n_cats))
  cat(sprintf("Reliability:   %.4f\n\n", x$reliability))
  cat("Theta grid: from", x$theta_range[1], "to", x$theta_range[2], "\n")
  cat("Peak test information:", round(max(x$test_info), 3),
      "at theta =", round(x$theta_grid[which.max(x$test_info)], 3), "\n\n")
  cat("Item Parameters (first 6):\n")
  print(head(x$item_params), row.names = FALSE)
  invisible(x)
}
