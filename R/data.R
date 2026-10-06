#' Simulated Attitude Scale Dataset
#'
#' A simulated dataset representative of a polytomous attitude measurement
#' instrument with 10 items on a 5-point Likert scale. Item parameters are
#' calibrated to represent a well-functioning attitude scale with varying
#' discrimination across items.
#'
#' @format A data frame with 500 rows and 11 variables:
#'   \describe{
#'     \item{A01--A10}{Item responses (1-5 Likert scale)}
#'     \item{group}{Group variable: "Group_1" or "Group_2"}
#'   }
#'
#' @examples
#' data(attitude_data)
#' head(attitude_data)
"attitude_data"


#' Simulated Character Scale Dataset
#'
#' A simulated dataset representative of a polytomous character assessment
#' instrument with 12 items on a 5-point Likert scale measuring three
#' character dimensions: cognitive (C01-C04), affective (A01-A04), and
#' behavioral (B01-B04).
#'
#' @format A data frame with 400 rows and 13 variables:
#'   \describe{
#'     \item{C01--C04}{Cognitive character items (1-5)}
#'     \item{A01--A04}{Affective character items (1-5)}
#'     \item{B01--B04}{Behavioral character items (1-5)}
#'     \item{group}{Group variable: "Group_1" or "Group_2"}
#'   }
#'
#' @examples
#' data(character_data)
#' head(character_data)
"character_data"


# ── Generate built-in datasets ────────────────────────────────────────────────
.generate_csem_datasets <- function() {

  set.seed(2025)

  sim_grm <- function(theta, a, b_mat) {
    n <- length(theta); k <- ncol(b_mat)+1
    out <- matrix(0L, n, nrow(b_mat))
    for (j in seq_len(nrow(b_mat))) {
      cum_p <- cbind(1,
        sapply(b_mat[j,], function(bk) 1/(1+exp(-a[j]*(theta-bk)))),
        0)
      cat_p <- t(apply(cum_p, 1, function(r) pmax(diff(-r), 0)))
      cat_p <- cat_p/rowSums(cat_p)
      out[,j] <- apply(cat_p, 1, function(p) sample(seq_len(k),1,prob=p))
    }
    as.data.frame(out)
  }

  # ── Attitude dataset (10 items, 5-point, N=500) ─────────────────────────────
  N1     <- 500
  theta1 <- rnorm(N1, 0, 1)
  a1 <- c(1.8,1.5,2.0,1.6,1.3,1.7,1.9,1.4,1.6,1.8)
  b1 <- matrix(c(
    -2.0,-1.0, 0.0, 1.0,
    -1.8,-0.8, 0.2, 1.2,
    -2.2,-1.2,-0.2, 0.8,
    -1.6,-0.6, 0.4, 1.4,
    -1.5,-0.5, 0.5, 1.5,
    -1.7,-0.7, 0.3, 1.3,
    -2.1,-1.1,-0.1, 0.9,
    -1.4,-0.4, 0.6, 1.6,
    -1.9,-0.9, 0.1, 1.1,
    -1.6,-0.6, 0.4, 1.4
  ), nrow=10, byrow=TRUE)

  att <- sim_grm(theta1, a1, b1)
  names(att) <- paste0("A", sprintf("%02d", 1:10))
  att$group  <- ifelse(theta1 > 0, "Group_1", "Group_2")

  # ── Character dataset (12 items, 5-point, N=400) ─────────────────────────────
  N2     <- 400
  theta2 <- rnorm(N2, 0, 1)
  a2 <- c(1.7,1.5,1.9,1.6, 1.8,1.4,1.7,1.5, 1.6,1.8,1.5,1.7)
  b2 <- matrix(rep(c(-2.0,-1.0, 0.0, 1.0), 12), nrow=12, byrow=TRUE) +
        matrix(rnorm(12*4, 0, 0.3), nrow=12)

  char <- sim_grm(theta2, a2, b2)
  names(char) <- c(paste0("C", sprintf("%02d",1:4)),
                   paste0("A", sprintf("%02d",1:4)),
                   paste0("B", sprintf("%02d",1:4)))
  char$group  <- ifelse(theta2 > 0, "Group_1", "Group_2")

  list(attitude_data = att, character_data = char)
}
