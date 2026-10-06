#' Visualize CSEM Profiles for GRM
#'
#' Produces visualizations of conditional standard error of measurement (CSEM)
#' profiles derived from a fitted Graded Response Model. Four plot types are
#' available: CSEM profile, test information function, conditional reliability,
#' and item information functions.
#'
#' @param x An object of class \code{"csemgrm_fit"},
#'   \code{"csemgrm_analytical"}, \code{"csemgrm_bootstrap"}, or
#'   \code{"csemgrm_reliability"}.
#' @param type Character. Plot type:
#'   \code{"csem"} (default) — CSEM profile across theta;
#'   \code{"info"} — test information function;
#'   \code{"reliability"} — conditional reliability profile;
#'   \code{"items"} — item information functions (requires
#'     \code{"csemgrm_fit"} object);
#'   \code{"bootstrap"} — CSEM with bootstrap CI (requires
#'     \code{"csemgrm_bootstrap"} object).
#' @param title Character. Custom plot title. Default: auto-generated.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' data(attitude_data)
#' items <- attitude_data[, 1:10]
#' fit   <- csem_estimate(items)
#' csem  <- csem_analytical(fit)
#' rel   <- csem_reliability(fit, csem)
#' csem_plot(csem, type = "csem")
#' csem_plot(csem, type = "info")
#' csem_plot(rel,  type = "reliability")
#' csem_plot(fit,  type = "items")
#'
#' @export
csem_plot <- function(x, type = "csem", title = NULL) {

  valid_types <- c("csem", "info", "reliability", "items", "bootstrap")
  if (!type %in% valid_types)
    stop("'type' must be one of: ", paste(valid_types, collapse=", "))

  # ── CSEM profile ─────────────────────────────────────────────────────────────
  if (type == "csem") {

    if (inherits(x, "csemgrm_analytical")) {
      df  <- x$csem_df
      ttl <- title %||% "Conditional Standard Error of Measurement (GRM)"
      p <- ggplot2::ggplot(df, ggplot2::aes(x=theta, y=csem)) +
        ggplot2::geom_line(color="#1a237e", linewidth=1.2) +
        ggplot2::geom_hline(yintercept=0.3, linetype="dashed",
                            color="#2e7d32", linewidth=0.8) +
        ggplot2::geom_hline(yintercept=0.5, linetype="dashed",
                            color="#e65100", linewidth=0.8) +
        ggplot2::annotate("text", x=max(df$theta)*0.85, y=0.30,
                          label="CSEM = 0.30", color="#2e7d32", size=3.2) +
        ggplot2::annotate("text", x=max(df$theta)*0.85, y=0.50,
                          label="CSEM = 0.50", color="#e65100", size=3.2) +
        ggplot2::labs(title=ttl,
                      subtitle="Dashed lines: CSEM = 0.30 (good) and 0.50 (poor)",
                      x=expression(theta~"(Latent Trait)"),
                      y="CSEM") +
        ggplot2::theme_minimal(base_size=13)

    } else if (inherits(x, "csemgrm_fit")) {
      anal <- csem_analytical(x)
      return(csem_plot(anal, type="csem", title=title))
    } else {
      stop("For type='csem', provide a 'csemgrm_fit' or 'csemgrm_analytical'.")
    }

  # ── Test information ──────────────────────────────────────────────────────────
  } else if (type == "info") {

    if (inherits(x, "csemgrm_analytical")) {
      df  <- x$csem_df
      ttl <- title %||% "Test Information Function (GRM)"
      p <- ggplot2::ggplot(df, ggplot2::aes(x=theta, y=test_info)) +
        ggplot2::geom_line(color="#1a237e", linewidth=1.2) +
        ggplot2::geom_vline(xintercept=x$theta_max_info,
                            linetype="dashed", color="#b71c1c", linewidth=0.8) +
        ggplot2::annotate("text",
                          x=x$theta_max_info + 0.15,
                          y=max(df$test_info)*0.5,
                          label=paste0("Peak: ", round(x$max_info,2),
                                       "\nat \u03B8=", x$theta_max_info),
                          color="#b71c1c", size=3.2, hjust=0) +
        ggplot2::labs(title=ttl,
                      x=expression(theta~"(Latent Trait)"),
                      y="Test Information I(\u03B8)") +
        ggplot2::theme_minimal(base_size=13)

    } else if (inherits(x, "csemgrm_fit")) {
      anal <- csem_analytical(x)
      return(csem_plot(anal, type="info", title=title))
    } else {
      stop("For type='info', provide a 'csemgrm_fit' or 'csemgrm_analytical'.")
    }

  # ── Conditional reliability ────────────────────────────────────────────────
  } else if (type == "reliability") {

    if (!inherits(x, "csemgrm_reliability"))
      stop("For type='reliability', provide a 'csemgrm_reliability' object.")

    ttl <- title %||% "Conditional Reliability Profile (GRM)"
    df  <- x$rel_df
    p <- ggplot2::ggplot(df, ggplot2::aes(x=theta, y=crel)) +
      ggplot2::geom_line(color="#1a237e", linewidth=1.2) +
      ggplot2::geom_hline(yintercept=x$marginal_reliability,
                          linetype="dashed", color="#e65100", linewidth=0.8) +
      ggplot2::geom_hline(yintercept=0.7,
                          linetype="dotted", color="#2e7d32", linewidth=0.8) +
      ggplot2::annotate("text", x=max(df$theta)*0.7,
                        y=x$marginal_reliability + 0.02,
                        label=paste0("Marginal \u03C1 = ",
                                     x$marginal_reliability),
                        color="#e65100", size=3.2) +
      ggplot2::labs(title=ttl,
                    subtitle="Dashed: marginal reliability | Dotted: 0.70 threshold",
                    x=expression(theta~"(Latent Trait)"),
                    y="Conditional Reliability") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_minimal(base_size=13)

  # ── Item information ──────────────────────────────────────────────────────────
  } else if (type == "items") {

    if (!inherits(x, "csemgrm_fit"))
      stop("For type='items', provide a 'csemgrm_fit' object.")

    df_long <- data.frame(
      theta = rep(x$theta_grid, ncol(x$item_info)),
      info  = as.vector(x$item_info),
      item  = rep(colnames(x$item_info), each=length(x$theta_grid))
    )
    ttl <- title %||% "Item Information Functions (GRM)"
    p <- ggplot2::ggplot(df_long,
           ggplot2::aes(x=theta, y=info, color=item, group=item)) +
      ggplot2::geom_line(linewidth=0.8, alpha=0.8) +
      ggplot2::labs(title=ttl,
                    x=expression(theta~"(Latent Trait)"),
                    y="Item Information",
                    color="Item") +
      ggplot2::theme_minimal(base_size=13)

  # ── Bootstrap CSEM ────────────────────────────────────────────────────────────
  } else if (type == "bootstrap") {

    if (!inherits(x, "csemgrm_bootstrap"))
      stop("For type='bootstrap', provide a 'csemgrm_bootstrap' object.")

    ttl <- title %||%
      paste0("CSEM with ", round(x$conf_level*100), "% Bootstrap CI (GRM)")
    df <- x$csem_df
    p <- ggplot2::ggplot(df, ggplot2::aes(x=theta)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=csem_lower, ymax=csem_upper),
                           fill="#1a237e", alpha=0.15) +
      ggplot2::geom_line(ggplot2::aes(y=csem_analytical),
                         color="#1a237e", linewidth=1.2) +
      ggplot2::geom_line(ggplot2::aes(y=csem_boot_mean),
                         color="#b71c1c", linewidth=0.8, linetype="dashed") +
      ggplot2::labs(title=ttl,
                    subtitle=paste0("Solid: analytical | Dashed: bootstrap mean | ",
                                    "Shaded: ", round(x$conf_level*100), "% CI"),
                    x=expression(theta~"(Latent Trait)"),
                    y="CSEM") +
      ggplot2::theme_minimal(base_size=13)
  }

  return(p)
}

# Null coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b
