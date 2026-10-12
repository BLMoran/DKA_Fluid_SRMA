# prior_dominance() for bayesma (cmdstanr) fits

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @noRd
.bma_draws <- function(x) {
  if (inherits(x, "draws")) return(posterior::as_draws_df(x))
  if (is.data.frame(x)) return(posterior::as_draws_df(x))
  f <- x$fit %||% x$stanfit %||% x
  if (inherits(f, "CmdStanMCMC")) {
    d <- tryCatch(f$draws(format = "df"), error = function(e) NULL)
    if (!is.null(d)) return(posterior::as_draws_df(d))
    if (!is.null(x$draws)) return(posterior::as_draws_df(x$draws))
    stop("CmdStan CSV files are gone (temp dir cleared) and no `$draws` stored. ",
         "Refit, or save fits with `fit$fit$save_object()` / `output_dir =`.")
  }
  if (!is.null(x$draws)) return(posterior::as_draws_df(x$draws))
  posterior::as_draws_df(f)
}

#' @noRd
.bma_detect <- function(d, mu_par = NULL, tau_par = NULL) {
  nm <- names(d)
  list(
    mu  = mu_par  %||% intersect(c("mu", "pooled", "theta", "d"), nm)[1],
    tau = tau_par %||% intersect(c("tau", "sigma_re", "sd_re"), nm)[1]
  )
}


#' @noRd
.bma_refit <- function(x, ..., call = NULL, env = globalenv()) {
  args <- x$meta$call_args
  if (is.null(args)) stop("`meta$call_args` is missing from the bayesma output.")
  if (is.call(args)) args <- as.list(args)[-1]
  
  new <- list(...)
  args[names(new)] <- new          # replace (not merge) priors / settings
  
  if (is.null(args$data))
    stop("`meta$call_args` has no `data`; pass `data = <your data>` via `...`.")
  
  do.call(bayesma::bayesma, args, envir = env)
}


#' @noRd
.bma_summ <- function(v, prefix) {
  tibble::tibble(
    !!paste0(prefix, "_median") := stats::median(v),
    !!paste0(prefix, "_lo")     := stats::quantile(v, .025, names = FALSE),
    !!paste0(prefix, "_hi")     := stats::quantile(v, .975, names = FALSE)
  )
}


#' @export
prior_dominance <- function(fit,
                            prior_draws   = NULL,
                            mu_par        = NULL,
                            tau_par       = NULL,
                            tau_grid      = c(0.5, 0.75, 1.5),
                            tau_prior_fun = function(s) bayesma::half_normal(0, s),
                            call          = NULL,
                            env           = globalenv(),
                            ...) {
  
  post <- .bma_draws(fit)
  p    <- .bma_detect(post, mu_par, tau_par)
  if (is.na(p$mu) || is.na(p$tau))
    stop("Could not locate mu / tau draws. Pass mu_par / tau_par explicitly.")
  
  post_mu  <- post[[p$mu]]
  post_tau <- post[[p$tau]]
  
  if (is.null(prior_draws))
    prior_draws <- .bma_refit(fit, sample_prior = TRUE, ..., call = call, env = env)
  prior_draws <- .bma_draws(prior_draws)
  if (!all(c(p$mu, p$tau) %in% names(prior_draws)))
    stop("`prior_draws` must contain '", p$mu, "' and '", p$tau, "'.")
  prior_mu  <- prior_draws[[p$mu]]
  prior_tau <- prior_draws[[p$tau]]
  
  # ---- 1. contraction ----
  contraction_of <- function(prior, post) {
    pv <- stats::var(prior, na.rm = TRUE)
    if (!is.finite(pv) || pv == 0) return(NA_real_)
    1 - stats::var(post, na.rm = TRUE) / pv
  }
  interp <- function(c) dplyr::case_when(
    is.na(c)  ~ "improper prior - undefined",
    c >= 0.75 ~ "data-dominated",
    c >= 0.30 ~ "informed by both",
    TRUE      ~ "prior-dominated")
  
  cmu  <- contraction_of(prior_mu,  post_mu)
  ctau <- contraction_of(prior_tau, post_tau)
  contraction <- tibble::tibble(
    parameter      = c("mu", "tau"),
    prior_var      = c(stats::var(prior_mu), stats::var(prior_tau)),
    posterior_var  = c(stats::var(post_mu),  stats::var(post_tau)),
    contraction    = c(cmu, ctau),
    interpretation = interp(c(cmu, ctau)))
  
  # ---- 2. prior percentile at posterior median ----
  pct  <- function(prior, x) mean(prior <= x, na.rm = TRUE)
  pmu  <- pct(prior_mu,  stats::median(post_mu))
  ptau <- pct(prior_tau, stats::median(post_tau))
  percentile <- tibble::tibble(
    parameter        = c("mu", "tau"),
    posterior_median = c(stats::median(post_mu), stats::median(post_tau)),
    prior_percentile = c(pmu, ptau),
    flag             = dplyr::case_when(
      c(pmu, ptau) > 0.95 ~ "upper tail - possible prior-data conflict",
      c(pmu, ptau) < 0.05 ~ "lower tail - possible prior-data conflict",
      TRUE                ~ "interior - no conflict"))
  
  # ---- 3. tau grid ----
  tau_sensitivity <- NULL
  if (!is.null(tau_grid)) {
    tau_sensitivity <- tau_scale_sweep(fit, scales = tau_grid, mu_par = p$mu,
                                       tau_par = p$tau, tau_prior_fun = tau_prior_fun,
                                       call = call, env = env, ...) |>
      dplyr::mutate(tau_prior = sprintf("HN(0, %g)", tau_prior_sd), .before = 1) |>
      dplyr::select(-tau_prior_sd, -tau_lo, -tau_hi)
  }
  
  # ---- verdict ----
  mu_range <- if (!is.null(tau_sensitivity)) diff(range(tau_sensitivity$mu_median)) else NA_real_
  verdict <- sprintf(
    "mu contraction = %s (%s); tau contraction = %s (%s).%s",
    ifelse(is.na(cmu),  "NA", sprintf("%.2f", cmu)),  interp(cmu),
    ifelse(is.na(ctau), "NA", sprintf("%.2f", ctau)), interp(ctau),
    if (!is.na(mu_range)) sprintf(" Pooled mu shifts by %.3f across the tau-prior grid.", mu_range) else "")
  
  if (!is.null(tau_sensitivity)) {
    bad <- tau_sensitivity$tau_prior[tau_sensitivity$conflict != "compatible"]
    verdict <- paste0(verdict, if (length(bad))
      sprintf(" Prior-data conflict flagged for: %s (data want more heterogeneity than the prior allows).",
              paste(bad, collapse = ", "))
      else " No prior-data conflict across the tau-prior grid.")
  }
  
  structure(list(contraction = contraction, percentile = percentile,
                 tau_sensitivity = tau_sensitivity, verdict = verdict,
                 mu_par = p$mu, tau_par = p$tau),
            class = "prior_dominance")
}


#' @export
prior_scale_sweep <- function(fit,
                              scales       = c(0.1, 0.25, 0.5, 1, 2, 5),
                              prior_mean   = 0,
                              mu_par       = NULL,
                              mu_prior_fun = function(m, s) bayesma::normal(m, s),
                              call         = NULL,
                              env          = globalenv(),
                              ...) {
  mu_par <- .bma_detect(.bma_draws(fit), mu_par)$mu
  
  purrr::map_dfr(scales, function(s) {
    mm <- .bma_refit(fit, mu_prior = mu_prior_fun(prior_mean, s), ...,
                     call = call, env = env) |>
      .bma_draws() |>
      purrr::pluck(mu_par)
    dplyr::bind_cols(tibble::tibble(prior_sd = s),
                     .bma_summ(mm, "mu"),
                     tibble::tibble(pr_benefit = mean(mm < 0)))
  })
}


#' @export
tau_scale_sweep <- function(fit,
                            scales        = c(0.25, 0.5, 0.75, 1, 1.5, 3),
                            mu_par        = NULL,
                            tau_par       = NULL,
                            tau_prior_fun = function(s) bayesma::half_normal(0, s),
                            call          = NULL,
                            env           = globalenv(),
                            ...) {
  p <- .bma_detect(.bma_draws(fit), mu_par, tau_par)
  
  purrr::map_dfr(scales, function(s) {
    dd <- .bma_refit(fit, tau_prior = tau_prior_fun(s), ..., call = call, env = env) |>
      .bma_draws()
    mm <- dd[[p$mu]]
    tt <- dd[[p$tau]]
    tau_pctile <- 2 * stats::pnorm(stats::median(tt) / s) - 1   # half-normal CDF
    
    dplyr::bind_cols(
      tibble::tibble(tau_prior_sd = s),
      .bma_summ(mm, "mu"),
      tibble::tibble(pr_benefit = mean(mm < 0)),
      .bma_summ(tt, "tau"),
      tibble::tibble(
        tau_prior_pctile = tau_pctile,
        conflict = dplyr::if_else(tau_pctile > 0.95,
                                  "conflict: prior too tight (data want more heterogeneity)",
                                  "compatible")))
  })
}


#' @export
prior_ess <- function(fit, data = NULL, se_col = "sei", typical_se = NULL,
                      prior_draws = NULL, mu_par = NULL,
                      call = NULL, env = globalenv(), ...) {
  
  mu_par <- .bma_detect(.bma_draws(fit), mu_par)$mu
  
  if (is.null(prior_draws))
    prior_draws <- .bma_refit(fit, sample_prior = TRUE, ..., call = call, env = env)
  prior_sd <- stats::sd(.bma_draws(prior_draws)[[mu_par]], na.rm = TRUE)
  
  k <- NA_integer_
  if (is.null(typical_se)) {
    if (is.null(data) || is.null(data[[se_col]]))
      stop("Provide `data` (with `se_col`) or `typical_se`.")
    se <- data[[se_col]]
    typical_se <- sqrt(1 / mean(1 / se^2, na.rm = TRUE))
    k <- sum(is.finite(se))
  }
  
  prior_prec  <- 1 / prior_sd^2
  study_prec  <- 1 / typical_se^2
  data_prec   <- if (!is.na(k)) k * study_prec else NA_real_
  
  structure(list(prior_sd    = prior_sd,
                 typical_se  = typical_se,
                 prior_ess   = prior_prec / study_prec,
                 k           = k,
                 prior_share = if (!is.na(data_prec)) prior_prec / (prior_prec + data_prec) else NA_real_),
            class = "prior_ess")
}


#' @export
print.prior_dominance <- function(x, ...) {
  cat("Prior dominance assessment\n==========================\n\n")
  cat("Posterior contraction (1 = data-dominated, 0 = prior-dominated):\n")
  print(x$contraction); cat("\n")
  cat("Posterior median within its prior:\n")
  print(x$percentile); cat("\n")
  if (!is.null(x$tau_sensitivity)) {
    cat("Pooled effect across tau-prior grid:\n")
    print(x$tau_sensitivity); cat("\n")
  }
  cat("Verdict:\n ", x$verdict, "\n")
  invisible(x)
}

#' @export
print.prior_ess <- function(x, ...) {
  cat(sprintf(
    "Prior effective sample size (pooled effect)\n  prior SD           = %.3f\n  typical study SE   = %.3f\n  prior ESS          = %.2f study-equivalents\n",
    x$prior_sd, x$typical_se, x$prior_ess))
  if (!is.na(x$k))
    cat(sprintf("  actual studies (k) = %d\n  prior share of total precision = %.1f%%\n",
                x$k, 100 * x$prior_share))
  invisible(x)
}

#' @export
gt_prior_dominance <- function(x, decimals = 3) {
  fmt <- paste0("%.", decimals, "f")
  
  contr <- x$contraction |>
    dplyr::transmute(Section = "Contraction", Item = parameter,
                     Value = ifelse(is.na(contraction), "\u2014", sprintf(fmt, contraction)),
                     Detail = interpretation)
  
  perc <- x$percentile |>
    dplyr::transmute(Section = "Prior percentile of posterior median", Item = parameter,
                     Value = sprintf("%.1f%%", 100 * prior_percentile),
                     Detail = flag)
  
  tau <- if (!is.null(x$tau_sensitivity)) {
    x$tau_sensitivity |>
      dplyr::transmute(
        Section = "Pooled effect across tau priors", Item = tau_prior,
        Value = sprintf(paste0(fmt, " [", fmt, ", ", fmt, "]"), mu_median, mu_lo, mu_hi),
        Detail = sprintf("P(benefit) = %.2f; tau = %.2f (%.0f%% of prior; %s)",
                         pr_benefit, tau_median, 100 * tau_prior_pctile, conflict))
  }
  
  dplyr::bind_rows(contr, perc, tau) |>
    gt::gt(groupname_col = "Section", rowname_col = "Item") |>
    gt::tab_header(title = gt::md("**Prior dominance assessment**"), subtitle = x$verdict) |>
    gt::cols_label(Value = "Estimate", Detail = "") |>
    gt::tab_options(row_group.font.weight = "bold") |>
    gt::sub_missing(missing_text = "\u2014")
}

#' @export
plot_prior_scale_sweep <- function(sweep, prior_mean = 0,
                                   x_label = "mu-prior SD (log scale)",
                                   y_label = "Pooled effect",
                                   font = NULL) {
  ggplot2::ggplot(sweep, ggplot2::aes(x = prior_sd, y = mu_median)) +
    ggplot2::geom_hline(yintercept = prior_mean, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mu_lo, ymax = mu_hi), fill = "dodgerblue", alpha = 0.15) +
    ggplot2::geom_line(colour = "blue", linewidth = 0.8) +
    ggplot2::geom_point(colour = "blue", size = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::annotate("text", x = min(sweep$prior_sd), y = prior_mean, label = "prior mean",
                      hjust = 0, vjust = -0.5, size = 3, colour = "grey40", family = font) +
    ggplot2::labs(x = x_label, y = y_label,
                  subtitle = "Flat curve => estimate is data-driven; bending toward the prior mean as the prior tightens => prior influence") +
    ggplot2::theme_light() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9, family = font),
                   axis.title = ggplot2::element_text(family = font),
                   axis.text = ggplot2::element_text(family = font))
}

#' @export
plot_tau_scale_sweep <- function(sweep,
                                 mu_label = "Pooled effect (\u03bc)",
                                 tau_label = "Heterogeneity (\u03c4)",
                                 null_value = 0,
                                 font = NULL) {
  lv <- c(mu_label, tau_label)
  
  df <- purrr::map2_dfr(c("mu", "tau"), lv, \(p, lab)
                        sweep |>
                          dplyr::transmute(tau_prior_sd, panel = lab,
                                           est = .data[[paste0(p, "_median")]],
                                           lo  = .data[[paste0(p, "_lo")]],
                                           hi  = .data[[paste0(p, "_hi")]])) |>
    dplyr::mutate(panel = factor(panel, levels = lv))
  
  hlines <- tibble::tibble(panel = factor(mu_label, levels = lv), yint = null_value)
  
  ggplot2::ggplot(df, ggplot2::aes(x = tau_prior_sd, y = est)) +
    ggplot2::geom_hline(data = hlines, ggplot2::aes(yintercept = yint),
                        linetype = "dashed", colour = "grey50") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "dodgerblue", alpha = 0.15) +
    ggplot2::geom_line(colour = "blue", linewidth = 0.8) +
    ggplot2::geom_point(colour = "blue", size = 2) +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "\u03c4-prior scale (half-normal SD, log scale)", y = NULL,
                  subtitle = "Flat \u03bc across \u03c4-prior scales => pooled effect robust to the heterogeneity prior") +
    ggplot2::theme_light() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", family = font),
                   plot.subtitle = ggplot2::element_text(size = 9, family = font),
                   axis.title = ggplot2::element_text(family = font),
                   axis.text = ggplot2::element_text(family = font),
                   panel.grid.minor = ggplot2::element_blank())
}

# Function to build dataset required for each model

#' @export
build_dataset <- function(data,
                          outcome,
                          n_int_col,
                          n_ctrl_col,
                          convert_mean = FALSE,
                          mean_method  = c("wan", "luo", "qe", "bc", "shi", "yang")) {
  
  mean_methods <- c(wan = "hozo/wan/bland", luo = "luo", qe = "qe",
                    bc = "bc", shi = "shi_lognormal", yang = "yang")
  mean_method  <- mean_methods[[match.arg(mean_method)]]
  
  col <- function(suffix) paste0(outcome, suffix)
  opt_col <- function(x) if (x %in% names(data)) data[[x]] else NA_real_
  
  base <- data |>
    dplyr::mutate(
      mean.g1 = .data[[col("_mean_int")]],
      sd.g1   = .data[[col("_sd_int")]],
      med.g1  = .data[[col("_median_int")]],
      iqr.g1  = split_iqr(.data[[col("_iqr_int")]]),
      n.g1    = dplyr::coalesce(.data[[n_int_col]], n_analysed_int),
      mean.g2 = .data[[col("_mean_ctrl")]],
      sd.g2   = .data[[col("_sd_ctrl")]],
      med.g2  = .data[[col("_median_ctrl")]],
      iqr.g2  = split_iqr(.data[[col("_iqr_ctrl")]]),
      n.g2    = dplyr::coalesce(.data[[n_ctrl_col]], n_analysed_ctrl),
      rho_w   = dplyr::coalesce(opt_col(col("_rho_w")), rho_w),
      rho_b   = dplyr::coalesce(opt_col(col("_rho_b")), rho_b)) |>
    tidyr::unpack(c(iqr.g1, iqr.g2), names_sep = "_") |>
    dplyr::rename(
      q1.g1 = iqr.g1_q1,
      q3.g1 = iqr.g1_q3,
      q1.g2 = iqr.g2_q1,
      q3.g2 = iqr.g2_q3) |>
    dplyr::mutate(
      mean.g1 = if (convert_mean) mean.g1 else NA_real_,
      sd.g1   = if (convert_mean) sd.g1   else NA_real_,
      mean.g2 = if (convert_mean) mean.g2 else NA_real_,
      sd.g2   = if (convert_mean) sd.g2   else NA_real_,
      min.g1  = NA_real_,
      max.g1  = NA_real_,
      min.g2  = NA_real_,
      max.g2  = NA_real_,
      m       = (n.g1 + n.g2) / (n_clusters * n_periods),
      deff    = dplyr::coalesce(1 + (m - 1) * rho_w - m * dplyr::coalesce(rho_b, 0), 1),
      ok.g1   = !is.na(n.g1) & ((!is.na(mean.g1) & !is.na(sd.g1)) | (!is.na(med.g1) & !is.na(q1.g1) & !is.na(q3.g1))),
      ok.g2   = !is.na(n.g2) & ((!is.na(mean.g2) & !is.na(sd.g2)) | (!is.na(med.g2) & !is.na(q1.g2) & !is.na(q3.g2)))) |>
    dplyr::filter(ok.g1, ok.g2)
  
  dat_mm <- base |>
    dplyr::select(
      study_id, Author, year,
      min.g1, q1.g1, med.g1, q3.g1, max.g1, n.g1, mean.g1, sd.g1,
      min.g2, q1.g2, med.g2, q3.g2, max.g2, n.g2, mean.g2, sd.g2)
  
  df <- if (isTRUE(convert_mean)) {
    
    # Per arm: use reported mean/SD, otherwise convert with mean_method
    arm_est <- c(int = "g1", ctrl = "g2") |>
      purrr::imap(\(g, arm) {
        d <- dat_mm |>
          dplyr::select(study_id, dplyr::ends_with(paste0(".", g))) |>
          dplyr::rename_with(\(x) sub(paste0("\\.", g, "$"), ".g1", x))
        
        has_ms <- !is.na(d$mean.g1) & !is.na(d$sd.g1)
        
        reported <- tibble::tibble(
          study_id = d$study_id[has_ms],
          mean     = d$mean.g1[has_ms],
          sd       = d$sd.g1[has_ms])
        
        converted <- if (any(!has_ms)) {
          d_conv <- d[!has_ms, ]
          r <- metamedian::metamean(d_conv, mean_method = mean_method, pool_studies = FALSE)
          tibble::tibble(
            study_id = d_conv$study_id,
            mean     = r$yi,
            sd       = r$sei * sqrt(d_conv$n.g1))
        } else {
          tibble::tibble(study_id = d$study_id[0], mean = numeric(), sd = numeric())
        }
        
        dplyr::bind_rows(reported, converted) |>
          dplyr::rename(
            !!paste0("mean_", arm) := mean,
            !!paste0("sd_", arm)   := sd)
      }) |>
      purrr::reduce(dplyr::full_join, by = "study_id")
    
    base |>
      dplyr::inner_join(arm_est, by = "study_id") |>
      dplyr::filter(!is.na(mean_int), !is.na(sd_int), !is.na(mean_ctrl), !is.na(sd_ctrl)) |>
      dplyr::mutate(
        n_int     = n.g1,
        n_ctrl    = n.g2,
        age_group = as.character(age_group),
        i_summary = sprintf("%.2f (%.2f)", mean_int, sd_int),
        c_summary = sprintf("%.2f (%.2f)", mean_ctrl, sd_ctrl),
        # Design effect
        sd_int    = sd_int  * sqrt(deff),
        sd_ctrl   = sd_ctrl * sqrt(deff),
        yi        = mean_int - mean_ctrl,
        sei       = sqrt(sd_int^2 / n_int + sd_ctrl^2 / n_ctrl),
        vi        = sei^2) |>
      dplyr::select(
        study_id, Author, year, bmes_type, age_group, n_int, n_ctrl, yi, sei, vi,
        i_summary, c_summary, mean_int, sd_int, mean_ctrl, sd_ctrl, deff,
        D1:Overall)
    
  } else {
    
    res <- metamedian::metamedian(dat_mm, median_method = "qe", pool_studies = FALSE)
    
    tibble::tibble(
      study_id = dat_mm$study_id,
      yi       = res$yi,
      sei      = res$sei) |>
      dplyr::left_join(base, by = "study_id") |>
      dplyr::filter(!is.na(yi), !is.na(sei)) |>
      dplyr::mutate(
        sei         = sei * sqrt(deff),
        vi          = sei^2,
        n_int       = n.g1,
        n_ctrl      = n.g2,
        age_group   = as.character(age_group),
        i_summary   = fmt_summary(mean.g1, sd.g1, med.g1, q1.g1, q3.g1),
        c_summary   = fmt_summary(mean.g2, sd.g2, med.g2, q1.g2, q3.g2),
        median_int  = med.g1,
        q1_int      = q1.g1,
        q3_int      = q3.g1,
        iqr_int     = sprintf("%.2f\u2013%.2f", q1_int, q3_int),
        median_ctrl = med.g2,
        q1_ctrl     = q1.g2,
        q3_ctrl     = q3.g2,
        iqr_ctrl    = sprintf("%.2f\u2013%.2f", q1_ctrl, q3_ctrl)) |>
      dplyr::select(
        study_id, Author, year, bmes_type, age_group, n_int, n_ctrl, yi, sei, vi,
        i_summary, c_summary,
        median_int, q1_int, q3_int, iqr_int,
        median_ctrl, q1_ctrl, q3_ctrl, iqr_ctrl, deff,
        D1:Overall)
  }
}

#' @export
split_iqr <- function(x) {
  x |>
    stringr::str_split(",\\s*") |>
    purrr::map(\(p) readr::parse_double(stringr::str_trim(p))) |>
    purrr::map(\(p) if (length(p) == 2L) p else c(NA_real_, NA_real_)) |>
    purrr::list_transpose() |>
    rlang::set_names(c("q1", "q3")) |>
    tibble::as_tibble()
}

#' @export
design_effect <- function(design, m, rho_w, rho_b, rho_wp) {
  dplyr::case_when(
    design == "Cluster RCT" ~ 1 + (m - 1) * rho_w,
    design == "Cluster-crossover RCT" ~ 1 + (m - 1) * rho_w - m * rho_b,
    design == "Crossover RCT" ~ (1 - rho_wp) / 2,
    .default = 1)
}

#' @export
fmt_summary <- function(mean, sd, med, q1, q3) {
  dplyr::if_else(
    !is.na(mean),
    sprintf("%.2f (%.2f)", mean, sd),
    sprintf("%.2f (%.2f\u2013%.2f)", med, q1, q3))
}

#' @export
summary_row <- function(fit, outcome, measure, units, ratio = FALSE) {
  mu <- fit$draws$mu
  q <- stats::quantile(mu, c(0.5, 0.025, 0.975), names = FALSE)
  if (ratio) q <- exp(q)
  
  tibble::tibble(
    Outcome = outcome,
    Measure = measure,
    Units = units,
    `Estimate [95% CrI]` = sprintf("%.2f [%.2f, %.2f]", q[1], q[2], q[3]),
    `Pr(Benefit)` = mean(mu < 0))
}


#' @export
fmt_n <- function(e, n) sprintf("%d/%d (%.1f%%)", e, n, 100 * e / n)


#' @export
fit_binary_rr <- function(data, studyvar, event_int, event_ctrl, n_int, n_ctrl,
                          mu_prior = bayesma::normal(0, 1),
                          tau_prior = bayesma::half_normal(0, 0.75)) {
  bayesma::bayesma(
    data       = data,
    studyvar   = studyvar,
    event_int  = event_int,
    event_ctrl = event_ctrl,
    n_int      = n_int,
    n_ctrl     = n_ctrl,
    likelihood = "binomial",
    estimand   = "RR",
    model_type = "random_effect",
    stage      = "two_stage",
    mu_prior   = mu_prior,
    tau_prior  = tau_prior)
}


#' @export
binary_rr_row <- function(data, fit, label,
                          studyvar, event_int, event_ctrl, n_int, n_ctrl) {
  mu <- fit$summary$median[1]
  lo <- fit$summary$q2.5[1]
  hi <- fit$summary$q97.5[1]

  tibble::tibble(
    Outcome          = label,
    k                = dplyr::n_distinct(data[[studyvar]]),
    `Intervention`   = fmt_n(sum(data[[event_int]]),  sum(data[[n_int]])),
    `Control`        = fmt_n(sum(data[[event_ctrl]]), sum(data[[n_ctrl]])),
    `RR [95% CrI]`   = sprintf("%.2f [%.2f, %.2f]", exp(mu), exp(lo), exp(hi)),
    `Pr(Benefit)`    = mean(fit$draws$mu < 0))
}

#' @export
style_gt <- function(tbl) {
  tbl |>
    gt::tab_options(
      table.border.top.width = gt::px(2),
      table.border.top.color = "black",
      column_labels.border.top.width = gt::px(2),
      column_labels.border.top.color = "black",
      table.border.bottom.width = gt::px(2),
      table.border.bottom.color = "black",
      table_body.border.bottom.width = gt::px(2),
      table_body.border.bottom.color = "black") |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold"),
        gt::cell_borders(sides = "bottom", weight = gt::px(2), color = "black")),
      locations = gt::cells_column_labels(columns = dplyr::everything()))
}

#' @export
fmt_n <- function(e, n) sprintf("%d/%d (%.1f%%)", sum(e), sum(n), 100 * sum(e) / sum(n))

#' @export
ae_row <- function(fit, event) {
  d <- dat_ae[[event]]
  mu <- fit$draws$mu
  q <- exp(stats::quantile(mu, c(0.5, 0.025, 0.975), names = FALSE))
  
  tibble::tibble(
    `Adverse Event` = event,
    `Number of Trials` = nrow(d),
    BMES = fmt_n(d$event_int_raw, d$n_int_raw),
    Saline = fmt_n(d$event_ctrl_raw, d$n_ctrl_raw),
    `RR [95% CrI]` = sprintf("%.2f [%.2f, %.2f]", q[1], q[2], q[3]),
    `Pr(Benefit)` = mean(mu < 0))
}