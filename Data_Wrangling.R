library(tidyverse)

# Set up folders
out_dir <- here::here("Data", "Clean")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#model_dir   <- here::here("Data", "Models")
#dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

# Functions
save_rds <- function(x, name) saveRDS(x, file.path(out_dir, paste0(name, ".rds")))
load_rds <- function(name) readRDS(file.path(derived_dir, paste0(name, ".rds")))
#save_fit <- function(x, name) saveRDS(x, file.path(model_dir, paste0(name, ".rds")))

split_iqr <- function(x) {
  x |>
    stringr::str_split(",\\s*") |>
    purrr::map(\(p) readr::parse_double(stringr::str_trim(p))) |>
    purrr::map(\(p) if (length(p) == 2L) p else c(NA_real_, NA_real_)) |>
    purrr::list_transpose() |>
    rlang::set_names(c("q1", "q3")) |>
    tibble::as_tibble()
}

design_effect <- function(design, m, rho_w, rho_b, rho_wp) {
  dplyr::case_when(
    design == "Cluster RCT" ~ 1 + (m - 1) * rho_w,
    design == "Cluster-crossover RCT" ~ 1 + (m - 1) * rho_w - m * rho_b,
    design == "Crossover RCT" ~ (1 - rho_wp) / 2,
    .default = 1)
}

fmt_summary <- function(mean, sd, med, q1, q3) {
  dplyr::if_else(
    !is.na(mean),
    sprintf("%.2f (%.2f)", mean, sd),
    sprintf("%.2f (%.2f\u2013%.2f)", med, q1, q3))
}

# Function to build dataset required for each model
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

# Import Data
dat <- readxl::read_excel("Data/DKA_MA_Combined_ROB2.xlsx", sheet = "Sheet1", skip  = 1) |>
  filter(!is.na(Author)) |>
  select(-(ph_baseline_mean_int:hco3_units)) |> 
  mutate(
    study_id = str_c(Author, year, sep = " "),
    age_group = fct_relevel(as_factor(age_group), "Adult", "Paediatric"),
    bss_type = as_factor(bss_type),
    design = recode(
      trial_design,
      "rct" = "Individual-patient RCT",
      "XRCT" = "Crossover RCT",
      "cRCT" = "Cluster RCT",
      "CXRCT" = "Cluster-crossover RCT"),
    bmes_type = dplyr::case_when(
      bss_type %in% c("plasmalyte_a", "plasmalyte_ringers") ~ "Plasmalyte",
      bss_type == "ringers" ~ "Ringer's Solution",
      bss_type == "hartmanns" ~ "Hartmanns Solution",
      bss_type == "ringerfundin" ~ "Ringerfundin",
      bss_type == "sterofundin" ~ "Sterofundin",
      TRUE ~ NA_character_))

# Create datasets for each outcome (primary estimand: mean difference, Wan conversion)
dat_hosp_los <- build_dataset(dat, "hosp_los", "n_analysed_int", "n_analysed_ctrl", convert_mean = TRUE, mean_method = "wan")
dat_icu_los <- build_dataset(dat, "icu_los", "n_icu_los_int", "n_icu_los_ctrl", convert_mean = TRUE, mean_method = "wan")
dat_time_res <- build_dataset(dat, "time_res", "n_analysed_int", "n_analysed_ctrl", convert_mean = TRUE, mean_method = "wan")
dat_insulin_dur <- build_dataset(dat, "insulin_dur", "N_insulin_dur_int", "N_insulin_dur_ctrl", convert_mean = TRUE, mean_method = "wan")

# Mortality outcome (binary)
opt_col <- function(data, x) if (x %in% names(data)) data[[x]] else NA_real_

dat_mortality <- dat |>
  dplyr::mutate(
    event_int  = dplyr::coalesce(hosp_mort_int,  `28d_mort_int`),
    event_ctrl = dplyr::coalesce(hosp_mort_ctrl, `28d_mort_ctrl`),
    n_int      = n_analysed_int,
    n_ctrl     = n_analysed_ctrl,
    age_group  = as.character(age_group),
    i_summary  = sprintf("%d/%d (%.1f%%)", event_int,  n_int,  100 * event_int  / n_int),
    c_summary  = sprintf("%d/%d (%.1f%%)", event_ctrl, n_ctrl, 100 * event_ctrl / n_ctrl),
    m          = (n_int + n_ctrl) / (n_clusters * n_periods),
    deff       = dplyr::coalesce(1 + (m - 1) * rho_w - m * dplyr::coalesce(rho_b, 0), 1),
    # Effective sample size (Cochrane Handbook 23.1.4)
    event_int_eff  = round(event_int  / deff),
    n_int_eff      = round(n_int      / deff),
    event_ctrl_eff = round(event_ctrl / deff),
    n_ctrl_eff     = round(n_ctrl     / deff)) |>
  dplyr::filter(!is.na(event_int), !is.na(event_ctrl))

# Sensitivity: difference in medians (QE)
dat_hosp_los_dom <- build_dataset(dat, "hosp_los", "n_analysed_int", "n_analysed_ctrl")

# Sensitivity: ratio of medians
dat_rom <- dat_hosp_los_dom |>
  dplyr::mutate(
    log_sd_int  = (log(q3_int)  - log(q1_int))  / 1.349,
    log_sd_ctrl = (log(q3_ctrl) - log(q1_ctrl)) / 1.349,
    yi = log(median_int) - log(median_ctrl),
    vi = (pi / 2) * (log_sd_int^2 / n_int + log_sd_ctrl^2 / n_ctrl) * deff,
    sei = sqrt(vi),
    mean_int = yi,
    mean_ctrl = 0,
    sd_int = sei / sqrt(1 / n_int + 1 / n_ctrl),
    sd_ctrl = sd_int) |>
  dplyr::filter(is.finite(yi), is.finite(sei)) |>
  dplyr::select(
    study_id, Author, year, n_int, n_ctrl, yi, sei, vi, i_summary, c_summary,
    mean_int, sd_int, mean_ctrl, sd_ctrl,
    median_int, q1_int, q3_int, iqr_int,
    median_ctrl, q1_ctrl, q3_ctrl, iqr_ctrl, deff,
    D1:Overall)

# Sensitivity: mean difference by conversion method
# Yan 2024 excluded where the conversion is implausible (bimodal control arm)
dat_wan <- build_dataset(dat, "hosp_los", "n_analysed_int", "n_analysed_ctrl", convert_mean = TRUE, mean_method = "wan") |>
  rename(yi_wan = yi, sei_wan = sei, vi_wan = vi)
dat_luo <- build_dataset(dat, "hosp_los", "n_analysed_int", "n_analysed_ctrl", convert_mean = TRUE, mean_method = "luo") |>
  rename(yi_luo = yi, sei_luo = sei, vi_luo = vi)
dat_qe <- build_dataset(dat, "hosp_los", "n_analysed_int", "n_analysed_ctrl", convert_mean = TRUE, mean_method = "qe") |>
  rename(yi_qe = yi, sei_qe = sei, vi_qe = vi) |>
  filter_out(study_id %in% "Yan 2024")
dat_bc <- build_dataset(dat, "hosp_los", "n_analysed_int", "n_analysed_ctrl", convert_mean = TRUE, mean_method = "bc")  |>
  rename(yi_bc = yi, sei_bc = sei, vi_bc = vi) |>
  filter_out(study_id %in% "Yan 2024")
dat_shi <- build_dataset(dat, "hosp_los", "n_analysed_int", "n_analysed_ctrl", convert_mean = TRUE, mean_method = "shi")|>
  rename(yi_shi = yi, sei_shi = sei, vi_shi = vi) |>
  filter_out(study_id %in% "Yan 2024")
dat_yang <- build_dataset(dat, "hosp_los", "n_analysed_int", "n_analysed_ctrl", convert_mean = TRUE, mean_method = "yang")|>
  rename(yi_yang = yi, sei_yang = sei, vi_yang = vi)

dat_mean <- dat_wan |>
  select(study_id:vi_wan) |>
  left_join(dat_luo |> select(study_id, yi_luo, sei_luo, vi_luo), by = "study_id") |>
  left_join(dat_qe |> select(study_id, yi_qe, sei_qe, vi_qe), by = "study_id") |>
  left_join(dat_bc |> select(study_id, yi_bc, sei_bc, vi_bc), by = "study_id") |>
  left_join(dat_shi |> select(study_id, yi_shi, sei_shi, vi_shi), by = "study_id") |>
  left_join(dat_yang |> select(study_id, yi_yang, sei_yang, vi_yang), by = "study_id")

# Trial Counts
counts <- list(
  n_trials = nrow(dat),
  n_participants = sum(dat$n_analysed_int + dat$n_analysed_ctrl, na.rm = TRUE),
  n_primary = nrow(dat_hosp_los),
  n_primary_pts = sum(dat_hosp_los$n_int + dat_hosp_los$n_ctrl, na.rm = TRUE))

# Adverse Events
adverse_events <- c(
  "Mortality" = "Mortality",
  "Hypokalaemia" = "Hypokalaemia",
  "Hypoglycaemia" = "Hypoglycaemia",
  "Mechanical ventilation" = "Mechanical ventilation",
  "AKI" = "AKI")

# Binary-outcome ICCs (none reported; BEST-DKA hospital LOS ICCs assumed)
icc_ae <- tibble::tribble(
  ~study_id, ~rho_w_ae, ~rho_b_ae,
  "Ramanan 2026", 0.0108, 0)

design_ae <- dat |>
  dplyr::distinct(study_id, design, n_clusters, n_periods, rho_w, rho_b)

dat_ae_raw <- readxl::read_excel("Data/Adverse.xlsx") |>
  dplyr::mutate(
    study_id = paste(Study, Year),
    Event = dplyr::if_else(grepl("^Death", Outcome), "Mortality", Outcome)) |>
  dplyr::filter(Event %in% adverse_events)

build_ae_dataset <- function(data, event) {
  data |>
    dplyr::filter(Event == event) |>
    dplyr::left_join(design_ae, by = "study_id") |>
    dplyr::left_join(icc_ae, by = "study_id") |>
    dplyr::mutate(
      rho_w = dplyr::coalesce(rho_w_ae, rho_w),
      rho_b = dplyr::coalesce(rho_b_ae, rho_b, 0),
      m = (Balanced_N + Saline_N) / (n_clusters * n_periods),
      deff = dplyr::coalesce(design_effect(design, m, rho_w, rho_b, NA_real_), 1),
      i_summary = sprintf("%d/%d (%.1f%%)", Balanced_events, Balanced_N, 100 * Balanced_events / Balanced_N),
      c_summary = sprintf("%d/%d (%.1f%%)", Saline_events, Saline_N, 100 * Saline_events / Saline_N),
      # DEFF-adjusted effective counts (Cochrane Handbook 23.1.4)
      event_int = round(Balanced_events / deff),
      n_int = round(Balanced_N / deff),
      event_ctrl = round(Saline_events / deff),
      n_ctrl = round(Saline_N / deff)) |>
    dplyr::filter(!(event_int == 0 & event_ctrl == 0)) |>
    dplyr::select(
      study_id, Author = Study, year = Year, Event, design, deff, i_summary, c_summary,
      event_int, n_int, event_ctrl, n_ctrl, event_int_raw = Balanced_events, n_int_raw = Balanced_N,
      event_ctrl_raw = Saline_events, n_ctrl_raw = Saline_N)
}

dat_ae <- purrr::map(adverse_events, \(ev) build_ae_dataset(dat_ae_raw, ev))

save_rds(dat, "dat")
save_rds(dat_hosp_los, "dat_hosp_los")
save_rds(dat_icu_los, "dat_icu_los")
save_rds(dat_time_res, "dat_time_res")
save_rds(dat_insulin_dur, "dat_insulin_dur")
save_rds(dat_mortality, "dat_mortality")
save_rds(dat_hosp_los_dom, "dat_hosp_los_dom")
save_rds(dat_rom, "dat_rom")
save_rds(dat_mean, "dat_mean")
save_rds(counts, "counts")
save_rds(dat_ae, "dat_ae")

