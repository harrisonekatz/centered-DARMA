## =========================================================
##  Weekly Finance Application: B-DARMA (ALR) — with Rolling Origins
##  Data: Fed H.8 (weekly, SA) — bank asset composition
##  Components: Cash, Securities, Loans, Other (= Total - Cash - Securities - Loans)
##  Models: Raw-MA vs Centered-MA (ALR space), Dirichlet likelihood
##  Outputs: EDA, ADF validation, diagnostics, fixed holdout metrics,
##           ROLLING 1-step metrics over a 2–3 year window (avg ELPD diff)
## =========================================================

## -------------------- USER CONFIG --------------------
CONFIG <- list(
  sim_name      = "h8_weekly_bank_assets_composition",
  base_dir      = "results",
  seed          = 20251015,
  
  # Trim weekly history to keep HMC stable (set to Inf to keep full)
  keep_years    = 10,
  
  # Fixed holdout (kept for comparability)
  test_weeks    = 104,
  
  # Stan MCMC (safer defaults for this model)
  chains        = 4,
  iter          = 2000,
  warmup        = 1000,
  adapt_delta   = 0.9,
  max_treedepth = 12,
  
  # Auto-refit policy (if diagnostics fail)
  auto_refit    = TRUE,
  rhat_max      = 1.01,
  ess_bulk_min  = 400,
  max_refits    = 1,
  refit_factor  = 2,
  
  # Floors for numerics (Dirichlet / simplex)
  eps_prob      = 1e-10,
  eps_shape     = 1e-10,
  
  # Weekly H.8 FRED series candidates (first that works is used)
  fred_ids = list(
    total = c("TLAACBW027SBOG"),
    cash  = c("CASACBW027SBOG"),
    securities = c("SBCACBW027SBOG", "SBCACBW027NBOG"),  # NSA fallback if SA unavailable
    loans = c("TOTLL", "LLBACBW027SBOG")
  ),
  
  # Forecast-draws subsample for ELPD / coverage (speed)
  S_sub = 400,
  
  # --------------- Rolling 1-step evaluation settings ---------------
  rolling = list(
    do            = TRUE,     # turn ON/OFF
    weeks         = 104,      # evaluate last ~2 years (1-step ahead)
    step          = 1,        # 1=weekly origins, 4=monthly origins
    min_train     = 104,      # minimum training weeks per origin (>= 2 years recommended)
    
    # Lighter sampler for the many re-fits (adjust to taste)
    chains        = 2,
    iter          = 1200,
    warmup        = 600,
    adapt_delta   = 0.95,
    max_treedepth = 12,
    S_sub         = 200       # fewer draws for rolling ELPD/Coverage speed
  )
)

## ---- Safety: single-threaded math on shared servers ----
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

## -------------------- Packages -----------------------
suppressPackageStartupMessages({
  need <- c(
    "rstan","dplyr","tibble","purrr","matrixStats","readr","fs","stringr","tidyr",
    "lubridate","ggplot2","scales","urca","stats"
  )
  have <- rownames(installed.packages())
  miss <- setdiff(need, have)
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
  lapply(need, library, character.only = TRUE)
})
rstan_options(auto_write = TRUE)
options(mc.cores = 1L, warn = 1)

## -------------------- Results directory ----------------
set.seed(CONFIG$seed)
fs::dir_create(CONFIG$base_dir)
root_outdir <- file.path(CONFIG$base_dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", CONFIG$sim_name))
fs::dir_create(root_outdir); fs::dir_create(file.path(root_outdir, "figs")); fs::dir_create(file.path(root_outdir, "tables"))
writeLines(capture.output(str(CONFIG)), file.path(root_outdir, "config.txt"))

## =========================================================
## 1) DATA: Pull weekly H.8 series from FRED, build composition
## =========================================================

fred_csv <- function(series_id) {
  url <- paste0("https://fred.stlouisfed.org/graph/fredgraph.csv?id=", series_id)
  df <- tryCatch(
    readr::read_csv(url, show_col_types = FALSE, progress = FALSE),
    error = function(e) stop("Failed to download/parse FRED CSV for ", series_id, "\nURL: ", url, "\n", e$message)
  )
  nms <- names(df)
  date_col <- intersect(nms, c("DATE","date","observation_date"))
  if (length(date_col) == 0L) stop("No date column found for series ", series_id, ". Columns: ", paste(nms, collapse=", "))
  val_col <- setdiff(nms, date_col)
  if (length(val_col) != 1L) stop("Ambiguous value column for ", series_id, ". Columns: ", paste(nms, collapse=", "))
  df |>
    dplyr::transmute(date = as.Date(.data[[date_col]]),
                     value = suppressWarnings(as.numeric(.data[[val_col]]))) |>
    dplyr::filter(!is.na(date))
}

try_first_working <- function(id_vec, label) {
  for (sid in id_vec) {
    cat("  Trying", label, "series:", sid, "\n")
    ok <- try(fred_csv(sid), silent = TRUE)
    if (!inherits(ok, "try-error")) {
      cat("    -> using", sid, "\n")
      ok$series <- label; ok$id <- sid; return(ok)
    }
  }
  stop("None of the candidate FRED IDs worked for ", label, ". Please replace with a valid weekly H.8 ID.")
}

message("Downloading weekly H.8 series from FRED ...")
df_total <- try_first_working(CONFIG$fred_ids$total,      "total")
df_cash  <- try_first_working(CONFIG$fred_ids$cash,       "cash")
df_secr  <- try_first_working(CONFIG$fred_ids$securities, "securities")
df_loans <- try_first_working(CONFIG$fred_ids$loans,      "loans")

# Merge to a single weekly panel (common dates only), then optionally trim to last N years
df_w <- df_total |>
  dplyr::select(date, total = value) |>
  dplyr::inner_join(df_cash  |> dplyr::select(date, cash = value),       by = "date") |>
  dplyr::inner_join(df_secr  |> dplyr::select(date, securities = value), by = "date") |>
  dplyr::inner_join(df_loans |> dplyr::select(date, loans = value),      by = "date") |>
  dplyr::arrange(date)

if (is.finite(CONFIG$keep_years)) {
  cutoff <- max(df_w$date, na.rm = TRUE) - lubridate::years(CONFIG$keep_years)
  df_w <- dplyr::filter(df_w, date >= cutoff)
  message("Keeping last ", CONFIG$keep_years, " years: T = ", nrow(df_w))
}

# 'Other' bucket, sanity check
df_w <- df_w |>
  dplyr::mutate(other = total - cash - securities - loans)

neg_other <- sum(df_w$other < 0, na.rm = TRUE)
if (neg_other > 0) {
  frac <- round(100*neg_other/nrow(df_w), 2)
  warning(neg_other, " rows (", frac, "%) have negative 'other' (ID mismatch or SA/NSA mix).")
  if (neg_other > 0.05 * nrow(df_w)) {
    stop("Too many negative 'other' rows (", frac, "%). Check the loans/securities IDs or SA/NSA consistency.")
  }
}

# Build composition (shares)
.EPSP <- CONFIG$eps_prob
.enforce_simplex_eps <- function(Y, eps = .EPSP) {
  Y <- as.matrix(Y); Y[!is.finite(Y)] <- eps; Y[Y < eps] <- eps
  rs <- rowSums(Y); rs[rs <= 0] <- 1
  sweep(Y, 1, rs, "/")
}

Y <- df_w |>
  dplyr::transmute(
    date,
    cash_share       = pmax(cash/total, 0),
    securities_share = pmax(securities/total, 0),
    loans_share      = pmax(loans/total, 0),
    other_share      = pmax(1 - cash_share - securities_share - loans_share, 0)
  )

Y_mat <- .enforce_simplex_eps(as.matrix(Y[, -1]))
stopifnot(all(abs(rowSums(Y_mat) - 1) < 1e-10))

readr::write_csv(
  dplyr::bind_cols(Y["date"], tibble::as_tibble(Y_mat, .name_repair = ~c("cash","securities","loans","other"))),
  file.path(root_outdir, "h8_weekly_composition.csv")
)

## =========================================================
## 2) EDA: Plots & basic diagnostics
## =========================================================

eda_long <- dplyr::bind_cols(Y["date"], tibble::as_tibble(Y_mat, .name_repair = ~c("cash","securities","loans","other"))) |>
  tidyr::pivot_longer(-date, names_to = "component", values_to = "share")

p_stack <- ggplot2::ggplot(eda_long, aes(x = date, y = share, fill = component)) +
  ggplot2::geom_area(alpha = 0.85, color = "white", size = 0.1) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(title = "H.8 Bank Assets — Composition (shares of total, weekly SA)",
                x = NULL, y = "Share of total assets",
                caption = "Sources: FRB H.8 via FRED") +
  ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(root_outdir, "figs", "h8_weekly_stack_shares.png"), p_stack, width = 10, height = 5, dpi = 200)

summ_tbl <- eda_long |>
  dplyr::group_by(component) |>
  dplyr::summarise(
    mean   = mean(share),
    sd     = sd(share),
    p01    = quantile(share, 0.01),
    p99    = quantile(share, 0.99),
    min    = min(share),
    max    = max(share),
    n      = dplyr::n(), .groups = "drop"
  )
readr::write_csv(summ_tbl, file.path(root_outdir, "tables", "summary_shares_weekly.csv"))

## ALR transform & ADF: choose j_star = loans (largest/persistent)
alr <- function(y, j_star) {
  J <- length(y)
  idx <- setdiff(seq_len(J), j_star)
  ly <- log(pmax(y, .EPSP))
  ly[idx] - ly[j_star]
}
J <- 4L; j_star <- 3L
eta <- t(apply(Y_mat, 1L, alr, j_star = j_star))
colnames(eta) <- c("eta_cash","eta_securities","eta_other")

# ADF with drift, lag=6
adf_res <- lapply(1:ncol(eta), function(k) urca::ur.df(eta[,k], type = "drift", lags = 6))
adf_tbl <- tibble::tibble(
  component = colnames(eta),
  teststat  = sapply(adf_res, function(u) u@teststat["tau2"]),
  crit_5pct = sapply(adf_res, function(u) u@cval["tau2","5pct"]),
  reject_5pct = teststat < crit_5pct
)
readr::write_csv(adf_tbl, file.path(root_outdir, "tables", "adf_alr_weekly.csv"))

## =========================================================
## 3) z_t (precision regressors): lagged, smoothed ALR volatility
## =========================================================

# Realized ALR "volatility": sqrt(mean((eta_t - eta_{t-1})^2))
deta <- rbind(eta[1, , drop=FALSE], eta[-1, , drop=FALSE] - eta[-nrow(eta), , drop=FALSE])
rv <- sqrt(rowMeans(deta^2))

# 4-week trailing mean; then lag by 1 (no look-ahead)
w <- rep(1/4, 4)
rv_smooth <- as.numeric(stats::filter(rv, filter = w, sides = 1))
if (anyNA(rv_smooth)) rv_smooth[is.na(rv_smooth)] <- rv_smooth[which(!is.na(rv_smooth))[1]]
z_vol <- c(rv_smooth[1], rv_smooth[-length(rv_smooth)])  # lag by 1

## =========================================================
## 4) Stan models (Raw-MA vs Centered-MA)
## =========================================================

stan_code_centered <- "
functions {
  vector alr_simplex(vector y, int j_star) {
    int J = num_elements(y); int K = J - 1; vector[K] out; int k = 1;
    for (j in 1:J) if (j != j_star) { out[k] = log(y[j]) - log(y[j_star]); k += 1; }
    return out;
  }
  vector alr_inv(vector eta, int J, int j_star) {
    int K = J - 1; vector[J] mu; vector[K] e = exp(eta); real denom = 1 + sum(e);
    mu[j_star] = 1 / denom;
    { int k = 1; for (j in 1:J) if (j != j_star) { mu[j] = e[k] / denom; k += 1; } }
    return mu;
  }
  vector E_alr_dirichlet(vector mu, real phi, int j_star) {
    int J = num_elements(mu); vector[J] a = phi * mu; vector[J] dg = digamma(a); real dg_ref = dg[j_star];
    int K = J - 1; vector[K] out; int k = 1;
    for (j in 1:J) if (j != j_star) { out[k] = dg[j] - dg_ref; k += 1; }
    return out;
  }
}
data {
  int<lower=1> T; int<lower=2> J; int<lower=1,upper=J> j_star; int<lower=1> K;
  int<lower=0> P; int<lower=0> Q; int<lower=1> R; int<lower=1> R_phi;
  matrix[T, J] y; matrix[T, R] x; matrix[T, R_phi] z;
}
parameters {
  array[P] matrix[K, K] A; array[Q] matrix[K, K] B; matrix[K, R] beta; vector[R_phi] gamma;
}
transformed parameters {
  matrix[T, K] eta; matrix[T, J] mu; vector[T] phi; matrix[T, K] e_centered;
  for (t in 1:T) {
    vector[K] rhs = beta * x[t]';
    for (p in 1:P) if (t > p) rhs += A[p] * ( alr_simplex(y[t-p]', j_star) - beta * x[t-p]' );
    for (q in 1:Q) if (t > q) rhs += B[q] * ( e_centered[t-q]' );
    eta[t] = rhs';
    mu[t] = (alr_inv(rhs, J, j_star))';
    phi[t] = exp(z[t] * gamma);
    {
      vector[K] ec = alr_simplex(y[t]', j_star) - E_alr_dirichlet(mu[t]', phi[t], j_star);
      e_centered[t] = ec';
    }
  }
}
model {
  for (p in 1:P) to_vector(A[p]) ~ normal(0, 0.5);
  for (q in 1:Q) to_vector(B[q]) ~ normal(0, 0.5);
  to_vector(beta) ~ normal(0, 1); gamma ~ normal(0, 1);
  for (t in 1:T) target += dirichlet_lpdf( to_vector(y[t]) | phi[t] * to_vector(mu[t]) );
}
"

stan_code_raw <- "
functions {
  vector alr_simplex(vector y, int j_star) {
    int J = num_elements(y); int K = J - 1; vector[K] out; int k = 1;
    for (j in 1:J) if (j != j_star) { out[k] = log(y[j]) - log(y[j_star]); k += 1; }
    return out;
  }
  vector alr_inv(vector eta, int J, int j_star) {
    int K = J - 1; vector[J] mu; vector[K] e = exp(eta); real denom = 1 + sum(e);
    mu[j_star] = 1 / denom;
    { int k = 1; for (j in 1:J) if (j != j_star) { mu[j] = e[k] / denom; k += 1; } }
    return mu;
  }
}
data {
  int<lower=1> T; int<lower=2> J; int<lower=1,upper=J> j_star; int<lower=1> K;
  int<lower=0> P; int<lower=0> Q; int<lower=1> R; int<lower=1> R_phi;
  matrix[T, J] y; matrix[T, R] x; matrix[T, R_phi] z;
}
parameters {
  array[P] matrix[K, K] A; array[Q] matrix[K, K] B; matrix[K, R] beta; vector[R_phi] gamma;
}
transformed parameters {
  matrix[T, K] eta; matrix[T, J] mu; vector[T] phi;
  for (t in 1:T) {
    vector[K] rhs = beta * x[t]';
    for (p in 1:P) if (t > p) rhs += A[p] * ( alr_simplex(y[t-p]', j_star) - beta * x[t-p]' );
    for (q in 1:Q) if (t > q) rhs += B[q] * ( alr_simplex(y[t-q]', j_star) - eta[t-q]' );
    eta[t] = rhs';
    mu[t] = (alr_inv(rhs, J, j_star))';
    phi[t] = exp(z[t] * gamma);
  }
}
model {
  for (p in 1:P) to_vector(A[p]) ~ normal(0, 0.5);
  for (q in 1:Q) to_vector(B[q]) ~ normal(0, 0.5);
  to_vector(beta) ~ normal(0, 1); gamma ~ normal(0, 1);
  for (t in 1:T) target += dirichlet_lpdf( to_vector(y[t]) | phi[t] * to_vector(mu[t]) );
}
"

message("Compiling Stan models ...")
sm_centered <- rstan::stan_model(model_code = stan_code_centered)
sm_raw      <- rstan::stan_model(model_code = stan_code_raw)

## =========================================================
## 5) Helper functions (numerics, forecasting, metrics, diagnostics)
## =========================================================

.EPSS <- CONFIG$eps_shape

alr_inv_soft <- function(eta, J, j_star) {
  z <- numeric(J); k <- 1L
  for (j in 1:J) { if (j == j_star) z[j] <- 0 else { z[j] <- eta[k]; k <- k + 1L } }
  cst <- max(z); ez <- exp(pmin(z - cst, 700))
  mu <- ez / sum(ez); mu <- pmax(mu, .EPSP); mu / sum(mu)
}
E_alr_dirichlet_R <- function(mu, phi, j_star) {
  a <- pmax(phi * pmax(mu, .EPSP), .EPSS)
  dg <- digamma(a)
  idx <- setdiff(seq_along(mu), j_star)
  dg[idx] - dg[j_star]
}
validate_for_stan <- function(y) {
  if (any(!is.finite(y))) stop("non-finite in y")
  if (any(is.na(y))) stop("NA in y")
  if (any(y <= 0)) stop("Dirichlet requires strictly positive entries")
  if (any(abs(rowSums(y) - 1) > 1e-10)) stop("rows must sum to 1")
  invisible(TRUE)
}

hmc_diagnostics <- function(fit, rhat_max, ess_bulk_min, target_treedepth) {
  s <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  div <- sum(sapply(s, function(ch) sum(ch[, "divergent__"])))
  td_hits <- sum(sapply(s, function(ch) sum(ch[, "treedepth__"] >= target_treedepth)))
  summ <- summary(fit)$summary
  sel <- grepl("^(A\\[|B\\[|beta\\[|gamma\\[)", rownames(summ))
  sub <- summ[sel, , drop = FALSE]
  rhat_bad <- sum(sub[, "Rhat"] > rhat_max, na.rm = TRUE)
  ess_bad  <- sum(sub[, "n_eff"] < ess_bulk_min, na.rm = TRUE)
  list(divergences = div, max_treedepth_hits = td_hits,
       n_rhat_gt = rhat_bad, n_low_ess_bulk = ess_bad)
}

fit_bayesian_checked <- function(sm, y, x, z, J, j_star, P, Q,
                                 iter, warmup, chains, adapt_delta, max_treedepth,
                                 seed, label) {
  validate_for_stan(y)
  data_list <- list(T=nrow(y), J=J, j_star=j_star, K=J-1L, P=P, Q=Q,
                    R=ncol(x), R_phi=ncol(z), y=y, x=x, z=z)
  ctrl <- list(adapt_delta=adapt_delta, max_treedepth=max_treedepth)
  attempt <- 0
  repeat {
    attempt <- attempt + 1
    fit <- sampling(sm, data=data_list, iter=iter, warmup=warmup, chains=chains,
                    seed=seed + attempt, control=ctrl, refresh=200, init=0)
    dg <- hmc_diagnostics(fit, CONFIG$rhat_max, CONFIG$ess_bulk_min, max_treedepth)
    bad <- (dg$divergences > 0) || (dg$n_rhat_gt > 0) || (dg$n_low_ess_bulk > 0)
    if (!CONFIG$auto_refit || !bad || attempt > CONFIG$max_refits) {
      return(list(fit=fit, diag=dg, iter=iter, warmup=warmup, chains=chains,
                  adapt_delta=adapt_delta, max_treedepth=max_treedepth, attempts=attempt))
    }
    message(sprintf("[%s] diagnostics flagged; refitting with more iterations (attempt %d)", label, attempt))
    iter <- as.integer(iter * CONFIG$refit_factor)
    warmup <- as.integer(warmup * CONFIG$refit_factor)
    adapt_delta <- min(0.999, adapt_delta + 0.01)
    gc()
  }
}

metrics_overall <- function(y_test, mu_hat) {
  y_test <- as.matrix(y_test); mu_hat <- as.matrix(mu_hat)
  sq <- (y_test - mu_hat)^2; ab <- abs(y_test - mu_hat)
  tibble::tibble(
    component = c(paste0("y_", seq_len(ncol(y_test))), "Total"),
    RMSE = c(sqrt(colMeans(sq, na.rm = TRUE)), sqrt(mean(sq, na.rm = TRUE))),
    MAE  = c(colMeans(ab, na.rm = TRUE),       mean(ab, na.rm = TRUE))
  )
}
dirichlet_logpdf <- function(y, a, eps_y = CONFIG$eps_prob, eps_a = CONFIG$eps_shape) {
  a <- pmax(as.numeric(a), eps_a); yy <- pmax(as.numeric(y), eps_y)
  lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * log(yy))
}
elpd_one_step <- function(y_test, mu_draws, phi_draws) {
  Tt <- nrow(y_test); S <- ncol(phi_draws); lpd <- numeric(Tt)
  for (t in 1:Tt) {
    lp <- numeric(S)
    for (s in 1:S) lp[s] <- dirichlet_logpdf(y_test[t, ], pmax(phi_draws[t, s] * mu_draws[t, , s], CONFIG$eps_shape))
    lpd[t] <- matrixStats::logSumExp(lp) - log(S)
  }
  sum(lpd)
}
predictive_coverage_y <- function(y_test, mu_draws, phi_draws, level=0.95, seed=1) {
  set.seed(seed); Tt <- nrow(y_test); J <- ncol(y_test); S <- ncol(phi_draws)
  lo <- (1-level)/2; hi <- 1 - lo; inside <- 0L; total <- Tt*J
  for (t in 1:Tt) {
    yrep <- matrix(NA_real_, S, J)
    for (s in 1:S) {
      a <- pmax(phi_draws[t, s] * mu_draws[t, , s], CONFIG$eps_shape)
      g <- rgamma(J, shape = a, rate = 1); g[!is.finite(g) | g <= 0] <- CONFIG$eps_prob
      yr <- g / sum(g); yr <- pmax(yr, CONFIG$eps_prob); yr <- yr / sum(yr)
      yrep[s, ] <- yr
    }
    ql <- apply(yrep, 2, stats::quantile, probs=lo); qu <- apply(yrep, 2, stats::quantile, probs=hi)
    yrow <- pmax(y_test[t, ], CONFIG$eps_prob); yrow <- yrow / sum(yrow)
    inside <- inside + sum(yrow >= ql & yrow <= qu)
  }
  inside / total
}

forecast_1step_mean <- function(draws, y_all, x, z, J, j_star, P, Q,
                                model_type=c("raw","centered"), T_train, T_test) {
  model_type <- match.arg(model_type)
  T_total <- nrow(x); idx <- seq.int(T_train+1L, T_total); K <- J-1L
  S <- dim(draws$beta)[1]; mu_hat <- array(NA_real_, c(length(idx), J, S))
  for (s in 1:S) {
    beta_s <- draws$beta[s,,]; if (is.null(dim(beta_s))) beta_s <- matrix(beta_s, nrow=K, ncol=1)
    gamma_s <- as.numeric(draws$gamma[s,]); A_s <- array(0, c(K,K,P)); B_s <- array(0, c(K,K,Q))
    if (P>0) for (p in 1:P) A_s[,,p] <- draws$A[s,p,,]
    if (Q>0) for (q in 1:Q) B_s[,,q] <- draws$B[s,q,,]
    eta_s <- matrix(0, T_total, K)
    for (t in 1:T_total) {
      rhs <- drop(beta_s %*% x[t,])
      if (P>0) for (p in 1:P) if (t>p) rhs <- rhs + A_s[,,p] %*% (alr(y_all[t-p,], j_star) - drop(beta_s %*% x[t-p, ]))
      if (Q>0) for (q in 1:Q) if (t>q) {
        if (model_type=="raw") rhs <- rhs + B_s[,,q] %*% (alr(y_all[t-q,], j_star) - eta_s[t-q,])
        else {
          mu_tq  <- alr_inv_soft(eta_s[t-q,], J, j_star); phi_tq <- exp(drop(z[t-q,] %*% gamma_s))
          ecent  <- alr(y_all[t-q,], j_star) - E_alr_dirichlet_R(mu_tq, phi_tq, j_star)
          rhs <- rhs + B_s[,,q] %*% ecent
        }
      }
      eta_s[t,] <- rhs
      if (t %in% idx) mu_hat[which(idx==t),,s] <- alr_inv_soft(eta_s[t,], J, j_star)
    }
  }
  apply(mu_hat, c(1,2), mean, na.rm=TRUE)
}

forecast_1step_draws <- function(draws, y_all, x, z, J, j_star, P, Q,
                                 model_type=c("raw","centered"), T_train, T_test,
                                 S_sub=NULL, seed=1) {
  set.seed(seed); model_type <- match.arg(model_type)
  T_total <- nrow(x); idx <- seq.int(T_train+1L, T_total); K <- J-1L
  S <- dim(draws$beta)[1]; sel <- if (is.null(S_sub)||S_sub>=S) 1:S else sort(sample(1:S, S_sub))
  S2 <- length(sel); mu_draws <- array(NA_real_, c(length(idx), J, S2)); phi_draws <- matrix(NA_real_, length(idx), S2)
  for (ii in seq_along(sel)) {
    s <- sel[ii]
    beta_s <- draws$beta[s,,]; if (is.null(dim(beta_s))) beta_s <- matrix(beta_s, nrow=K, ncol=1)
    gamma_s <- as.numeric(draws$gamma[s,]); A_s <- array(0, c(K,K,P)); B_s <- array(0, c(K,K,Q))
    if (P>0) for (p in 1:P) A_s[,,p] <- draws$A[s,p,,]
    if (Q>0) for (q in 1:Q) B_s[,,q] <- draws$B[s,q,,]
    eta_s <- matrix(0, T_total, K)
    for (t in 1:T_total) {
      rhs <- drop(beta_s %*% x[t,])
      if (P>0) for (p in 1:P) if (t>p) rhs <- rhs + A_s[,,p] %*% (alr(y_all[t-p,], j_star) - drop(beta_s %*% x[t-p, ]))
      if (Q>0) for (q in 1:Q) if (t>q) {
        if (model_type=="raw") rhs <- rhs + B_s[,,q] %*% (alr(y_all[t-q,], j_star) - eta_s[t-q,])
        else {
          mu_tq  <- alr_inv_soft(eta_s[t-q,], J, j_star); phi_tq <- exp(drop(z[t-q,] %*% gamma_s))
          ecent  <- alr(y_all[t-q,], j_star) - E_alr_dirichlet_R(mu_tq, phi_tq, j_star)
          rhs <- rhs + B_s[,,q] %*% ecent
        }
      }
      eta_s[t,] <- rhs
      if (t %in% idx) { k <- which(idx==t); mu_draws[k,,ii] <- alr_inv_soft(eta_s[t,], J, j_star); phi_draws[k,ii] <- exp(drop(z[t,] %*% gamma_s)) }
    }
  }
  list(mu_draws=mu_draws, phi_draws=phi_draws)
}

## =========================================================
## 6) Fixed Train/Test split, build X and Z (standardized on TRAIN)
## =========================================================

T_total <- nrow(Y_mat)
T_test  <- min(CONFIG$test_weeks, floor(T_total * 0.25))
T_train <- T_total - T_test
stopifnot(T_train >= 52)

y_tr <- Y_mat[1:T_train, , drop=FALSE]
y_te <- Y_mat[(T_train+1):T_total, , drop=FALSE]

# X: intercept only
x <- matrix(1, nrow = T_total, ncol = 1)

# Z: intercept + standardized lagged ALR-vol (no look-ahead; standardize on TRAIN only)
z_mean <- mean(z_vol[1:T_train], na.rm = TRUE)
z_sd   <- sd(z_vol[1:T_train], na.rm = TRUE); if (!is.finite(z_sd) || z_sd <= 0) z_sd <- 1
z_std  <- as.numeric((z_vol - z_mean) / z_sd)
z <- cbind(1, z_std)  # R_phi = 2

P <- 1L; Q <- 1L; R <- ncol(x); R_phi <- ncol(z)
J <- ncol(Y_mat); K <- J - 1L

## =========================================================
## 7) Fit both models on fixed training data (baseline comparison)
## =========================================================

message("Fitting Centered-MA (fixed holdout) ...")
fit_c <- fit_bayesian_checked(
  sm_centered, y_tr, x[1:T_train,,drop=FALSE], z[1:T_train,,drop=FALSE],
  J=J, j_star=j_star, P=P, Q=Q,
  iter=CONFIG$iter, warmup=CONFIG$warmup, chains=CONFIG$chains,
  adapt_delta=CONFIG$adapt_delta, max_treedepth=CONFIG$max_treedepth,
  seed=CONFIG$seed, label="Centered MA"
)

message("Fitting Raw-MA (fixed holdout) ...")
fit_r <- fit_bayesian_checked(
  sm_raw, y_tr, x[1:T_train,,drop=FALSE], z[1:T_train,,drop=FALSE],
  J=J, j_star=j_star, P=P, Q=Q,
  iter=CONFIG$iter, warmup=CONFIG$warmup, chains=CONFIG$chains,
  adapt_delta=CONFIG$adapt_delta, max_treedepth=CONFIG$max_treedepth,
  seed=CONFIG$seed + 1000, label="Raw MA"
)

diag_tbl <- tibble::tibble(
  model = c("Centered MA","Raw MA"),
  divergences = c(fit_c$diag$divergences, fit_r$diag$divergences),
  max_treedepth_hits = c(fit_c$diag$max_treedepth_hits, fit_r$diag$max_treedepth_hits),
  n_rhat_gt = c(fit_c$diag$n_rhat_gt, fit_r$diag$n_rhat_gt),
  n_low_ess_bulk = c(fit_c$diag$n_low_ess_bulk, fit_r$diag$n_low_ess_bulk),
  attempts = c(fit_c$attempts, fit_r$attempts),
  iter = c(fit_c$iter, fit_r$iter),
  warmup = c(fit_c$warmup, fit_r$warmup)
)
readr::write_csv(diag_tbl, file.path(root_outdir, "tables", "hmc_diagnostics_weekly.csv"))

## =========================================================
## 8) Forecast 1-step ahead on fixed holdout & evaluate
## =========================================================

draws_c <- rstan::extract(fit_c$fit, permuted=TRUE)
draws_r <- rstan::extract(fit_r$fit, permuted=TRUE)

mu_hat_c <- forecast_1step_mean(draws_c, Y_mat, x, z, J, j_star, P, Q, "centered", T_train, T_test)
mu_hat_r <- forecast_1step_mean(draws_r, Y_mat, x, z, J, j_star, P, Q, "raw",      T_train, T_test)

m_c <- metrics_overall(y_te, mu_hat_c)
m_r <- metrics_overall(y_te, mu_hat_r)

dr_c <- forecast_1step_draws(draws_c, Y_mat, x, z, J, j_star, P, Q, "centered", T_train, T_test, S_sub = CONFIG$S_sub, seed=CONFIG$seed+3)
dr_r <- forecast_1step_draws(draws_r, Y_mat, x, z, J, j_star, P, Q, "raw",      T_train, T_test, S_sub = CONFIG$S_sub, seed=CONFIG$seed+4)

elpd_c  <- elpd_one_step(y_te, dr_c$mu_draws, dr_c$phi_draws)
elpd_r  <- elpd_one_step(y_te, dr_r$mu_draws, dr_r$phi_draws)
cov95_c <- predictive_coverage_y(y_te, dr_c$mu_draws, dr_c$phi_draws, 0.95, seed=CONFIG$seed+5)
cov95_r <- predictive_coverage_y(y_te, dr_r$mu_draws, dr_r$phi_draws, 0.95, seed=CONFIG$seed+6)

scalar <- tibble::tibble(
  model = c("Centered MA","Raw MA"),
  ELPD_1step = c(elpd_c, elpd_r),
  Pred_Cov95 = c(cov95_c, cov95_r)
)
readr::write_csv(m_c |> dplyr::mutate(model="Centered MA"),
                 file.path(root_outdir, "tables", "metrics_centered_weekly.csv"))
readr::write_csv(m_r |> dplyr::mutate(model="Raw MA"),
                 file.path(root_outdir, "tables", "metrics_raw_weekly.csv"))
readr::write_csv(scalar, file.path(root_outdir, "tables", "metrics_scalar_weekly.csv"))

# quick comparison figure for Total row only
totals <- dplyr::bind_rows(m_c, m_r) |>
  dplyr::filter(component == "Total") |>
  dplyr::mutate(model = rep(c("Centered MA","Raw MA"), each=1)) |>
  dplyr::left_join(scalar, by="model")

p_bar <- ggplot2::ggplot(totals, aes(x = model, y = RMSE)) +
  ggplot2::geom_col(alpha = 0.85) +
  ggplot2::geom_text(aes(label = sprintf("MAE=%.5f\nELPD=%.1f\nCov95=%.1f%%", MAE, ELPD_1step, 100*Pred_Cov95)),
                     vjust = -0.6, size = 3.2) +
  ggplot2::labs(title = "Weekly H.8 composition — Out-of-sample 1-step performance (Total)",
                subtitle = paste0("Train: ", Y$date[1], "–", Y$date[T_train], " | Test: ",
                                  Y$date[T_train+1], "–", Y$date[T_total]),
                y = "Overall RMSE", x = NULL,
                caption = "Dirichlet predictive ELPD and 95% coverage shown as annotations") +
  ggplot2::theme_minimal(base_size = 12)
ggplot2::ggsave(file.path(root_outdir, "figs", "weekly_comparison_total.png"), p_bar, width = 8, height = 5, dpi = 200)

## =========================================================
## 9) Rolling 1-step evaluation over the last N weeks (avg ELPD diff)
## =========================================================

if (isTRUE(CONFIG$rolling$do)) {
  message("Starting rolling 1-step evaluation ...")
  
  # Helper: build standardized z for a given origin t0 (train = 1..t0, predict t0+1)
  make_z_sub <- function(t0) {
    T_sub <- t0 + 1
    z_mean <- mean(z_vol[1:t0], na.rm = TRUE)
    z_sd   <- sd(z_vol[1:t0], na.rm = TRUE); if (!is.finite(z_sd) || z_sd <= 0) z_sd <- 1
    cbind(1, (z_vol[1:T_sub] - z_mean)/z_sd)
  }
  
  min_train <- max(52, CONFIG$rolling$min_train)
  end_idx   <- nrow(Y_mat) - 1
  start_idx <- max(min_train, end_idx - CONFIG$rolling$weeks + 1)
  origins   <- seq(start_idx, end_idx, by = max(1L, CONFIG$rolling$step))
  
  # Storage
  rows_tot <- list()
  rows_elpd <- list()
  rows_diag <- list()
  
  for (t0 in origins) {
    date_pred <- Y$date[t0 + 1]
    message(sprintf("  Origin t=%d -> predicting %s", t0, as.character(date_pred)))
    
    # Subsets
    y_tr0 <- Y_mat[1:t0, , drop=FALSE]
    y_te1 <- Y_mat[t0+1, , drop=FALSE]
    
    x_sub <- matrix(1, nrow = t0+1, ncol = 1)
    z_sub <- make_z_sub(t0)
    
    # Fit both (lighter sampler)
    fc <- fit_bayesian_checked(
      sm_centered, y_tr0, x_sub[1:t0,,drop=FALSE], z_sub[1:t0,,drop=FALSE],
      J=J, j_star=j_star, P=P, Q=Q,
      iter=CONFIG$rolling$iter, warmup=CONFIG$rolling$warmup, chains=CONFIG$rolling$chains,
      adapt_delta=CONFIG$rolling$adapt_delta, max_treedepth=CONFIG$rolling$max_treedepth,
      seed=CONFIG$seed + t0, label=sprintf("Centered@%s", as.character(date_pred))
    )
    fr <- fit_bayesian_checked(
      sm_raw, y_tr0, x_sub[1:t0,,drop=FALSE], z_sub[1:t0,,drop=FALSE],
      J=J, j_star=j_star, P=P, Q=Q,
      iter=CONFIG$rolling$iter, warmup=CONFIG$rolling$warmup, chains=CONFIG$rolling$chains,
      adapt_delta=CONFIG$rolling$adapt_delta, max_treedepth=CONFIG$rolling$max_treedepth,
      seed=CONFIG$seed + 1000 + t0, label=sprintf("Raw@%s", as.character(date_pred))
    )
    
    drc <- rstan::extract(fc$fit, permuted=TRUE)
    drr <- rstan::extract(fr$fit, permuted=TRUE)
    
    mu_c <- forecast_1step_mean(drc, Y_mat[1:(t0+1), , drop=FALSE], x_sub, z_sub, J, j_star, P, Q, "centered", T_train=t0, T_test=1)
    mu_r <- forecast_1step_mean(drr, Y_mat[1:(t0+1), , drop=FALSE], x_sub, z_sub, J, j_star, P, Q, "raw",      T_train=t0, T_test=1)
    
    m_c1 <- metrics_overall(y_te1, mu_c) |> dplyr::filter(component=="Total") |> dplyr::mutate(model="Centered MA")
    m_r1 <- metrics_overall(y_te1, mu_r) |> dplyr::filter(component=="Total") |> dplyr::mutate(model="Raw MA")
    rows_tot[[length(rows_tot)+1]] <- dplyr::bind_rows(m_c1, m_r1) |> dplyr::mutate(date = date_pred)
    
    # Draw-based ELPD and coverage (fast subsample)
    drC <- forecast_1step_draws(drc, Y_mat[1:(t0+1), , drop=FALSE], x_sub, z_sub, J, j_star, P, Q, "centered", T_train=t0, T_test=1, S_sub = CONFIG$rolling$S_sub, seed=CONFIG$seed+3+t0)
    drR <- forecast_1step_draws(drr, Y_mat[1:(t0+1), , drop=FALSE], x_sub, z_sub, J, j_star, P, Q, "raw",      T_train=t0, T_test=1, S_sub = CONFIG$rolling$S_sub, seed=CONFIG$seed+4+t0)
    elpd_c  <- elpd_one_step(y_te1, drC$mu_draws, drC$phi_draws)
    elpd_r  <- elpd_one_step(y_te1, drR$mu_draws, drR$phi_draws)
    cov95_c <- predictive_coverage_y(y_te1, drC$mu_draws, drC$phi_draws, 0.95, seed=CONFIG$seed+5+t0)
    cov95_r <- predictive_coverage_y(y_te1, drR$mu_draws, drR$phi_draws, 0.95, seed=CONFIG$seed+6+t0)
    
    rows_elpd[[length(rows_elpd)+1]] <- tibble::tibble(
      date = date_pred,
      ELPD_centered = elpd_c,
      ELPD_raw      = elpd_r,
      ELPD_diff     = elpd_c - elpd_r,
      Cov95_centered = cov95_c,
      Cov95_raw      = cov95_r
    )
    
    rows_diag[[length(rows_diag)+1]] <- tibble::tibble(
      date = date_pred,
      model = c("Centered MA","Raw MA"),
      divergences = c(fc$diag$divergences, fr$diag$divergences),
      max_treedepth_hits = c(fc$diag$max_treedepth_hits, fr$diag$max_treedepth_hits),
      n_rhat_gt = c(fc$diag$n_rhat_gt, fr$diag$n_rhat_gt),
      n_low_ess_bulk = c(fc$diag$n_low_ess_bulk, fr$diag$n_low_ess_bulk),
      attempts = c(fc$attempts, fr$attempts)
    )
    
    gc()
  }
  
  rolling_totals <- dplyr::bind_rows(rows_tot) |> dplyr::arrange(date)
  rolling_elpd   <- dplyr::bind_rows(rows_elpd) |> dplyr::arrange(date)
  rolling_diag   <- dplyr::bind_rows(rows_diag) |> dplyr::arrange(date)
  
  readr::write_csv(rolling_totals, file.path(root_outdir, "tables", "rolling_per_origin_totals.csv"))
  readr::write_csv(rolling_elpd,   file.path(root_outdir, "tables", "rolling_elpd_cov95.csv"))
  readr::write_csv(rolling_diag,   file.path(root_outdir, "tables", "rolling_diagnostics.csv"))
  
  # Summary: average ELPD diff, mean RMSE/MAE over origins (Total component)
  summary_elpd <- rolling_elpd |>
    dplyr::summarise(
      origins = dplyr::n(),
      ELPD_centered_sum = sum(ELPD_centered),
      ELPD_raw_sum      = sum(ELPD_raw),
      ELPD_diff_sum     = sum(ELPD_diff),
      ELPD_diff_mean    = mean(ELPD_diff),
      ELPD_diff_sd      = sd(ELPD_diff)
    )
  
  summary_point <- rolling_totals |>
    dplyr::mutate(model = gsub("\\s+", ".", model)) |>
    dplyr::select(date, model, RMSE, MAE) |>
    tidyr::pivot_wider(
      names_from = model,
      values_from = c(RMSE, MAE),
      names_sep = "_"
    ) |>
    dplyr::summarise(
      origins = dplyr::n(),
      RMSE_centered_mean = mean(RMSE_Centered.MA, na.rm = TRUE),
      RMSE_raw_mean      = mean(RMSE_Raw.MA,      na.rm = TRUE),
      MAE_centered_mean  = mean(MAE_Centered.MA,  na.rm = TRUE),
      MAE_raw_mean       = mean(MAE_Raw.MA,       na.rm = TRUE)
    )
  
  
  readr::write_csv(summary_elpd, file.path(root_outdir, "tables", "rolling_elpd_summary.csv"))
  readr::write_csv(summary_point, file.path(root_outdir, "tables", "rolling_point_summary.csv"))
  
  # Figure: cumulative ELPD diff (Centered - Raw)
  rolling_elpd$cum_diff <- cumsum(rolling_elpd$ELPD_diff)
  p_cum <- ggplot2::ggplot(rolling_elpd, aes(x = date, y = cum_diff)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Rolling 1-step: Cumulative ELPD Difference (Centered − Raw)",
                  y = "Cumulative ELPD diff", x = NULL,
                  caption = "Positive favors Centered-MA") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(root_outdir, "figs", "rolling_elpd_cumdiff.png"), p_cum, width = 10, height = 5, dpi = 200)
  
  # Figure: Total RMSE per origin
  rmse_plot <- rolling_totals |>
    dplyr::select(date, model, RMSE) |>
    ggplot2::ggplot(aes(x = date, y = RMSE, color = model)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Rolling 1-step: Total RMSE per origin",
                  y = "RMSE (Total share)", x = NULL) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(root_outdir, "figs", "rolling_total_rmse.png"), rmse_plot, width = 10, height = 5, dpi = 200)
  
  message("Rolling evaluation complete.")
}

message("Done. Outputs in: ", root_outdir)
