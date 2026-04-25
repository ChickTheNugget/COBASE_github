library(ggplot2)
library(dplyr)
library(tidyr)
library(copula)

if (!exists("config")) config <- config::get()

results_dir <- paste0(config$cobase_dir, config$results_folder)
fig_dir     <- paste0(results_dir, "Figures/ParamComparison/")
table_dir   <- paste0(results_dir, "Tables/ParamComparison/")
if (!dir.exists(fig_dir))   dir.create(fig_dir,   recursive = TRUE)
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

FAMILIES <- c("Clayton", "Frank", "Gumbel")

# Convert a vector of theta values to Kendall's tau using the copula package
theta_to_tau <- function(family, theta) {
    cop_fn <- switch(family,
        "Clayton" = claytonCopula,
        "Frank"   = frankCopula,
        "Gumbel"  = gumbelCopula
    )
    sapply(theta, function(p) {
        if (is.na(p)) return(NA_real_)
        tryCatch(tau(cop_fn(param = p)), error = function(e) NA_real_)
    })
}

# Build a long-format data frame of paired (obs, ens) thetas per day per family
build_pair_df <- function(score_file) {
    env <- new.env()
    load(score_file, envir = env)
    p <- env$params

    rows <- lapply(FAMILIES, function(fam) {
        obs <- p[[fam]]
        ens <- p[[paste0("Ens", fam)]]
        n   <- length(obs)
        data.frame(
            family    = fam,
            day       = seq_len(n),
            theta_obs = as.numeric(obs),
            theta_ens = as.numeric(ens)
        )
    })
    df <- do.call(rbind, rows)
    df$family <- factor(df$family, levels = FAMILIES)

    df$tau_obs <- NA_real_
    df$tau_ens <- NA_real_
    for (fam in FAMILIES) {
        sel <- df$family == fam
        df$tau_obs[sel] <- theta_to_tau(fam, df$theta_obs[sel])
        df$tau_ens[sel] <- theta_to_tau(fam, df$theta_ens[sel])
    }
    df
}

# ---------- Plot helpers ----------

plot_theme <- theme_bw(base_size = 12) +
    theme(
        strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank(),
        legend.position  = "top"
    )

scatter_plot <- function(df, dataset_label) {
    paired <- df[stats::complete.cases(df[, c("theta_obs", "theta_ens")]), ]
    cor_lab <- paired %>%
        group_by(family) %>%
        summarise(
            n   = dplyr::n(),
            r_p = cor(theta_obs, theta_ens, method = "pearson"),
            r_s = cor(theta_obs, theta_ens, method = "spearman"),
            .groups = "drop"
        ) %>%
        mutate(label = sprintf("n = %d\nPearson r = %.2f\nSpearman = %.2f", n, r_p, r_s))

    ggplot(paired, aes(x = theta_obs, y = theta_ens)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
        geom_point(alpha = 0.45, size = 1.1) +
        geom_text(
            data = cor_lab,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE
        ) +
        facet_wrap(~ family, scales = "free", nrow = 1) +
        labs(
            title = sprintf("Scatter: theta_ens vs theta_obs (%s)", dataset_label),
            x = expression(theta[obs]),
            y = expression(theta[ens])
        ) +
        plot_theme
}

bland_altman_plot <- function(df, dataset_label) {
    paired <- df[stats::complete.cases(df[, c("theta_obs", "theta_ens")]), ]
    paired$avg  <- (paired$theta_obs + paired$theta_ens) / 2
    paired$diff <- paired$theta_ens - paired$theta_obs

    band <- paired %>%
        group_by(family) %>%
        summarise(
            bias  = mean(diff),
            sd_d  = sd(diff),
            lower = mean(diff) - 1.96 * sd(diff),
            upper = mean(diff) + 1.96 * sd(diff),
            .groups = "drop"
        ) %>%
        mutate(label = sprintf("bias = %.3f\n+/- 1.96 SD = [%.3f, %.3f]", bias, lower, upper))

    ggplot(paired, aes(x = avg, y = diff)) +
        geom_hline(yintercept = 0, colour = "grey60") +
        geom_point(alpha = 0.45, size = 1.1) +
        geom_hline(data = band, aes(yintercept = bias),  colour = "red", linewidth = 0.6) +
        geom_hline(data = band, aes(yintercept = lower), colour = "red", linetype = "dashed") +
        geom_hline(data = band, aes(yintercept = upper), colour = "red", linetype = "dashed") +
        geom_text(
            data = band,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.2, size = 3, inherit.aes = FALSE
        ) +
        facet_wrap(~ family, scales = "free", nrow = 1) +
        labs(
            title    = sprintf("Bland-Altman: ensemble vs observation fits (%s)", dataset_label),
            subtitle = "Solid red = mean bias; dashed red = 95% limits of agreement (mean +/- 1.96 SD).",
            x = expression((theta[obs] + theta[ens]) / 2),
            y = expression(theta[ens] - theta[obs])
        ) +
        plot_theme
}

distribution_plot <- function(df, dataset_label) {
    long <- df %>%
        select(family, day, theta_obs, theta_ens) %>%
        pivot_longer(c(theta_obs, theta_ens), names_to = "method", values_to = "theta") %>%
        mutate(method = factor(
            recode(method, theta_obs = "obs-based", theta_ens = "ens-based"),
            levels = c("obs-based", "ens-based")
        )) %>%
        filter(!is.na(theta))

    ggplot(long, aes(x = method, y = theta, fill = method)) +
        geom_boxplot(outlier.alpha = 0.4, width = 0.55) +
        facet_wrap(~ family, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = c("obs-based" = "#1f77b4", "ens-based" = "#ff7f0e")) +
        labs(
            title    = sprintf("Distribution of theta over the test period (%s)", dataset_label),
            subtitle = "Each box is one method's distribution of theta across all forecast days.",
            x = NULL, y = expression(theta), fill = NULL
        ) +
        plot_theme
}

summary_table <- function(df) {
    df %>%
        group_by(family) %>%
        summarise(
            n_days        = dplyr::n(),
            na_obs_pct    = 100 * mean(is.na(theta_obs)),
            na_ens_pct    = 100 * mean(is.na(theta_ens)),
            mean_theta_obs = mean(theta_obs, na.rm = TRUE),
            mean_theta_ens = mean(theta_ens, na.rm = TRUE),
            sd_theta_obs   = sd(theta_obs,   na.rm = TRUE),
            sd_theta_ens   = sd(theta_ens,   na.rm = TRUE),
            mean_tau_obs   = mean(tau_obs,   na.rm = TRUE),
            mean_tau_ens   = mean(tau_ens,   na.rm = TRUE),
            sd_tau_obs     = sd(tau_obs,     na.rm = TRUE),
            sd_tau_ens     = sd(tau_ens,     na.rm = TRUE),
            mean_diff_theta = mean(theta_ens - theta_obs, na.rm = TRUE),
            mean_abs_diff_theta = mean(abs(theta_ens - theta_obs), na.rm = TRUE),
            .groups = "drop"
        )
}

# ---------- Gaussian copula helpers (params is a list of correlation vectors) ----------

# Vectorised theta -> tau for a Gaussian copula: tau = (2/pi) * arcsin(rho)
theta_to_tau_gaussian <- function(rho) {
    (2 / pi) * asin(pmin(pmax(rho, -1), 1))
}

# Build a long-format data frame of paired (obs, ens) correlations per (day, pair).
# CopGCA / EnsCopGCA store params as a list of length n; each element is a
# numeric vector of length d*(d-1)/2 (upper triangle of the correlation matrix)
# or NULL if the fit fell back to indepCopula.
build_gaussian_pair_df <- function(score_file) {
    env <- new.env()
    load(score_file, envir = env)
    p <- env$params

    obs <- p[["CopGCA"]]
    ens <- p[["EnsCopGCA"]]
    if (is.null(obs) || is.null(ens)) {
        return(data.frame(day = integer(0), pair_idx = integer(0),
                          rho_obs = numeric(0), rho_ens = numeric(0),
                          tau_obs = numeric(0), tau_ens = numeric(0)))
    }

    n <- max(length(obs), length(ens))
    rows <- vector("list", n)
    for (nn in seq_len(n)) {
        v_obs <- obs[[nn]]
        v_ens <- ens[[nn]]
        if (is.null(v_obs) || is.null(v_ens)) next
        if (length(v_obs) != length(v_ens)) next
        rows[[nn]] <- data.frame(
            day      = nn,
            pair_idx = seq_along(v_obs),
            rho_obs  = as.numeric(v_obs),
            rho_ens  = as.numeric(v_ens)
        )
    }
    df <- do.call(rbind, rows)
    if (is.null(df) || nrow(df) == 0) {
        return(data.frame(day = integer(0), pair_idx = integer(0),
                          rho_obs = numeric(0), rho_ens = numeric(0),
                          tau_obs = numeric(0), tau_ens = numeric(0)))
    }
    df$tau_obs <- theta_to_tau_gaussian(df$rho_obs)
    df$tau_ens <- theta_to_tau_gaussian(df$rho_ens)
    df
}

gaussian_scatter_plot <- function(df, dataset_label) {
    paired <- df[stats::complete.cases(df[, c("rho_obs", "rho_ens")]), ]
    cor_lab <- data.frame(
        n   = nrow(paired),
        r_p = cor(paired$rho_obs, paired$rho_ens, method = "pearson"),
        r_s = cor(paired$rho_obs, paired$rho_ens, method = "spearman")
    )
    cor_lab$label <- sprintf("n = %d\nPearson r = %.2f\nSpearman = %.2f",
                             cor_lab$n, cor_lab$r_p, cor_lab$r_s)

    ggplot(paired, aes(x = rho_obs, y = rho_ens)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
        geom_point(alpha = 0.08, size = 0.6) +
        geom_text(
            data = cor_lab,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE
        ) +
        labs(
            title = sprintf("Gaussian rho_ens vs rho_obs (%s)", dataset_label),
            x = expression(rho[obs]),
            y = expression(rho[ens])
        ) +
        plot_theme
}

gaussian_bland_altman_plot <- function(df, dataset_label) {
    paired <- df[stats::complete.cases(df[, c("rho_obs", "rho_ens")]), ]
    paired$avg  <- (paired$rho_obs + paired$rho_ens) / 2
    paired$diff <- paired$rho_ens - paired$rho_obs

    band <- data.frame(
        bias  = mean(paired$diff),
        sd_d  = sd(paired$diff),
        lower = mean(paired$diff) - 1.96 * sd(paired$diff),
        upper = mean(paired$diff) + 1.96 * sd(paired$diff)
    )
    band$label <- sprintf("bias = %.3f\n+/- 1.96 SD = [%.3f, %.3f]",
                          band$bias, band$lower, band$upper)

    ggplot(paired, aes(x = avg, y = diff)) +
        geom_hline(yintercept = 0, colour = "grey60") +
        geom_point(alpha = 0.08, size = 0.6) +
        geom_hline(data = band, aes(yintercept = bias),  colour = "red", linewidth = 0.6) +
        geom_hline(data = band, aes(yintercept = lower), colour = "red", linetype = "dashed") +
        geom_hline(data = band, aes(yintercept = upper), colour = "red", linetype = "dashed") +
        geom_text(
            data = band,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.2, size = 3, inherit.aes = FALSE
        ) +
        labs(
            title    = sprintf("Gaussian Bland-Altman (%s)", dataset_label),
            subtitle = "Solid red = mean bias; dashed red = 95% limits of agreement.",
            x = expression((rho[obs] + rho[ens]) / 2),
            y = expression(rho[ens] - rho[obs])
        ) +
        plot_theme
}

gaussian_distribution_plot <- function(df, dataset_label) {
    long <- df %>%
        select(day, pair_idx, rho_obs, rho_ens) %>%
        pivot_longer(c(rho_obs, rho_ens), names_to = "method", values_to = "rho") %>%
        mutate(method = factor(
            recode(method, rho_obs = "obs-based", rho_ens = "ens-based"),
            levels = c("obs-based", "ens-based")
        )) %>%
        filter(!is.na(rho))

    ggplot(long, aes(x = method, y = rho, fill = method)) +
        geom_boxplot(outlier.alpha = 0.4, width = 0.55) +
        scale_fill_manual(values = c("obs-based" = "#1f77b4", "ens-based" = "#ff7f0e")) +
        labs(
            title    = sprintf("Gaussian rho distribution (%s)", dataset_label),
            subtitle = "Pooled across all (day, station-pair) entries.",
            x = NULL, y = expression(rho), fill = NULL
        ) +
        plot_theme
}

# Per-day Frobenius distance between R_obs and R_ens.
# Only the d*(d-1)/2 unique off-diagonal entries are stored, so
# ||R_obs - R_ens||_F = sqrt(2 * sum (rho_obs - rho_ens)^2) (diagonal cancels).
gaussian_matrix_distance_plot <- function(df, dataset_label) {
    paired <- df[stats::complete.cases(df[, c("rho_obs", "rho_ens")]), ]
    per_day <- paired %>%
        group_by(day) %>%
        summarise(frob = sqrt(2 * sum((rho_obs - rho_ens)^2)), .groups = "drop")

    mean_frob <- mean(per_day$frob)

    ggplot(per_day, aes(x = day, y = frob)) +
        geom_line(colour = "grey50", alpha = 0.6) +
        geom_point(size = 0.6, alpha = 0.6) +
        geom_hline(yintercept = mean_frob, colour = "red", linewidth = 0.6) +
        annotate("text", x = -Inf, y = Inf,
                 label = sprintf("mean = %.3f", mean_frob),
                 hjust = -0.1, vjust = 1.2, size = 3) +
        labs(
            title    = sprintf("Gaussian per-day Frobenius distance (%s)", dataset_label),
            subtitle = "One value per forecast day; red line is the mean.",
            x = "Forecast day index",
            y = expression("|| R"[obs] * " - R"[ens] * " ||"[F])
        ) +
        plot_theme
}

gaussian_summary_table <- function(df) {
    paired <- df[stats::complete.cases(df[, c("rho_obs", "rho_ens")]), ]
    if (nrow(paired) == 0) return(data.frame())

    per_day <- paired %>%
        group_by(day) %>%
        summarise(frob = sqrt(2 * sum((rho_obs - rho_ens)^2)),
                  n_pairs = dplyr::n(), .groups = "drop")

    data.frame(
        n_days              = length(unique(paired$day)),
        n_pairs_per_day     = stats::median(per_day$n_pairs),
        mean_rho_obs        = mean(paired$rho_obs),
        mean_rho_ens        = mean(paired$rho_ens),
        sd_rho_obs          = sd(paired$rho_obs),
        sd_rho_ens          = sd(paired$rho_ens),
        mean_tau_obs        = mean(paired$tau_obs),
        mean_tau_ens        = mean(paired$tau_ens),
        mean_diff_rho       = mean(paired$rho_ens - paired$rho_obs),
        mean_abs_diff_rho   = mean(abs(paired$rho_ens - paired$rho_obs)),
        mean_frob_distance  = mean(per_day$frob)
    )
}

# ---------- Main loop ----------

run_for_dataset <- function(dataset, output_dim, label) {
    score_file <- paste0(results_dir, "Scores/score_env_", dataset, "_mout_", output_dim, ".RData")
    cat(sprintf("[ParamComparison] %s -> %s\n", label, score_file))

    df <- build_pair_df(score_file)

    p_sc <- scatter_plot(df, label)
    p_ba <- bland_altman_plot(df, label)
    p_di <- distribution_plot(df, label)
    tab  <- summary_table(df)

    ggsave(paste0(fig_dir, "scatter_",      dataset, ".pdf"), p_sc, width = 11, height = 4)
    ggsave(paste0(fig_dir, "blandaltman_",  dataset, ".pdf"), p_ba, width = 11, height = 4)
    ggsave(paste0(fig_dir, "distribution_", dataset, ".pdf"), p_di, width = 11, height = 4)

    write.csv(tab, paste0(table_dir, "summary_", dataset, ".csv"), row.names = FALSE)

    cat("Archimedean summary:\n")
    print(tab, width = Inf)

    # Gaussian
    df_g <- build_gaussian_pair_df(score_file)
    if (nrow(df_g) > 0) {
        ggsave(paste0(fig_dir, "gaussian_scatter_",      dataset, ".pdf"),
               gaussian_scatter_plot(df_g, label),         width = 6, height = 5)
        ggsave(paste0(fig_dir, "gaussian_blandaltman_",  dataset, ".pdf"),
               gaussian_bland_altman_plot(df_g, label),    width = 6, height = 5)
        ggsave(paste0(fig_dir, "gaussian_distribution_", dataset, ".pdf"),
               gaussian_distribution_plot(df_g, label),    width = 5, height = 4)
        ggsave(paste0(fig_dir, "gaussian_frobenius_",    dataset, ".pdf"),
               gaussian_matrix_distance_plot(df_g, label), width = 8, height = 4)
        gtab <- gaussian_summary_table(df_g)
        write.csv(gtab,
                  paste0(table_dir, "gaussian_summary_", dataset, ".csv"),
                  row.names = FALSE)
        cat("Gaussian summary:\n")
        print(gtab)
    } else {
        cat("[ParamComparison] No Gaussian params found in score env; skipping Gaussian plots.\n")
    }

    invisible(list(df = df, table = tab, df_g = df_g))
}

run_for_dataset("Mock_data", 10, "Mock_data")
run_for_dataset("KNMI",      51, "KNMI")
