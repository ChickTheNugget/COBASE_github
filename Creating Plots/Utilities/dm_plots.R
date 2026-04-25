library(ggplot2)
library(tidyr)
library(dplyr)

source(paste0(config$cobase_dir, "/Postprocessing/Utilities/load_data_util.R"))

# ---------------------------------------------------------------------------
# Pretty model names for the paper-selection figures.
# ---------------------------------------------------------------------------
# load_data_util.R already defines pretty_model_name() with a different mapping
# used by other plot utilities; this one is the paper-selection variant.
change_names <- function(model_name) {
    pretty_name <- switch(model_name,
        "ens"                    = "Raw Ensemble",
        "CopGCA"                 = "GCA",
        "SimSchaake-Q"           = "SimSSh",
        "SimSchaake-R"           = "SimSSh-R",
        "ECC-Q"                  = "ECC",
        "SSh-I14-Q"              = "SSh",
        "COBASE-CopGCA-Q"        = "COBASE-GCA",
        "COBASE-Clayton-Q"       = "COBASE-Clayton",
        "COBASE-Frank-Q"         = "COBASE-Frank",
        "COBASE-Gumbel-Q"        = "COBASE-Gumbel",
        "EnsCopGCA"              = "EnsGCA",
        "EnsClayton"             = "EnsClayton",
        "EnsFrank"               = "EnsFrank",
        "EnsGumbel"              = "EnsGumbel",
        "COBASE-EnsCopGCA-Q"     = "COBASE-EnsGCA",
        "COBASE-EnsClayton-Q"    = "COBASE-EnsClayton",
        "COBASE-EnsFrank-Q"      = "COBASE-EnsFrank",
        "COBASE-EnsGumbel-Q"     = "COBASE-EnsGumbel",
        model_name
    )
    return(pretty_name)
}

pretty_models <- function(x) vapply(as.character(x), change_names, character(1))

# ---------------------------------------------------------------------------
# Loaders.
# ---------------------------------------------------------------------------

load_and_prepare_data <- function(file_name, benchmark_name) {
    res <- new.env()
    load_score_env(file_name, res)
    dfmc <- load_dm_statistics(file_name)

    df <- subset(dfmc, benchmark == benchmark_name)
    bootstrap_columns <- grep("bootstrap_[0-9]+$", names(df), value = TRUE)
    df[bootstrap_columns] <- (-1) * df[bootstrap_columns]

    return(list(df = df, res = res))
}

filter_models <- function(df, model_names, input_scores) {
    if (!is.null(model_names)) {
        df <- subset(df, model %in% model_names)
    } else {
        model_names <- unique(df$model)
    }
    return(list(df = df, model_names = model_names, input_scores = input_scores))
}

# ---------------------------------------------------------------------------
# Paper-figure helpers.
# ---------------------------------------------------------------------------

# Build the long-form data frame for an ES/VS DM-statistic boxplot.
# `benchmarks` may be a single value (one benchmark for all models) or a vector
# the same length as `models` (one benchmark per model — paired).
prepare_dm_score_df <- function(group_info, benchmarks, models) {
    score_kinds <- c("es_list", "vs1_list")
    paired <- length(benchmarks) > 1

    df <- lapply(seq_along(group_info$group_names), function(gn) {
        file_name <- paste0(group_info$group_names[gn], "_mout_",
                            group_info$output_dim_standard[gn])

        lapply(score_kinds, function(score_kind) {
            df_raw <- if (paired) {
                bind_rows(lapply(seq_along(benchmarks), function(bb) {
                    data_env <- load_and_prepare_data(file_name, benchmarks[bb])
                    filter_models(data_env$df, models[bb], score_kind)$df
                }))
            } else {
                data_env <- load_and_prepare_data(file_name, benchmarks)
                filter_models(data_env$df, models, score_kind)$df
            }

            df_raw %>%
                subset(score %in% score_kind) %>%
                pivot_longer(cols = starts_with("bootstrap_"),
                             names_to = "bootstrap", values_to = "value") %>%
                mutate(group_name = file_name,
                       group_names_pretty = group_info$group_names_pretty[gn],
                       score = ifelse(score_kind == "es_list",
                                      "Energy Score", "Variogram Score"),
                       model = pretty_models(model),
                       benchmark = pretty_models(benchmark))
        }) %>% bind_rows()
    }) %>% bind_rows()

    df$group_names_pretty <- factor(
        gsub("_mout_[0-9][0-9]", "", df$group_names_pretty),
        levels = group_info$group_names_pretty
    )
    df
}

# Save a DM-statistic boxplot PDF using the recipe shared across paper figures.
save_dm_boxplot <- function(df, save_path, xlab, title,
                            vline_breaks  = NULL,
                            model_levels  = NULL,
                            x_aes         = "model",
                            facet_formula = ~score + group_names_pretty,
                            facet_nrow    = 2,
                            facet_ncol    = NULL,
                            width         = 12,
                            height        = 8) {
    alpha <- 0.25
    if (!is.null(model_levels)) {
        df$model <- factor(df$model, levels = model_levels)
    }

    p <- ggplot(df, aes(.data[[x_aes]], value)) +
        facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol) +
        geom_rect(aes(xmin = -Inf, xmax = Inf,
                      ymin = qnorm(alpha), ymax = qnorm(1 - alpha)),
                  fill = "gray75", color = "gray75", alpha = alpha) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed",
                   color = "gray10", linewidth = 0.8) +
        theme_bw() +
        xlab(xlab) + ylab("DM test statistic") +
        ggtitle(title) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
              axis.text.y = element_text(size = 12),
              text        = element_text(size = 16),
              plot.margin = margin(t = 10, r = 33, b = 10, l = 45))

    if (!is.null(vline_breaks)) {
        p <- p + geom_vline(xintercept = vline_breaks,
                            col = "grey50", lty = "dashed", linewidth = 0.5)
    }

    pdf(file = save_path, width = width, height = height)
    print(p)
    dev.off()
}
