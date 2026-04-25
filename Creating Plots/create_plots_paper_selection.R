source(paste0(config$cobase_dir, "/Creating Plots/Utilities/dm_plots.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/load_data_util.R"))

library(xtable)


fig_path <- function(name, suffix) {
    paste0(config$cobase_dir, config$results_folder,
           "Figures/", name, suffix, ".pdf")
}

tab_path <- function(name, suffix) {
    paste0(config$cobase_dir, config$results_folder,
           "Figures/", name, suffix, ".tex")
}


########################################
# Main: produce all figures + tables for one dataset.
########################################
create_all_results <- function(dataset_name, station_info, group_info) {
    suf   <- paste0("_", dataset_name)
    alpha <- 0.25

    ########################################
    # Figure 1 - dmcrps_t2m: CRPS by station, EMOS-Q vs EMOS-R
    ########################################
    print("working on: dmcrps_t2m")
    dmcrps_scores <- lapply(seq_along(group_info$group_names), function(gn) {
        file_name <- paste0(group_info$group_names[gn], "_mout_",
                            group_info$output_dim_standard[gn])
        res <- new.env(); load_score_env(file_name, res)

        crps_scores <- paste0("crps_", as.vector(outer(
            res$observation_columns, res$stations, paste, sep = "-")))

        data_env <- load_and_prepare_data(file_name, config$dmcrps_bm)
        filter_models(data_env$df, config$dmcrps_models, crps_scores)$df %>%
            subset(score %in% crps_scores) %>%
            pivot_longer(cols = starts_with("bootstrap_"),
                         names_to = "bootstrap", values_to = "value") %>%
            mutate(station = gsub("crps_", "", score), score = "crps") %>%
            separate(station, into = c("var", "station"), sep = "-") %>%
            mutate(var = ifelse(var %in% c("obs", "T_DRYB_10"), "T2m", "DPT"),
                   group = group_info$group_names[gn])
    }) %>% bind_rows() %>%
        left_join(station_info, by = c("station" = "st_num")) %>%
        mutate(st_name = factor(st_name, levels = station_info$st_name))

    n_st <- length(unique(dmcrps_scores$st_name))
    pdf(file = fig_path("flos_dmcrps_t2m", suf), width = 12, height = 10)
    print(
        ggplot(dmcrps_scores, aes(st_name, value)) +
            facet_wrap(~var, nrow = 2) +
            geom_rect(aes(xmin = -Inf, xmax = Inf,
                          ymin = qnorm(alpha), ymax = qnorm(1 - alpha)),
                      fill = "gray75", color = "gray75", alpha = alpha) +
            geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
            geom_hline(yintercept = 0, linetype = "dashed",
                       color = "gray10", linewidth = 0.8) +
            theme_bw() +
            xlab("Station") + ylab("DM test statistic") +
            ggtitle("CRPS: EMOS v EMOS-R") +
            theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
                  axis.text.y = element_text(size = 12),
                  text        = element_text(size = 16),
                  plot.margin = margin(t = 10, r = 33, b = 10, l = 45)) +
            (if (n_st > 1) geom_vline(xintercept = seq(1.5, n_st - 0.5, by = 1),
                                      col = "grey50", lty = "dashed",
                                      linewidth = 0.5))
    )
    dev.off()


    ########################################
    # Figure 2 - dmmvsimssh: SimSSh v SimSSh-R
    ########################################
    print("working on: dmmvsimssh")
    save_dm_boxplot(
        prepare_dm_score_df(group_info, config$dmmvsimssh_bm,
                            config$dmmvsimssh_models),
        save_path     = fig_path("flos_dmesvs_simssh", suf),
        xlab          = "Data group",
        title         = "SimSSh v SimSSh-R",
        x_aes         = "group_names_pretty",
        facet_formula = ~score,
        facet_nrow    = NULL,
        facet_ncol    = 2,
        height        = 6
    )


    ########################################
    # Figure 3 - dmesvs_mvpp: MVPP method v COBASE-GCA
    ########################################
    print("working on: dmesvs_mvpp")
    save_dm_boxplot(
        prepare_dm_score_df(group_info, config$dmesvs_mvpp_bm[2],
                            config$dmesvs_mvpp_models),
        save_path    = fig_path("flos_dmesvs_mvpp", suf),
        xlab         = "Reshuffling method",
        title        = "MVPP method v COBASE-GCA",
        vline_breaks = c(1.5, 2.5, 3.5),
        model_levels = c("SSh", "SimSSh", "ECC", "GCA")
    )


    ########################################
    # Figure 4 - dmesvs_parcop: COBASE-* v unshuffled (paired benchmark)
    ########################################
    print("working on: dmesvs_parcop")
    save_dm_boxplot(
        prepare_dm_score_df(group_info, config$dmesvs_parcop_bm,
                            config$dmesvs_parcop_models),
        save_path    = fig_path("flos_dmesvs_parcop", suf),
        xlab         = "Parametric copula method",
        title        = "COBASE-model v model",
        vline_breaks = c(1.5, 2.5, 3.5),
        model_levels = pretty_models(config$dmesvs_parcop_models)
    )


    ########################################
    # Figure 5 - dmesvs_enscop: COBASE-Ens* v unshuffled (paired benchmark)
    ########################################
    print("working on: dmesvs_enscop")
    save_dm_boxplot(
        prepare_dm_score_df(group_info, config$dmesvs_enscop_bm,
                            config$dmesvs_enscop_models),
        save_path    = fig_path("flos_dmesvs_enscop", suf),
        xlab         = "Ensemble copula method",
        title        = "COBASE-model v model (ensemble-based)",
        vline_breaks = c(1.5, 2.5, 3.5),
        model_levels = pretty_models(config$dmesvs_enscop_models)
    )


    ########################################
    # Figure 6 - dmesvs_ensvshist: ensemble-based v history-based COBASE
    ########################################
    print("working on: dmesvs_ensvshist")
    save_dm_boxplot(
        prepare_dm_score_df(group_info, config$dmesvs_ensvshist_bm,
                            config$dmesvs_ensvshist_models),
        save_path    = fig_path("flos_dmesvs_ensvshist", suf),
        xlab         = "Ensemble-based COBASE method",
        title        = "COBASE-Ens v COBASE (history-based)",
        vline_breaks = c(1.5, 2.5, 3.5),
        model_levels = pretty_models(config$dmesvs_ensvshist_models)
    )


    ########################################
    # Figure 7 - dmesvs_ecc: methods v ECC
    ########################################
    if (!is.null(config$dmesvs_ecc_bm)) {
        print("working on: dmesvs_ecc")
        save_dm_boxplot(
            prepare_dm_score_df(group_info, config$dmesvs_ecc_bm,
                                config$dmesvs_ecc_models),
            save_path    = fig_path("flos_dmesvs_ecc", suf),
            xlab         = "Method",
            title        = "Methods vs ECC",
            vline_breaks = c(2.5, 3.5, 4.5, 5.5, 6.5),
            model_levels = pretty_models(config$dmesvs_ecc_models)
        )
    }


    ########################################
    # Figure 7b - dmesvs_ecc_ens: ensemble-based COBASE v ECC
    ########################################
    if (!is.null(config$dmesvs_ecc_ens_bm)) {
        print("working on: dmesvs_ecc_ens")
        save_dm_boxplot(
            prepare_dm_score_df(group_info, config$dmesvs_ecc_ens_bm,
                                config$dmesvs_ecc_ens_models),
            save_path    = fig_path("flos_dmesvs_ecc_ens", suf),
            xlab         = "Method",
            title        = "Ensemble-based COBASE vs ECC",
            vline_breaks = c(1.5, 2.5, 3.5),
            model_levels = pretty_models(config$dmesvs_ecc_ens_models)
        )
    }


    ########################################
    # Figure 8 - dmesvs_ecc_combined: all methods v ECC
    ########################################
    if (!is.null(config$dmesvs_ecc_combined_bm)) {
        print("working on: dmesvs_ecc_combined")
        save_dm_boxplot(
            prepare_dm_score_df(group_info, config$dmesvs_ecc_combined_bm,
                                config$dmesvs_ecc_combined_models),
            save_path    = fig_path("flos_dmesvs_ecc_combined", suf),
            xlab         = "Method",
            title        = "All methods vs ECC",
            vline_breaks = c(2.5, 3.5, 7.5),
            model_levels = pretty_models(config$dmesvs_ecc_combined_models),
            width        = 16
        )
    }


    ########################################
    # Table B.3 - CRPS by station
    ########################################
    crps_models <- c("ens", config$crps_scores_models)

    out_crps <- lapply(seq_along(group_info$group_names), function(gn) {
        file_name <- paste0(group_info$group_names[gn], "_mout_",
                            group_info$output_dim_standard[gn])
        res <- new.env(); load_score_env(file_name, res)

        lapply(crps_models, function(m) {
            tmp <- colMeans(res[["crps_list"]][[m]])
            data.frame(crps = tmp, st_num = names(tmp), model = m,
                       score = "crps_list") %>%
                separate(st_num, into = c("var", "st_num"), sep = "-")
        }) %>% bind_rows() %>%
            mutate(model = pretty_models(model))
    }) %>% bind_rows()

    out_crps_df <- out_crps %>%
        filter(var %in% c("obs", "T_DRYB_10")) %>%
        merge(., station_info) %>%
        mutate(Station = st_name) %>%
        dplyr::select(Station, model, crps) %>%
        distinct(Station, model, .keep_all = TRUE) %>%
        pivot_wider(names_from = model, values_from = crps)

    out_crps_df <- out_crps_df[, c("Station", pretty_models(crps_models[-1]))]
    print(xtable(out_crps_df, digits = 4),
          include.rownames = FALSE, file = tab_path("crps_table", suf))


    ########################################
    # Table B.4 - ES/VS summary
    ########################################
    esvs_models <- c("ens", config$mvpp_scores_models)
    score_kinds <- c("es_list", "vs1_list")

    out_esvs <- lapply(seq_along(group_info$group_names), function(gn) {
        file_name <- paste0(group_info$group_names[gn], "_mout_",
                            group_info$output_dim_standard[gn])
        res <- new.env(); load_score_env(file_name, res)

        lapply(esvs_models, function(m) {
            lapply(score_kinds, function(s) {
                data.frame(value        = mean(res[[s]][[m]]),
                           group        = group_info$group_names[gn],
                           group_pretty = group_info$group_names_pretty[gn],
                           model        = m,
                           score        = ifelse(s == "es_list",
                                                 "Energy Score",
                                                 "Variogram Score"))
            }) %>% bind_rows()
        }) %>% bind_rows() %>%
            mutate(model = pretty_models(model))
    }) %>% bind_rows() %>% unique()

    out_esvs_formatted <- out_esvs %>%
        mutate(score_abbrev = ifelse(score == "Energy Score", "ESS", "VSS"),
               group_score  = paste(group_pretty, score_abbrev, sep = "_")) %>%
        dplyr::select(-score, -score_abbrev, -group, -group_pretty) %>%
        pivot_wider(names_from = group_score, values_from = value)

    print(xtable(out_esvs_formatted, digits = 4),
          include.rownames = FALSE, file = tab_path("esvs_table", suf))
}


########################################
# Run for both datasets.
########################################

# Mock data (synthetic)
mock_station_info <- data.frame(
    st_num  = c("1", "2", "3"),
    st_name = c("St1", "St2", "St3"),
    group   = "Mock_data",
    lat     = NA, lon = NA
)
mock_groups_info <- data.frame(
    group_names         = "Mock_data",
    group_names_pretty  = "Mock Data",
    output_dim_standard = 10
)
create_all_results("Mock_data", mock_station_info, mock_groups_info)


# KNMI data (real)
knmi_station_info <- data.frame(
    st_num  = c("235", "240", "260", "280", "310", "380"),
    st_name = c("De Kooy", "Schiphol", "De Bilt",
                "Leeuwarden", "Vlissingen", "Maastricht"),
    group   = "KNMI",
    lat     = c(52.93, 52.30, 52.10, 53.22, 51.44, 50.91),
    lon     = c(4.78, 4.77, 5.18, 5.75, 3.60, 5.77)
)
knmi_groups_info <- data.frame(
    group_names         = "KNMI",
    group_names_pretty  = "KNMI",
    output_dim_standard = 51
)
create_all_results("KNMI", knmi_station_info, knmi_groups_info)
