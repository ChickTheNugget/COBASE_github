# Recompute DM statistics (with ECC-Q as new benchmark) and regenerate plots
config <- config::get()
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/load_data_util.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/DM_util.R"))

# Recompute DM for both datasets
for (ds in list(
    list(name = "Mock_data", mout = 10),
    list(name = "KNMI", mout = 51)
)) {
    savename <- paste0(ds$name, "_mout_", ds$mout)
    cat(sprintf("Recomputing DM statistics for %s...\n", savename))
    benchmarks <- c(config$dmcrps_bm, config$dmmvsimssh_bm, config$dmesvs_mvpp_bm,
                    config$dmesvs_parcop_bm, config$dmesvs_enscop_bm,
                    config$dmesvs_ensvshist_bm, config$dmesvs_ecc_bm)
    benchmarks <- unique(benchmarks)
    compute_DM_scores(file_name = savename, benchmarks = benchmarks, parallelization = FALSE)
}

# Regenerate all plots
source(paste0(config$cobase_dir, "/Creating Plots/create_plots_paper_selection.R"))
