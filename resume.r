# Resume from where we left off:
# - Mock_data: fully done
# - KNMI: univariate done, need multivariate + DM + plots

config <- config::get()

# Skip univariate (already done for both datasets)

# Run KNMI multivariate only (Mock_data already done)
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/mvpp.R"))

postprocess_all(
    file_name               = "KNMI",
    observation_columns     = c("T_DRYB_10", "T_DEWP_10"),
    ensemble_regex          = c("^IFS_x2t_", "^IFS_x2d_"),
    output_dim_standard     = 51
)

# Generate plots for both datasets
source(paste0(config$cobase_dir, "/Creating Plots/create_plots_paper_selection.R"))
