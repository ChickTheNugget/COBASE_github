source(paste0(config$cobase_dir, "/Postprocessing/Utilities/mvpp.R"))

# Mock data
postprocess_all(
    file_name               = "Mock_data",
    observation_columns     = c("obs"),
    ensemble_regex          = c("^M_"),
    output_dim_standard     = 10
)

# KNMI data
postprocess_all(
    file_name               = "KNMI",
    observation_columns     = c("T_DRYB_10", "T_DEWP_10"),
    ensemble_regex          = c("^IFS_x2t_", "^IFS_x2d_"),
    output_dim_standard     = 51
)
