source(paste0(config$cobase_dir, "/Postprocessing/Utilities/load_data_util.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/ensfc_util.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/scores_util.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/mvpp_methods.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/DM_util.R"))

postprocess_all <- function(file_name, observation_columns, ensemble_regex, output_dim_standard)
{

    # Fixed settings
    trainingDays            <- config$UVPP_Training_Window + config$MVPP_Training_Window
    trainingWindow          <- config$MVPP_Training_Window

    # set random seed
    set.seed(1)

    # Stations and days from the data
    data            <- load_data(file_name)
    stations        <- unique(data$station)
    days            <- sort(unique(data$td))
    nout            <- length(days) - trainingDays
    d               <- length(stations) * length(observation_columns)
    get_dimension   <- load_dimension_transform(stations, observation_columns)

    # Ensemble members
    m <- sum(grepl(ensemble_regex[1], names(data))) # Each observation should have the same number of ens members

    # Score save path
    score_dir <- paste0(config$cobase_dir, config$results_folder, "/Scores/")
    if (!dir.exists(score_dir)) {
        dir.create(score_dir, recursive = TRUE)
    }
    savename <- paste0(file_name, "_mout_", output_dim_standard)
    score_file <- paste0(score_dir, "score_env_", savename, ".RData")

    # Environment to store scores in
    score.env <- new.env()

    # Try to load existing scores for resume
    if (file.exists(score_file)) {
        load(score_file, envir = score.env)
        cat(sprintf("[RESUME] Loaded existing scores from %s\n", score_file))
        existing_methods <- names(score.env$es_list)
        cat(sprintf("[RESUME] Already computed: %s\n", paste(existing_methods, collapse = ", ")))
    }

    # Save settings (always overwrite these)
    score.env$file_name             <- file_name
    score.env$d                     <- d
    score.env$stations              <- as.character(stations)
    score.env$nout                  <- nout
    score.env$trainingDays          <- trainingDays
    score.env$trainingWindow        <- trainingWindow
    score.env$observation_columns   <- observation_columns
    score.env$ensemble_regex        <- ensemble_regex

    # Load transformed data
    transformed_data                <- load_transformed_data(file_name)
    score.env$obs                   <- transformed_data$obs
    score.env$obs_init              <- transformed_data$obs
    score.env$ensfc                 <- transformed_data$ensfc
    score.env$ensfc_init            <- transformed_data$ensfc_init
    score.env$mvpp_list$ens         <- transformed_data$ensfc

    # Load precomputed datastructures
    sim_matrix                      <- load_sim_matrix(file_name)
    score.env$sim_matrix            <- sim_matrix

    # Helper: check if a method is already computed
    method_done <- function(method_name) {
        !is.null(score.env$es_list[[method_name]])
    }

    # Helper: save scores incrementally after each method
    save_checkpoint <- function() {
        save(
            list = ls(score.env),
            file = score_file,
            envir = score.env
        )
    }

    # Add the ens scores
    if (!method_done("ens")) {
        add_scores(score.env = score.env,
                   res = list(mvppout = transformed_data$ensfc),
                   method = "ens",
                   start_time = 0)
        save_checkpoint()
    }

    # General function to apply a MVPP method
    mvpp_adjusted <- function(...) {
        if ("output_dim" %in% names(list(...))) {
            return(mvpp(...,
                transformed_data  = transformed_data,
                score.env         = score.env
            ))
        }

        if ("trainingWindow" %in% names(list(...))) {
            return(mvpp(...,
                transformed_data  = transformed_data,
                score.env         = score.env,
                output_dim        = output_dim_standard,
                addTrainingDays   = TRUE
            ))
        }

        return(mvpp(...,
            transformed_data  = transformed_data,
            score.env         = score.env,
            trainingWindow    = trainingWindow,
            output_dim        = output_dim_standard
        ))
    }

    # Wrapper that skips already-computed methods and saves after each
    run_method <- function(label, run_fn) {
        if (method_done(label)) {
            cat(sprintf("[SKIP] %s already computed\n", label))
            return(invisible(NULL))
        }
        res <- run_fn()
        save_checkpoint()
        return(res)
    }


    ##################
    ## EMOS for ECC ##
    ##################
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Post-processing: %s (d=%d, n=%d, m=%d, output_dim=%d)\n", file_name, d, nout, m, output_dim_standard))
    cat(sprintf("========================================\n"))

    # EMOS-Q and EMOS-R at full ensemble size (for ECC) -- not scored, always recompute
    emos.q <- mvpp_adjusted(
        method            = "EMOS",
        variant           = "Q",
        output_dim        = m,
        saveScores        = FALSE
    )

    emos.r <- mvpp_adjusted(
        method            = "EMOS",
        variant           = "R",
        output_dim        = m,
        saveScores        = FALSE
    )

    #########
    ## ECC ##
    #########

    run_method("ECC-Q", function() mvpp_adjusted(
        method            = "ECC",
        variant           = "Q",
        EMOS_sample       = emos.q$mvppout
    ))

    run_method("ECC-R", function() mvpp_adjusted(
        method            = "ECC",
        variant           = "R",
        EMOS_sample       = emos.r$mvppout
    ))

    ##########
    ## EMOS ##
    ##########

    run_method("EMOS-Q", function() {
        emos.q <<- mvpp_adjusted(method = "EMOS", variant = "Q")
        emos.q
    })

    run_method("EMOS-R", function() {
        emos.r <<- mvpp_adjusted(method = "EMOS", variant = "R")
        emos.r
    })

    # Recompute EMOS-Q/R at output_dim if they were skipped (needed for downstream methods)
    if (is.null(emos.q) || length(dim(emos.q$mvppout)) == 0 || dim(emos.q$mvppout)[2] != output_dim_standard) {
        emos.q <- mvpp_adjusted(method = "EMOS", variant = "Q", saveScores = FALSE)
    }
    if (is.null(emos.r) || length(dim(emos.r$mvppout)) == 0 || dim(emos.r$mvppout)[2] != output_dim_standard) {
        emos.r <- mvpp_adjusted(method = "EMOS", variant = "R", saveScores = FALSE)
    }

    #########
    ## SSH ##
    #########

    run_method("SSh-I14-Q", function() mvpp_adjusted(
        method            = "SSh-I14",
        EMOS_sample       = emos.q$mvppout,
        variant           = "Q"
    ))

    run_method("SimSchaake-Q", function() mvpp_adjusted(
        method            = "SimSchaake",
        EMOS_sample       = emos.q$mvppout,
        sim_matrix        = sim_matrix,
        variant           = "Q"
    ))

    run_method("SimSchaake-R", function() mvpp_adjusted(
        method            = "SimSchaake",
        EMOS_sample       = emos.r$mvppout,
        sim_matrix        = sim_matrix,
        variant           = "R"
    ))

    #########
    ## GCA ##
    #########

    gca.cop <- NULL
    run_method("CopGCA", function() {
        gca.cop <<- mvpp_adjusted(method = "CopGCA")
        gca.cop
    })

    run_method("COBASE-CopGCA-Q", function() {
        if (is.null(gca.cop)) gca.cop <<- mvpp_adjusted(method = "CopGCA", saveScores = FALSE)
        mvpp_adjusted(
            method            = "CopGCA",
            EMOS_sample       = emos.q$mvppout,
            shuffle           = TRUE,
            MVPP_sample       = gca.cop$mvppout,
            variant           = "Q"
        )
    })

    #########################
    ## Archimedean Copulas ##
    #########################

    clayton <- NULL
    run_method("Clayton", function() {
        clayton <<- mvpp_adjusted(method = "Clayton")
        clayton
    })

    gumbel <- NULL
    run_method("Gumbel", function() {
        gumbel <<- mvpp_adjusted(method = "Gumbel")
        gumbel
    })

    frank <- NULL
    run_method("Frank", function() {
        frank <<- mvpp_adjusted(method = "Frank")
        frank
    })

    run_method("COBASE-Clayton-Q", function() {
        if (is.null(clayton)) clayton <<- mvpp_adjusted(method = "Clayton", saveScores = FALSE)
        mvpp_adjusted(
            method            = "Clayton",
            EMOS_sample       = emos.q$mvppout,
            shuffle           = TRUE,
            MVPP_sample       = clayton$mvppout,
            variant           = "Q"
        )
    })

    run_method("COBASE-Gumbel-Q", function() {
        if (is.null(gumbel)) gumbel <<- mvpp_adjusted(method = "Gumbel", saveScores = FALSE)
        mvpp_adjusted(
            method            = "Gumbel",
            EMOS_sample       = emos.q$mvppout,
            shuffle           = TRUE,
            MVPP_sample       = gumbel$mvppout,
            variant           = "Q"
        )
    })

    run_method("COBASE-Frank-Q", function() {
        if (is.null(frank)) frank <<- mvpp_adjusted(method = "Frank", saveScores = FALSE)
        mvpp_adjusted(
            method            = "Frank",
            EMOS_sample       = emos.q$mvppout,
            shuffle           = TRUE,
            MVPP_sample       = frank$mvppout,
            variant           = "Q"
        )
    })

    #####################################
    ## Ensemble-based Copula Variants  ##
    #####################################

    ens.cop.gca <- NULL
    run_method("EnsCopGCA", function() {
        ens.cop.gca <<- mvpp_adjusted(method = "EnsCopGCA")
        ens.cop.gca
    })

    ens.clayton <- NULL
    run_method("EnsClayton", function() {
        ens.clayton <<- mvpp_adjusted(method = "EnsClayton")
        ens.clayton
    })

    ens.gumbel <- NULL
    run_method("EnsGumbel", function() {
        ens.gumbel <<- mvpp_adjusted(method = "EnsGumbel")
        ens.gumbel
    })

    ens.frank <- NULL
    run_method("EnsFrank", function() {
        ens.frank <<- mvpp_adjusted(method = "EnsFrank")
        ens.frank
    })

    # COBASE shuffled variants
    run_method("COBASE-EnsCopGCA-Q", function() {
        if (is.null(ens.cop.gca)) ens.cop.gca <<- mvpp_adjusted(method = "EnsCopGCA", saveScores = FALSE)
        mvpp_adjusted(
            method            = "EnsCopGCA",
            EMOS_sample       = emos.q$mvppout,
            shuffle           = TRUE,
            MVPP_sample       = ens.cop.gca$mvppout,
            variant           = "Q"
        )
    })

    run_method("COBASE-EnsClayton-Q", function() {
        if (is.null(ens.clayton)) ens.clayton <<- mvpp_adjusted(method = "EnsClayton", saveScores = FALSE)
        mvpp_adjusted(
            method            = "EnsClayton",
            EMOS_sample       = emos.q$mvppout,
            shuffle           = TRUE,
            MVPP_sample       = ens.clayton$mvppout,
            variant           = "Q"
        )
    })

    run_method("COBASE-EnsGumbel-Q", function() {
        if (is.null(ens.gumbel)) ens.gumbel <<- mvpp_adjusted(method = "EnsGumbel", saveScores = FALSE)
        mvpp_adjusted(
            method            = "EnsGumbel",
            EMOS_sample       = emos.q$mvppout,
            shuffle           = TRUE,
            MVPP_sample       = ens.gumbel$mvppout,
            variant           = "Q"
        )
    })

    run_method("COBASE-EnsFrank-Q", function() {
        if (is.null(ens.frank)) ens.frank <<- mvpp_adjusted(method = "EnsFrank", saveScores = FALSE)
        mvpp_adjusted(
            method            = "EnsFrank",
            EMOS_sample       = emos.q$mvppout,
            shuffle           = TRUE,
            MVPP_sample       = ens.frank$mvppout,
            variant           = "Q"
        )
    })

    ###########################
    ## Compute DM statistics ##
    ###########################
    print("In Postprocessing/Utilities/mvpp.R: computing DM scores...")
    benchmarks = c(config$dmcrps_bm, config$dmmvsimssh_bm, config$dmesvs_mvpp_bm, config$dmesvs_parcop_bm, config$dmesvs_enscop_bm, config$dmesvs_ensvshist_bm, config$dmesvs_ecc_bm)
    compute_DM_scores(file_name = savename,
                      benchmarks = benchmarks,
                      parallelization = TRUE)

}
