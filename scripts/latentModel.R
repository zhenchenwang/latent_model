# This script performs the latent variable model analysis.
# It reads the ground truth data, performs discretization, learns a network structure,
# enriches the network with latent variables, and generates synthetic data.

# Load libraries
library(bnlearn)
library(pcalg)
library(LaplacesDemon)
library(Rgraphviz)
library(ggplot2)
library(gridExtra)
library(pracma)
library(missForest)
library(gRain)
library(cluster)
library(arules)

# Load configuration
source("scripts/config.R")

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

# Define a function to format variables to their correct types
format_variables <- function(df) {
  # Factor variables
  factor_vars <- c("strokeha", "af", "atyantip", "steroid", "impot", "migr",
                   "ra", "ckidney", "semi", "sle", "treathyp", "type1",
                   "type2", "ethr", "smoking", "fh_cad", "gender", "region")
  for (var in factor_vars) {
    if (var %in% names(df)) df[[var]] <- as.factor(df[[var]])
  }

  # Numeric variables
  numeric_vars <- c("age", "bmi", "choleratio", "sbp", "sbps")
  for (var in numeric_vars) {
    if (var %in% names(df)) df[[var]] <- as.numeric(df[[var]])
  }

  return(df)
}

# Function that transforms the dataset into integer values (from 0..(n-1))
to_integer_values <- function(df) {
  for (k in 1:ncol(df)) {
    df[, k] <- as.integer(df[, k]) - 1
  }
  return(df)
}

# Function to count the number of levels in each factor variable
count_levels <- function(df) {
  sapply(df, nlevels)
}

# Function for discretization of continuous variables using equal intervals
discretize_data <- function(df, discretization_bins) {
  df_discrete <- df

  for (col_name in names(discretization_bins)) {
    if (col_name %in% names(df_discrete) && is.numeric(df_discrete[[col_name]])) {
      num_bins <- discretization_bins[[col_name]]
      df_discrete[[col_name]] <- as.factor(
        arules::discretize(
          df_discrete[[col_name]],
          method = 'interval',
          breaks = num_bins,
          labels = c(0:(num_bins - 1)),
          include.lowest = TRUE,
          right = TRUE
        )
      )
    }
  }

  return(df_discrete)
}

# -----------------------------------------------------------------------------
# Main Function to Generate Synthetic Data
# -----------------------------------------------------------------------------

generate_synthetic_data <- function(ground_truth_sample) {

  # --- Initial Preprocessing ---
  formatted_sample <- format_variables(ground_truth_sample)

  discretization_bins <- list(bmi = 6, age = 4, choleratio = 5, sbp = 5, sbps = 5)
  discretized_sample <- discretize_data(formatted_sample, discretization_bins)

  imputed_sample <- missForest(discretized_sample)$ximp

  # --- Structural Learning from Data (without latent vars) ---
  sem_no_latent <- structural.em(
    formatted_sample,
    maximize = "hc",
    fit = "mle",
    return.all = TRUE,
    start = NULL,
    max.iter = 5
  )

  # --- Enrich the network structure with latent variables ---
  observed_vars <- names(formatted_sample)
  latent_var_names <- paste("L", 1:num_latent_vars, sep = "")

  amat_enriched <- matrix(
    0,
    nrow = length(observed_vars) + num_latent_vars,
    ncol = length(observed_vars) + num_latent_vars,
    dimnames = list(c(observed_vars, latent_var_names), c(observed_vars, latent_var_names))
  )
  amat_enriched[observed_vars, observed_vars] <- amat(sem_no_latent$dag)

  latent_connections <- list(
    L1 = c("age", "af", "treathyp"),
    L2 = c("steroid", "treathyp"),
    L3 = c("impot", "gender"),
    L4 = c("migr", "gender", "choleratio"),
    L5 = c("strokeha", "ckidney", "type2", "choleratio", "sbps"),
    L6 = c("strokeha", "ckidney", "type2")
  )

  for (latent_var in names(latent_connections)) {
    for (observed_var in latent_connections[[latent_var]]) {
      amat_enriched[latent_var, observed_var] <- 1
    }
  }

  dag_enriched <- empty.graph(c(observed_vars, latent_var_names))
  amat(dag_enriched) <- amat_enriched

  # --- Initialize Latent Variables ---
  # Using clustering (k-means) on the imputed, discretized sample
  df_with_latent <- formatted_sample
  for (i in 1:length(num_latent_states)) {
    current_num_states <- num_latent_states[i]

    for (latent_var in names(latent_connections)) {
      children <- latent_connections[[latent_var]]
      kmeans_data <- imputed_sample[, children, drop = FALSE]
      for(col in 1:ncol(kmeans_data)) {
        if(is.factor(kmeans_data[,col])) {
          kmeans_data[,col] <- as.numeric(kmeans_data[,col])
        }
      }

      kmeans_result <- kmeans(
        kmeans_data,
        centers = current_num_states,
        iter.max = 10,
        nstart = 20
      )
      df_with_latent[[latent_var]] <- as.factor(kmeans_result$cluster[sample(nrow(imputed_sample), nrow(df_with_latent))])
      levels(df_with_latent[[latent_var]]) <- 0:(current_num_states - 1)
    }
  }

  # --- Run Structural EM with the enriched DAG ---
  sem_with_latent <- structural.em(
    df_with_latent,
    maximize = "hc",
    maximize.args = list(score = "aic-cg"),
    fit = "mle",
    fit.args = list(replace.unidentifiable = TRUE),
    return.all = TRUE,
    start = dag_enriched,
    max.iter = 5,
    debug = FALSE
  )

  # --- Generate and Post-process Synthetic Data ---
  syn_sample <- rbn(sem_with_latent$fitted, n = nrow(ground_truth_sample))

  syn_sample[, 'bmi'] <- round(syn_sample[, 'bmi'], 1)
  syn_sample[, 'choleratio'] <- round(syn_sample[, 'choleratio'], 1)
  syn_sample[, 'sbp'] <- round(syn_sample[, 'sbp'], 0)
  syn_sample[, 'sbps'] <- abs(round(syn_sample[, 'sbps'], 2))
  syn_sample[, 'age'] <- round(syn_sample[, 'age'], 0)

  syn_sample <- syn_sample[!(syn_sample$gender == "F" & syn_sample$impot == 1), ]

  return(syn_sample)
}

# -----------------------------------------------------------------------------
# Confidence Analysis of Edges (Helper Functions)
# -----------------------------------------------------------------------------

# This part of the code is for analysis and is kept separate from the main generation function.

run_rfci_bootstrap <- function(formatted_sample) {
  discretization_bins <- list(bmi = 6, age = 4, choleratio = 5, sbp = 5, sbps = 5)
  discretized_sample <- discretize_data(formatted_sample, discretization_bins)
  imputed_sample <- missForest(discretized_sample)$ximp
  num_levels <- count_levels(imputed_sample)
  integer_sample <- to_integer_values(imputed_sample)

  bootstrap_samples <- list()
  for (i in 1:num_bootstrap_samples) {
    index <- sample(1:nrow(integer_sample), size = nrow(integer_sample), replace = TRUE)
    bootstrap_samples[[i]] <- integer_sample[index, ]
  }

  rfci_results <- list()
  for (i in 1:num_bootstrap_samples) {
    suffStat <- list(dm = bootstrap_samples[[i]], nlev = num_levels, adaptDF = FALSE)
    rfci_results[[i]] <- rfci(
      suffStat,
      indepTest = disCItest,
      alpha = 0.9,
      skel.method = 'stable',
      labels = colnames(discretized_sample),
      verbose = TRUE,
      m.max = 3
    )
  }
  return(rfci_results)
}

find_edge_confidence <- function(rfci_results, var_names) {
  num_models <- length(rfci_results)
  confidence_matrix <- matrix(0, nrow = length(var_names), ncol = length(var_names), dimnames = list(var_names, var_names))

  for (k in 1:num_models) {
    amat <- rfci_results[[k]]@amat
    for (i in 1:length(var_names)) {
      for (j in 1:length(var_names)) {
        if (i == j) next
        if (amat[i, j] == 2 && amat[j, i] == 2) { # Bidirected edge
          confidence_matrix[i, j] <- confidence_matrix[i, j] + 1
        } else if ((amat[i, j] == 2 && amat[j, i] == 1) || (amat[i, j] == 1 && amat[j, i] == 2)) { # o--> edge
          confidence_matrix[i, j] <- confidence_matrix[i, j] + 0.5
        } else if (amat[i, j] == 1 && amat[j, i] == 1) { # o--o edge
          confidence_matrix[i, j] <- confidence_matrix[i, j] + 0.33
        }
      }
    }
  }

  confidence_matrix <- confidence_matrix / num_models

  confidence_df <- as.data.frame(as.table(confidence_matrix))
  names(confidence_df) <- c("variable1", "variable2", "confidenceLevel")
  return(confidence_df)
}

plot_edge_confidence <- function(confidence_df, formatted_sample) {
  variables <- unique(confidence_df$variable1)
  plot_list <- list()

  for (i in 1:length(variables)) {
    var_name <- as.character(variables[i])
    subset_df <- subset(confidence_df, variable1 == var_name)

    threshold <- ifelse(is.numeric(formatted_sample[[var_name]]), 0.7, 0.9)

    p <- ggplot(subset_df, aes(x = variable2, y = confidenceLevel)) +
      geom_bar(stat = "identity", fill = "chartreuse3", colour = "black") +
      ggtitle(paste("Edge Confidence for:", var_name)) +
      ylim(0, 1) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    plot_list[[i]] <- p
  }

  grid.arrange(grobs = plot_list, ncol = 3)
}


# -----------------------------------------------------------------------------
# Main Execution Block
# -----------------------------------------------------------------------------

# This block will only run when the script is executed directly
if (sys.nframe() == 0) {

  # Load ground truth data
  ground_truth_data <- read.delim(ground_truth_file)

  # Main loop to generate multiple synthetic datasets
  for (n in 1:num_iterations) {

    cat(paste("\n--- Starting Iteration", n, "of", num_iterations, "---\n"))

    # Take a new sample for the iteration
    gt_sample <- ground_truth_data[sample(nrow(ground_truth_data), 100000), ]

    # Write the ground truth sample for this iteration
    gt_sample_path <- file.path(results_dir, paste0(n, gt_sample_prefix, ".txt"))
    write.table(gt_sample, gt_sample_path, sep = "\t", row.names = FALSE)

    # Generate the synthetic sample
    syn_sample <- generate_synthetic_data(gt_sample)

    # Write the synthetic sample for this iteration
    syn_sample_path <- file.path(results_dir, paste0(n, syn_sample_prefix, ".txt"))
    write.table(syn_sample, syn_sample_path, sep = "\t", row.names = FALSE)

    cat(paste("\n--- Iteration", n, "Complete ---\n"))
  }

  # --- Perform Confidence Analysis ---
  cat("\n--- Performing Edge Confidence Analysis ---\n")
  initial_sample <- ground_truth_data[sample(nrow(ground_truth_data), sample_size), ]
  formatted_sample_for_analysis <- format_variables(initial_sample)
  rfci_results <- run_rfci_bootstrap(formatted_sample_for_analysis)
  discretized_sample_for_analysis <- discretize_data(formatted_sample_for_analysis, list(bmi = 6, age = 4, choleratio = 5, sbp = 5, sbps = 5))
  edge_confidence <- find_edge_confidence(rfci_results, colnames(discretized_sample_for_analysis))
  print(edge_confidence)
  # plot_edge_confidence(edge_confidence, formatted_sample_for_analysis)

}
