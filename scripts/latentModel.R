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
# Data Loading and Preprocessing
# -----------------------------------------------------------------------------

# Load ground truth data
ground_truth_data <- read.delim(ground_truth_file)

# Take a sample from the ground truth data
initial_sample <- ground_truth_data[sample(nrow(ground_truth_data), sample_size), ]

# Define a function to format variables to their correct types
format_variables <- function(df) {
  # Factor variables
  factor_vars <- c("strokeha", "af", "atyantip", "steroid", "impot", "migr",
                   "ra", "ckidney", "semi", "sle", "treathyp", "type1",
                   "type2", "ethr", "smoking", "fh_cad", "gender", "region")
  for (var in factor_vars) {
    df[[var]] <- as.factor(df[[var]])
  }

  # Numeric variables
  numeric_vars <- c("age", "bmi", "choleratio", "sbp", "sbps")
  for (var in numeric_vars) {
    df[[var]] <- as.numeric(df[[var]])
  }

  return(df)
}

# Apply variable formatting
formatted_sample <- format_variables(initial_sample)

# Function that transforms the dataset into integer values (from 0..(n-1))
to_integer_values <- function(df) {
  for (k in 1:ncol(df)) {
    df[, k] <- as.integer(df[, k]) - 1
  }
  return(df)
}

# Function to count the number of levels in each factor variable
count_levels <- function(df) {
  num_levels <- sapply(df, nlevels)
  return(num_levels)
}

# Function for discretization of continuous variables using equal intervals
discretize_data <- function(df, discretization_bins) {
  df_discrete <- df

  # Discretize specified numeric columns
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
# Latent Variable Experiment
# -----------------------------------------------------------------------------

# Discretize continuous variables
discretization_bins <- list(bmi = 6, age = 4, choleratio = 5, sbp = 5, sbps = 5)
discretized_sample <- discretize_data(formatted_sample, discretization_bins)

# Impute missing values using Random Forest
imputed_sample <- missForest(discretized_sample)$ximp

# Get the number of levels for each variable
num_levels <- count_levels(imputed_sample)

# Convert dataframe to integer values starting from 0
integer_sample <- to_integer_values(imputed_sample)

# Bootstrap resampling to learn network structure
bootstrap_samples <- list()
for (i in 1:num_bootstrap_samples) {
  index <- sample(1:nrow(integer_sample), size = nrow(integer_sample), replace = TRUE)
  bootstrap_samples[[i]] <- integer_sample[index, ]
}

# Run the RFCI algorithm on bootstrap samples
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

# Plot the resulting PAGs (Partial Ancestral Graphs)
# for (i in 1:length(rfci_results)) {
#   plot(rfci_results[[i]])
# }

# -----------------------------------------------------------------------------
# Structural EM for Model with Latent Variables
# -----------------------------------------------------------------------------

original_sample_df <- formatted_sample

for (n in 1:num_iterations) {

  # Take a new sample for the iteration
  iter_sample <- original_sample_df[sample(nrow(original_sample_df), 100000), ]

  # Write the ground truth sample for this iteration
  gt_sample_path <- file.path(results_dir, paste0(n, gt_sample_prefix, ".txt"))
  write.table(iter_sample, gt_sample_path, sep = "\t", row.names = FALSE)

  # Learn structure without latent variables
  sem_no_latent <- structural.em(
    iter_sample,
    maximize = "hc",
    fit = "mle",
    return.all = TRUE,
    start = NULL,
    max.iter = 5
  )

  # --- Enrich the network structure with latent variables ---

  observed_vars <- names(iter_sample)
  latent_var_names <- paste("L", 1:num_latent_vars, sep = "")

  # Create an empty adjacency matrix for the enriched graph
  amat_enriched <- matrix(
    0,
    nrow = length(observed_vars) + num_latent_vars,
    ncol = length(observed_vars) + num_latent_vars,
    dimnames = list(c(observed_vars, latent_var_names), c(observed_vars, latent_var_names))
  )

  # Add edges from the original DAG
  amat_enriched[observed_vars, observed_vars] <- amat(sem_no_latent$dag)

  # Define connections from latent to observed variables
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

  # Create the enriched DAG
  dag_enriched <- empty.graph(c(observed_vars, latent_var_names))
  amat(dag_enriched) <- amat_enriched

  # --- Compare different numbers of states for the latent variables ---

  df_with_latent_clustering <- list()
  df_with_latent_random <- list()

  for (i in 1:length(num_latent_states)) {
    current_num_states <- num_latent_states[i]
    df_clust <- iter_sample
    df_rand <- iter_sample

    # Initialize latent variables using clustering (k-means)
    for (latent_var in names(latent_connections)) {
      children <- latent_connections[[latent_var]]
      kmeans_data <- imputed_sample[, children, drop = FALSE]
      # Convert factors to numeric for kmeans
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
      df_clust[[latent_var]] <- as.factor(kmeans_result$cluster[sample(nrow(imputed_sample), nrow(df_clust))])
      levels(df_clust[[latent_var]]) <- 0:(current_num_states - 1)

      # Initialize latent variables with random values
      df_rand[[latent_var]] <- as.factor(sample(0:(current_num_states - 1), size = nrow(df_rand), replace = TRUE))
    }

    df_with_latent_clustering[[i]] <- df_clust
    df_with_latent_random[[i]] <- df_rand
  }

  # --- Run Structural EM with the enriched DAG ---

  sem_results_random <- list()
  synthetic_data_random <- list()

  for (j in 1:length(df_with_latent_clustering)) {
    sem_results_random[[j]] <- structural.em(
      df_with_latent_clustering[[j]],
      maximize = "hc",
      maximize.args = list(score = "aic-cg"),
      fit = "mle",
      fit.args = list(replace.unidentifiable = TRUE),
      return.all = TRUE,
      start = dag_enriched,
      max.iter = 5,
      debug = FALSE
    )

    # Generate synthetic data
    syn_sample <- rbn(sem_results_random[[j]]$fitted, n = nrow(iter_sample))

    # Post-process synthetic data
    syn_sample[, 'bmi'] <- round(syn_sample[, 'bmi'], 1)
    syn_sample[, 'choleratio'] <- round(syn_sample[, 'choleratio'], 1)
    syn_sample[, 'sbp'] <- round(syn_sample[, 'sbp'], 0)
    syn_sample[, 'sbps'] <- abs(round(syn_sample[, 'sbps'], 2))
    syn_sample[, 'age'] <- round(syn_sample[, 'age'], 0)

    # Apply biological checks (e.g., remove impotent females)
    syn_sample <- syn_sample[!(syn_sample$gender == "F" & syn_sample$impot == 1), ]

    # Write the synthetic sample for this iteration
    syn_sample_path <- file.path(results_dir, paste0(n, syn_sample_prefix, ".txt"))
    write.table(syn_sample, syn_sample_path, sep = "\t", row.names = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Confidence Analysis of Edges
# -----------------------------------------------------------------------------

# Function to calculate the confidence of edges between variables from RFCI results
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

  # Convert matrix to a data frame for easier plotting
  confidence_df <- as.data.frame(as.table(confidence_matrix))
  names(confidence_df) <- c("variable1", "variable2", "confidenceLevel")
  return(confidence_df)
}

# Calculate edge confidence
edge_confidence <- find_edge_confidence(rfci_results, colnames(discretized_sample))
print(edge_confidence)

# Function to plot the confidence levels for each variable
plot_edge_confidence <- function(confidence_df) {
  variables <- unique(confidence_df$variable1)
  plot_list <- list()

  for (i in 1:length(variables)) {
    var_name <- variables[i]
    subset_df <- subset(confidence_df, variable1 == var_name)

    # Define a threshold for highlighting confident edges
    threshold <- ifelse(is.numeric(formatted_sample[[var_name]]), 0.7, 0.9)

    p <- ggplot(subset_df, aes(x = variable2, y = confidenceLevel)) +
      geom_bar(stat = "identity", fill = "chartreuse3", colour = "black") +
      ggtitle(paste("Edge Confidence for:", var_name)) +
      ylim(0, 1) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    plot_list[[i]] <- p
  }

  # Arrange and print plots
  # This part may need adjustment depending on the number of variables
  grid.arrange(grobs = plot_list, ncol = 3)
}

# Plot the results
# plot_edge_confidence(edge_confidence)
