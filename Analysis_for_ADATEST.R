# -------------------------------
# Simulation run: process all settings
# -------------------------------
# Load data
path = "~/Sim_datasets/Dietswap_list"
files <- list.files(path = path, pattern = "\\.RData$", full.names = TRUE)
for(file in files) {load(file)}

# List of simulated dataset groups by setting
datasets_per_setting <- list(
  Setting1_n50_b1_rnt20 = Setting1_n50_b1_rnt20,
  Setting2_n200_b1_rnt20 = Setting2_n200_b1_rnt20,
  Setting3_n50_b1.5_rnt20 = Setting3_n50_b1.5_rnt20,
  Setting4_n200_b1.5_rnt20 = Setting4_n200_b1.5_rnt20,
  Setting5_n50_b3_rnt20 = Setting5_n50_b3_rnt20,
  Setting6_n200_b3_rnt20 = Setting6_n200_b3_rnt20,
  Setting7_n50_b5_rnt20 = Setting7_n50_b5_rnt20,
  Setting8_n200_b5_rnt20 = Setting8_n200_b5_rnt20,
  Setting9_n50_b1_rnt50 = Setting9_n50_b1_rnt50,
  Setting10_n200_b1_rnt50 = Setting10_n200_b1_rnt50,
  Setting11_n50_b1.5_rnt50 = Setting11_n50_b1.5_rnt50,
  Setting12_n200_b1.5_rnt50 = Setting12_n200_b1.5_rnt50,
  Setting13_n50_b3_rnt50 = Setting13_n50_b3_rnt50,
  Setting14_n200_b3_rnt50 = Setting14_n200_b3_rnt50,
  Setting15_n50_b5_rnt50 = Setting15_n50_b5_rnt50,
  Setting16_n200_b5_rnt50 = Setting16_n200_b5_rnt50
)

# Vector of setting names (keys in datasets_per_setting)
settings_list <- c("Setting1_n50_b1_rnt20", "Setting2_n200_b1_rnt20",
                   "Setting3_n50_b1.5_rnt20", "Setting4_n200_b1.5_rnt20",
                   "Setting5_n50_b3_rnt20", "Setting6_n200_b3_rnt20",
                   "Setting7_n50_b5_rnt20", "Setting8_n200_b5_rnt20",
                   "Setting9_n50_b1_rnt50", "Setting10_n200_b1_rnt50",
                   "Setting11_n50_b1.5_rnt50", "Setting12_n200_b1.5_rnt50",
                   "Setting13_n50_b3_rnt50", "Setting14_n200_b3_rnt50",
                   "Setting15_n50_b5_rnt50", "Setting16_n200_b5_rnt50"
                 )


# Function to extract setting info from the name
extract_setting_info <- function(setting_name) {
  parts <- strsplit(setting_name, "_")[[1]]
  setting_number <- gsub("Setting", "", parts[1])
  sample_size <- as.numeric(gsub("n", "", parts[2]))
  beta <- as.numeric(gsub("b", "", parts[3]))
  causal_taxa <- as.numeric(gsub("rnt", "", parts[4]))
  
  list(Setting = paste("Setting", setting_number),
       SampleSize = sample_size,
       Beta = beta,
       CausalTaxa = causal_taxa)
}

# Initialize list to store full simulation outputs (sim_res) for every dataset
sim_all <- list()
final_summary <- list()

#set.seed(19)

# Loop over each setting
for (setting_name in settings_list) {
  
  datasets <- datasets_per_setting[[setting_name]]
  if (!is.null(names(datasets))) {
    datasets <- datasets[order(as.numeric(names(datasets)))]
  }
  
  setting_info <- extract_setting_info(setting_name)
  
  sim_all[[setting_name]] <- list()
  eval_list <- list()
  
  # Loop over each dataset in the current setting
  for (j in seq_along(datasets)) {
    sim_res <- ADATEST(datasets[[j]])
    sim_all[[setting_name]][[j]] <- sim_res
    eval_df <- do.call(cbind, sim_res$Results)
    eval_list[[j]] <- eval_df
  }
  
  # Combine results for the 100 simulations and Calculate average sensitivity, FDR, and type1error
  combined_results <- do.call(rbind, eval_list)
  
  avg_metrics <- colMeans(combined_results[, c("FDR", "TPR", "Type_1_error")], na.rm = TRUE)
  
  # Create a summary row for the current setting
  summary_row <- data.frame(
    SimMethod = "MIDASim",
    Method = "ADATEST",
    FDR = avg_metrics["FDR"],
    Sensitivity = avg_metrics["TPR"],
    Type1Error = avg_metrics["Type_1_error"],
    SampleSize = setting_info$SampleSize,
    Beta = setting_info$Beta,
    CausalTaxa = setting_info$CausalTaxa,
    Setting = setting_info$Setting,
    stringsAsFactors = FALSE
  )
  
  final_summary[[length(final_summary) + 1]] <- summary_row
}

final_summary_df <- do.call(rbind, final_summary)

View(final_summary_df)
ADATEST_dietswap_res = sim_all
ADATEST_dietswap_eval_res = final_summary_df
