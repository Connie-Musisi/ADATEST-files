# Load data
path = "~/Sim_datasets/Dietswap_list"
files <- list.files(path = path, pattern = "\\.RData$", full.names = TRUE)
for(file in files) {load(file)}


# Results for all the methods fro all settings
methods_list <- list(
 "ADAPT" = ADAPT_method,
 "WilcoxCLR" = WilcoxTest_clr,
 "WilcoxTSS" = WilcoxTest_TSS,
 "ANCOMBC2" = ancombc2_method,
 "DESeq2" = deseq2_method,
 "ALDEx2" = aldex_method,
 "LinDA" = LinDA_test,
 "ZicoSeq" = ZicoSeq_method,
 "Corncob" = corncob_method,
 "camp" = camp_method
)


settings_list <- c(#"Setting16_highcomp" 
"Setting1_n50_b1_rnt20", "Setting2_n200_b1_rnt20",
"Setting3_n50_b1.5_rnt20", "Setting4_n200_b1.5_rnt20",
"Setting5_n50_b3_rnt20", "Setting6_n200_b3_rnt20",
"Setting7_n50_b5_rnt20", "Setting8_n200_b5_rnt20",
"Setting9_n50_b1_rnt50", "Setting10_n200_b1_rnt50",
"Setting11_n50_b1.5_rnt50", "Setting12_n200_b1.5_rnt50",
"Setting13_n50_b3_rnt50", "Setting14_n200_b3_rnt50",
"Setting15_n50_b5_rnt50", "Setting16_n200_b5_rnt50"
)

extract_setting_info <- function(setting_name) {
  parts <- strsplit(setting_name, "_")[[1]]
  setting_number <- gsub("Setting", "", parts[1])
  sample_size <- gsub("n", "", parts[2])
  beta <- gsub("b", "", parts[3])
  causal_taxa <- gsub("rnt", "", parts[4])
  
  list(Setting = paste("Setting", setting_number),
       SampleSize = as.numeric(sample_size),
       Beta = as.numeric(beta),
       CausalTaxa = as.numeric(causal_taxa))
}

datasets_per_setting <- list(
  #Setting16_highcomp=Setting16_highcomp
  Setting1_n50_b1_rnt20=Setting1_n50_b1_rnt20,
  Setting2_n200_b1_rnt20=Setting2_n200_b1_rnt20,
  Setting3_n50_b1.5_rnt20=Setting3_n50_b1.5_rnt20,
  Setting4_n200_b1.5_rnt20=Setting4_n200_b1.5_rnt20,
  Setting5_n50_b3_rnt20=Setting5_n50_b3_rnt20,
  Setting6_n200_b3_rnt20=Setting6_n200_b3_rnt20,
  Setting7_n50_b5_rnt20=Setting7_n50_b5_rnt20,
  Setting8_n200_b5_rnt20=Setting8_n200_b5_rnt20,
  Setting9_n50_b1_rnt50=Setting9_n50_b1_rnt50,
  Setting10_n200_b1_rnt50=Setting10_n200_b1_rnt50,
  Setting11_n50_b1.5_rnt50=Setting11_n50_b1.5_rnt50,
  Setting12_n200_b1.5_rnt50=Setting12_n200_b1.5_rnt50,
  Setting13_n50_b3_rnt50=Setting13_n50_b3_rnt50,
  Setting14_n200_b3_rnt50=Setting14_n200_b3_rnt50,
  Setting15_n50_b5_rnt50=Setting15_n50_b5_rnt50,
  Setting16_n200_b5_rnt50=Setting16_n200_b5_rnt50
)


results_list <- list()
for (setting_name in settings_list) {
    datasets <- datasets_per_setting[[setting_name]]
  
  # Extract setting information (Setting X, sample size, beta, causal taxa)
  setting_info <- extract_setting_info(setting_name)

  for (method_name in names(methods_list)) {
        method_function <- methods_list[[method_name]]
        method_results <- list()
    
    for (dataset in datasets) {
      
      res <- method_function(dataset)
      
      if (method_name == "ALDEx2") {
        res <- res[["ALDEx2_w"]]
      }
      
      tax_table_data <- tax_table(dataset)
      isDA <- tax_table_data[rownames(res), "isDA"]
      isDA <- as.vector(isDA)
      
      result <- eval(alpha=0.05, pval=res$adjP, n.taxa=length(isDA), de.ind=isDA)
      method_results[[length(method_results) + 1]] <- result
    }
    
    combined_results <- do.call(rbind, method_results)
    avg_result <- colMeans(combined_results[, c("sensitivity", "FDR", "type1error")])
    
    results_list[[length(results_list) + 1]] <- data.frame(
      SimMethod = "MIDASim",
      Method = method_name,
      FDR = avg_result["FDR"],
      Sensitivity = avg_result["sensitivity"],
      Type1Error = avg_result["type1error"],
      SampleSize = setting_info$SampleSize,
      Beta = setting_info$Beta,
      CausalTaxa = setting_info$CausalTaxa,
      Setting = setting_info$Setting
    )
  }
}

final_results <- do.call(rbind, results_list)
View(final_results)

competitor_highcomp_results = final_results
