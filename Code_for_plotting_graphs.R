library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(dplyr)

##############################################################
# 1. Result plots : Distance from largest sensitivity
##############################################################

# Load data
allresults = read.csv("all_countibd_results.csv", sep=";")
allres <- allresults %>% filter(Beta != 1) # no effect
allres <- allres %>% filter(Setting != "High comp")

# Results from one sample size group
all_res<-allres[allres$SampleSize==200,]
all_res$Sensitivity <- as.numeric(all_res$Sensitivity)
all_res$FDR <- as.numeric(all_res$FDR)
all_res$scen.set<-paste(all_res$Setting, all_res$Scenario,
                        all_res$samp, all_res$fc)

Results.Summary<-function(all_res, method="Adaptive") {
  ressum<-data.frame(
    scenario=unique(all_res$scen.set),
    rank.sens=NA,  # rank of sensitivity among competitors
    sens.diff=NA,  # difference with largest sens
    ctrl.05=NA,    # in FDR control (<=0.05)
    ctrl.10=NA,    # FDR<=0.10
    sens.diff.05=NA,  # difference with largest sens (among FDR<0.05)
    sens.diff.10=NA   # difference with largest sens (among FDR<0.10)
  )
  
  cnt<-1
  for(s in ressum$scenario) {
    tmp<-all_res[all_res$scen.set==s,]
    c05<-tmp$FDR<=0.05
    c10<-tmp$FDR<=0.10
    
    ressum$rank.sens[cnt]<-rank(-tmp$Sensitivity)[tmp$Method==method]  
    ressum$sens.diff[cnt]<-tmp$Sensitivity[tmp$Method==method]-max(tmp$Sensitivity)
    ressum$ctrl.05[cnt]<-tmp$FDR[tmp$Method==method]<=0.05
    ressum$ctrl.10[cnt]<-tmp$FDR[tmp$Method==method]<=0.10
    
    ressum$sens.diff.05[cnt]<-tmp$Sensitivity[tmp$Method==method]-max(tmp$Sensitivity[c05])
    ressum$sens.diff.05[cnt]<-ifelse(
      ressum$ctrl.05[cnt],ressum$sens.diff.05[cnt],NA
    )
    ressum$sens.diff.10[cnt]<-tmp$Sensitivity[tmp$Method==method]-max(tmp$Sensitivity[c10])
    ressum$sens.diff.10[cnt]<-ifelse(
      ressum$ctrl.10[cnt],ressum$sens.diff.10[cnt],NA
    )
    
    cnt<-cnt+1
  }
  
  return(ressum)
}

Results.Summary(all_res = all_res, method = "ADATEST")

AllResults<-list()
for(m in unique(all_res$Method)) {
  AllResults[[m]]<-Results.Summary(all_res = all_res, method = m)
}

# Sensitivity difference
SensDiff<-data.frame(Method=unique(all_res$Method),
                     se.df.max=NA,
                     se.df.Q2=NA,
                     se.df.min=NA,
                     se.df.max.05=NA,
                     se.df.Q2.05=NA,
                     se.df.min.05=NA,
                     se.df.max.10=NA,
                     se.df.Q2.10=NA,
                     se.df.min.10=NA,
                     ctrl.05=NA,
                     ctrl.10=NA)

cnt<-1
for(m in unique(all_res$Method)) {
  res<-AllResults[[m]]
  q<-quantile(-res$sens.diff)
  SensDiff$se.df.max[cnt]<-q[5]
  SensDiff$se.df.Q2[cnt]<-q[3]
  SensDiff$se.df.min[cnt]<-q[1]
  
  q<-quantile(-res$sens.diff[res$ctrl.05])
  SensDiff$se.df.max.05[cnt]<-q[5]
  SensDiff$se.df.Q2.05[cnt]<-q[3]
  SensDiff$se.df.min.05[cnt]<-q[1]
  
  q<-quantile(-res$sens.diff[res$ctrl.10])
  SensDiff$se.df.max.10[cnt]<-q[5]
  SensDiff$se.df.Q2.10[cnt]<-q[3]
  SensDiff$se.df.min.10[cnt]<-q[1]
  
  SensDiff$ctrl.05[cnt]<-round(mean(res$ctrl.05),2)
  SensDiff$ctrl.10[cnt]<-round(mean(res$ctrl.10),2)
  
  cnt<-cnt+1
}

SensDiff[,-1]<-round(SensDiff[,-1],3)
pd<-SensDiff

# Order the dataframe pd by se.df.max.10 (descending)
pd <- pd[order(-pd$se.df.max.10), ][]
pd$Method.nr <- 1:nrow(pd)

par(mar = c(5, 10, 4, 2))

# Identify the best competitor (smallest sensitivity difference)
non_adaptive_methods <- pd$Method[pd$Method != "ADATEST"]
best_method <- non_adaptive_methods[which.min(pd$se.df.max.10[pd$Method != "ADATEST"])]

# Plot the results
plot(pd$se.df.max.10, pd$Method.nr,
     xlim = c(0, 1), ylim = range(pd$Method.nr),
     ylab = "", xlab = "Difference for largest sensitivity",
     #main = "Results for n=25 per group (FDR<10%)",
     yaxt = "n", cex.lab = 1.5, cex.axis = 1.0)

axis(2, at = pd$Method.nr, labels = pd$Method, las = 2)

for (i in 1:nrow(pd)) {
  method_name <- pd$Method[i]
    line_color <- ifelse(method_name == "ADATEST", "red",
                       ifelse(method_name == best_method, "blue", "black"))
  
  text_color <- line_color
    lines(c(pd$se.df.Q2.10[i], pd$se.df.max.10[i]), c(i, i), col = line_color, lwd = 2)
    text(pd$se.df.max.10[i], i - 0.2, pd$ctrl.10[i], col = text_color, cex = 0.8, adj = 0) # Shifted down by 0.2
}





##############################################################
# 2. Scores distributions
##############################################################
# Set the path to your folder containing the .RData files

# a) Using MIDASim data results
load("~/Results/adaptive_res_set15.RData")
parest_list <- list()
for (i in 1:length(adaptive_res_set15[["Setting15_n50_b5_rnt50"]])) {
  parest <- adaptive_res_set15[["Setting15_n50_b5_rnt50"]][[i]][["Train_results"]][["Train_parest"]]
  parest_list[[i]] <- parest
}
parest <- do.call(rbind, parest_list)


# b) Using SPSimSeq/NB data results
setwd("~/Results/NB_3.3A_qval")
files <- list.files(pattern = "^Results_/d+/.RData$")
parest_list <- list()
for (i in seq_along(files)) {
  load(files[i])
  parest_list[[i]] <- as.numeric(output[["Train_results"]][["Train_parest"]])
}
parest <- do.call(rbind, parest_list)


avg_parest <- colMeans(parest)
n <- length(avg_parest)
avg_parest1  <- avg_parest[1:(n/2)]
avg_parest2 <- avg_parest[((n/2) + 1):n]

# same axis limits for the plots
ymin <- min(avg_parest1, avg_parest2)
ymax <- max(avg_parest1, avg_parest2)

par(mfrow = c(2,1),   
    oma = c(0, 4, 0, 0), 
    mar = c(4, 2, 2, 1)) 

# First group
plot(avg_parest1, xlab = "Ranks", ylab = "", cex.lab = 1.7, cex.axis = 1.2,
     ylim = c(ymin, ymax))

# Second group
plot(avg_parest2, xlab = "Ranks", ylab = "", cex.lab = 1.7, cex.axis = 1.2,
     ylim = c(ymin, ymax))

mtext("Scores", side = 2, outer = TRUE, line = 2.5, cex = 2)




##############################################################
# 3. FDR and Sensitivity plots
##############################################################
# Set the path to your folder containing the .RData files
all_results=read.csv("all_countibd_results.csv", sep = ";")
View(all_results)

all_res <- all_results %>% filter(Setting != "High comp")
all_res <- all_res %>% filter(Beta != 1)
View(all_res)

all_res <- all_res %>%
  group_by(Method, CausalTaxa, SampleSize) %>%
  arrange(Beta, .by_group = TRUE) %>%
  ungroup()

all_res$Sensitivity <- as.numeric(all_res$Sensitivity)
all_res$FDR <- as.numeric(all_res$FDR)

method_sizes <- c(
  "ADAPT" = 0.7,
  "ADATEST" = 1,  
  "ALDEx2" = 0.7,
  "ANCOMBC2" = 0.7,
  "Corncob" = 0.7,
  "DESeq2" = 0.7,
  "LinDA" = 0.7,
  "WilcoxCLR" = 0.7,
  "WilcoxTSS"=0.7,
  "ZicoSeq" = 0.7,
  "camp" = 0.7
)

# Define the method colors explicitly (including zicoseq if needed)
method_colors <- c(
  "ADAPT" = "purple", 
  "ALDEx2" = "blue", 
  "corncob" = "green4", 
  "LinDA" = "yellow", 
  "WilcoxCLR" = "darkorange", 
  "ANCOMBC2" = "brown", 
  "DESeq2" = "pink", 
  "WilcoxTSS" = "darkgrey", 
  "ADATEST" = "red", 
  "Zicoseq" = "cyan",
  "camp" = "gold" # Ensure zicoseq has a visible color
)

# FDR
fdr_plot <- ggplot(all_res, aes(x = Beta, y = FDR, color = Method, 
                                group = Method, size = Method)) +
  geom_line(lineend="round") +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "solid", color = "black") + 
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = c(1.5, 3, 5)) +
  labs(x = "Beta", y = "FDR", color = "Method") +
  scale_color_manual(values = method_colors) +
  scale_size_manual(values = method_sizes, guide = "none") +
  facet_wrap(~ CausalTaxa + SampleSize, 
             labeller = labeller(CausalTaxa = function(x) paste("CausalTaxa =", x))) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5),
    panel.border = element_blank(),
    axis.ticks = element_line(size = 0.5)
  )

print(fdr_plot)

# Sensitivity
sens_plot <- ggplot(all_res, aes(x = Beta, y = Sensitivity, color = Method,
                                 group = Method, size = Method)) +
  geom_line(lineend="round") +
  geom_point() +
  scale_x_continuous(breaks = c(1.5, 3, 5)) +
  labs(
    x = "Beta",
    y = "Sensitivity",
    color = "Method"  
  ) +
  scale_color_manual(values = method_colors) +  
  scale_size_manual(values = method_sizes, guide = "none") +
  facet_wrap(~ CausalTaxa + SampleSize, 
             labeller = labeller(CausalTaxa = function(x) paste("CausalTaxa =", x))) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "bottom",  
    panel.grid = element_blank(),       
    axis.line = element_line(size = 0.5), 
    panel.border = element_blank(),     
    axis.text.x = element_text(),       
    axis.text.y = element_text(),       
    axis.ticks = element_line(size = 0.5)
  )

print(sens_plot)




##############################################################
# 4. QQ-plots
##############################################################
# Set the path to your folder containing the .RData files
load("~/Results/adaptive_res_set15.RData")
# Load your data
original = adaptive_res_set8[["Setting8_n200_b5_rnt20"]][[1]][["Original_results"]][["org_scaled"]]
training = adaptive_res_set8[["Setting8_n200_b5_rnt20"]][[1]][["Train_results"]][["Train_scaled"]]

# Define sample groups
group1 <- 1:(nrow(original)/2) # samples in group 1
group2 <- (nrow(original)/2 + 1):(nrow(original)) # samples in group 2

# Compute LFC from original data
min.without.zero.orig <- apply(original, 1, function(x) { min(x[x != 0]) })
lfc_orig <- log2((rowMeans(original[, group2]) + min.without.zero.orig)/(rowMeans(original[, group1]) + min.without.zero.orig))

# Compute LFC for training data
min.without.zero.train <- apply(training, 1, function(x) { min(x[x != 0]) })
lfc_train <- log2((rowMeans(training[, group2]) + min.without.zero.train)/(rowMeans(training[, group1]) + min.without.zero.train))

# Split taxa into 5 clusters using quantiles
abs_lfc_orig <- abs(lfc_orig)
abs_lfc_train <- abs(lfc_train)
breaks_abs <- quantile(abs_lfc_orig, probs = seq(0, 1, length.out = 6), na.rm = TRUE)

cluster_orig <- cut(
  abs_lfc_orig,
  breaks = breaks_abs,
  labels = 1:5,
  include.lowest = TRUE
)

cluster_train <- cut(
  abs_lfc_train,
  breaks = breaks_abs,
  labels = 1:5,
  include.lowest = TRUE
)

# Make the plots for each cluster
all_abund_orig_g2 <- as.vector(original[, group2])
all_abund_train_g2 <- as.vector(training[, group2])

global_min <- min(all_abund_orig_g2, all_abund_train_g2)
global_max <- max(all_abund_orig_g2, all_abund_train_g2)

par(mfrow = c(1, 1))

for (i in 1:5) {
  orig_taxa_i  <- which(cluster_orig == i)
  train_taxa_i <- which(cluster_train == i)
    abund_orig_g1 <- as.vector(original[orig_taxa_i, group1])
  abund_train_g1 <- as.vector(training[train_taxa_i, group1])
    if (length(abund_orig_g1) == 0 || length(abund_train_g1) == 0) {
    plot(1, 1,
         type = "n",
         xlim = c(global_min, global_max),
         ylim = c(global_min, global_max),
         xlab = "Original data",
         ylab = "Training data",
         main = paste("Cluster", i, "(empty)"))
    next
  }
  
  median_abs_lfc <- median(abs_lfc_orig[orig_taxa_i])
  
  qqplot(
    abund_orig_g1, abund_train_g1,
    main = paste0("|LFC|=", round(median_abs_lfc, 2)),
    xlab = "Original data",
    ylab = "Training data",
    xlim = c(global_min, global_max),
    ylim = c(global_min, global_max),
    pch = 19, col = "blue"
  )
  abline(0, 1, col = "red", lty = 2)
}



