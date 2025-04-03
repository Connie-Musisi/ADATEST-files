# Libraries
library(DESeq2)
library(ANCOMBC)
library(MicrobiomeStat)
library(ALDEx2)
library(signtrans)
library(RioNorm2)
library(corncob)
library(GUniFrac)
library(phyloseq)
library(camp)
library(LDM)

## Statistical Tests to compare results
# Input: Dataset in phyloseq object, with grouping variable coded as 'group'
# Output: Data frame of results: test statistic (if provided), raw p-value, adjusted p-value 

# Method: camp
camp_method = function(data){
  OTU = as.data.frame(otu_table(data))
  group = sample_data(data)$group
  camp_res = camp::camp(as.matrix(OTU),group)
  pVals = camp_res
  adjPVals <- p.adjust(pVals, method = "BH")
  
  out = data.frame("pval" = pVals, "adjP" = adjPVals)
  rownames(out) = rownames(tax_table(data))
  return(out)
}

# Method: ADAPT
ADAPT_method = function(data){
  adapt_res = ADAPT::adapt(data, cond.var = "group")
  pVals = adapt_res@details[["pval"]]
  adjPVals <- adapt_res@details[["adjusted_pval"]]
  adjPVals_BH <- p.adjust(pVals, method = "BH")
  
  teststat = adapt_res@details[["teststat"]]
  out = data.frame("pval" = pVals, "adjP" = adjPVals, "adjP_BH" = adjPVals_BH, 
                   "teststat"=teststat)
  rownames(out) = adapt_res@details[["Taxa"]]
  return(out)
}

# Method: DESEQ2
deseq2_method = function(data){

  if (!taxa_are_rows(data))
  {
    data <- t(data)
  } else {
    data <- data
  }
  
  data@sam_data$group =  factor(data@sam_data$group) 
  
  dds <- phyloseq::phyloseq_to_deseq2(data, design = as.formula("~ group"))
  
  dds<-DESeq2::estimateSizeFactors(dds, type = 'poscounts')
  
  test_deseq2 = DESeq2::DESeq(dds)
  
  res = DESeq2::results(test_deseq2,independentFiltering=FALSE)
  
  res_pval = res$pvalue
  
  res_adjP = p.adjust(as.numeric(res$pvalue),method="BH")
  
  res_stat = res$stat
  
  out = data.frame(stat = res_stat, pval = res_pval, adjP = res_adjP) 
  rownames(out) = res@rownames
    
  return(out)
}    
    
# Method: ALDex2
aldex_method = function(data){
  
  ## Adjust data &format for INPUT 
  # Counts should be in data frame/ matrix format
  if(taxa_are_rows(data)){
    counts = data.frame(otu_table(data))
  }else{
    counts = data.frame(t(otu_table(data)))
  }
  
  # Vector with grouping labels
  group = sample_data(data)$group
  
  ## Run ALDEx2
  res = ALDEx2::aldex(counts,group, test='t')
  
  # Effect size
  res_stat = res$effect
  rownames(res) -> names(res_stat)
  
  
  # raw p-values  
  res_pvalue_t = res$we.ep   # t-test 
  rownames(res) -> names(res_pvalue_t)
  
  res_pvalue_w = res$wi.ep   # wilcoxon
  rownames(res) -> names(res_pvalue_w)
  
  res_adjP_t = p.adjust(res_pvalue_t,method="BH")
  res_adjP_w = p.adjust(res_pvalue_w,method="BH")
  
  
  # grouping results in data frames
  out_t = data.frame( stat=res_stat,pval=res_pvalue_t, adjP = res_adjP_t)
  out_w  = data.frame( stat=res_stat,pval=res_pvalue_w, adjP = res_adjP_w )
  
  out = list(ALDEx2_t = out_t,
             ALDEx2_w = out_w,
             total = res)
  return(out)
}

# Method: Corncob
corncob_method = function(data){
  if(!taxa_are_rows(data)){
data <- t(data)
  }else{}
  
  # Filter out taxa with zero total counts
  otuMat <- as(otu_table(data), "matrix")
  validTaxa <- rowSums(otuMat) > 0
  data <- prune_taxa(validTaxa, data)
  
  otu_table_data <- as(otu_table(data), "matrix")
  variances <- apply(otu_table_data, 1, var)
  non_zero_variance <- variances > 0
  
  # Filter out taxa with zero variance
  data2 <- prune_taxa(non_zero_variance, data)
  # Run Corncob
  res = corncob::differentialTest(formula = ~ group,
                                  phi.formula = ~ group,
                                  formula_null = ~ 1,
                                  phi.formula_null = ~ group,
                                  test = "LRT", boot = FALSE,
                                  data = data2,
                                  #fdr = "BH",
                                  fdr_cutoff = 0.05)
  
  # Save results
  #res_stat = c()
  #for (i in (1:ntaxa(data))){
  #  res_stat = c(res_stat,
  #               
  #               if (length(res[["all_models"]][[i]])!=1){
  #                 res[["all_models"]][[i]][["coefficients"]][2, 3]
  #               }else("NA")
  #               
  #  )
  #}
  res_pval = res[["p"]]
  #res_adjP = p.adjust(as.numeric(res_pval),method="BH")
  #res_adjP = (qvalue(res_pval))$qvalues
  res_adjP = res[["p_fdr"]]
  names(res_pval) -> names(res_adjP)
  #names(res_pval) -> names(res_stat)
  #out = data.frame(stat = res_stat, pval = res_pval, adjP = res_adjP) 
  out = data.frame(pval = res_pval, adjP = res_adjP) 
  
  return(out)
}

# Method: ANCOM-BC2
ancombc2_method <- function(data){
  # Ensure that taxa are rows by transposing the OTU table if needed
  if (!taxa_are_rows(data)){
    otu_table(data) <- t(otu_table(data))
  }
  
  # Convert group variable to factor
  sample_data(data)$group <- as.factor(sample_data(data)$group)
  
  # Remove taxa with zero variance
  otu_table_data <- as(otu_table(data), "matrix")
  variances <- apply(otu_table_data, 1, var)
  non_zero_variance <- variances > 0
  
  # Filter out taxa with zero variance
  data2 <- prune_taxa(non_zero_variance, data)
  
  # Since we now have taxa as rows, set taxa_are_rows = TRUE
  res <- ancombc2(data2, taxa_are_rows = TRUE, 
                  fix_formula = "group", group = "group",
                  prv_cut = 0.20)
  
  pVals <- res[["res"]][["p_group1"]]
  adjPVals <- res[["res"]][["q_group1"]]
  adjPVals_BH <- p.adjust(pVals, method = "BH")
  
  out <- data.frame("pval" = pVals, "adjP" = adjPVals, "adjP_BH" = adjPVals_BH)
  rownames(out) <- res[["res"]][["taxon"]]
  return(out)
}

# Method: LinDA (microbiomeStat package)
LinDA_test = function(data){
  if(taxa_are_rows(data)){
    otu_matrix <- as.data.frame(otu_table(data))
  }else{
otu_matrix <- as.data.frame(t(otu_table(data)))
}
# Extract grouping variable
group <- data.frame(group=sample_data(data)$group)


linda.obj <- MicrobiomeStat::linda(otu_matrix, group, formula = '~group', alpha = 0.05)
pVals = linda.obj[["output"]][["group"]][["pvalue"]]
#adjPVals <- p.adjust(pVals, method = "BH")
adjPVals = linda.obj[["output"]][["group"]][["padj"]]


stat = linda.obj[["output"]][["group"]][["stat"]]
out = data.frame("pval" = pVals, "adjP" = adjPVals, "stat"=stat)
rownames(out) = rownames(linda.obj[["feature.dat.use"]])
return(out)
}

# Method: Wilcoxon rank sum test : raw results - tss & clr normalization
WilcoxTest_TSS = function(data){
  if(!taxa_are_rows(data)){
    data = t(data)
  }else{}
  data_norm = transform_sample_counts(data, function(x) x / sum(x))
  counts <- as(otu_table(data_norm), "matrix")
  groupVar = 'group'
  groupBinary <- get_variable(data_norm, groupVar)
  groupBinary <- groupBinary == groupBinary[1L]
  stat <- apply(counts, MARGIN = 1L, 
                function(x, ind) wilcox.test(x = x[ind], y = x[!ind], exact = FALSE)$statistic, ind = groupBinary)
  
  pVals <- apply(counts, MARGIN = 1L, 
                 function(x, ind) wilcox.test(x = x[ind], y = x[!ind], exact = FALSE)$p.value, ind = groupBinary)
  
  
  adjPVals <- p.adjust(pVals, method = "BH")
  
  out = data.frame("pval" = pVals, "adjP" = adjPVals, "stat"=stat)
  return(out)
}

WilcoxTest_clr = function(data){
  if(!taxa_are_rows(data)){
    data = t(data)
  }else{
    data=data
  }
  data_norm = microbiome::transform(data, transform = "clr")
  counts <- as(otu_table(data_norm), "matrix")
  groupVar = 'group'
  groupBinary <- get_variable(data_norm, groupVar)
  groupBinary <- groupBinary == groupBinary[1L]
  stat <- apply(counts, MARGIN = 1L, 
                function(x, ind) wilcox.test(x = x[ind], y = x[!ind], exact = FALSE)$statistic, ind = groupBinary)
  
  pVals <- apply(counts, MARGIN = 1L, 
                 function(x, ind) wilcox.test(x = x[ind], y = x[!ind], exact = FALSE)$p.value, ind = groupBinary)
  
  
  adjPVals <- p.adjust(pVals, method = "BH")
  
  out = data.frame("pval" = pVals, "adjP" = adjPVals, "stat"=stat)
  return(out)
}

# Method: ZicoSeq
ZicoSeq_method =function(data){
  if(taxa_are_rows(data)){
    otu_matrix <- as.matrix(otu_table(data))
  }else{
    otu_matrix <- as.matrix(t(otu_table(data)))
  }
  feature_data <- otu_matrix[rowSums(otu_matrix) > 0, ]
  meta_data <- data.frame(sample_data(data))
  res = ZicoSeq(meta.dat = meta_data, feature.dat = feature_data,
                grp.name = "group", feature.dat.type = "count",
                prev.filter = 0.2, mean.abund.filter = 0,
                max.abund.filter = 0.002, min.prop = 0)
  res_pval = res$p.raw
  
  res_adjP2 =res$p.adj.fdr
  
  res_adjP = p.adjust(as.numeric(res$p.raw),method="BH")
  
  res_stat = res$F0
  
  out = data.frame(stat = res_stat, pval = res_pval, adjP = res_adjP, adjP2 = res_adjP2) 
  
  #rownames(out) = rownames(res[["feature.dat"]])
  return(out)
}


## Evaluate results based on thresholds 
# Input: adjusted p-values, vector with index of I0 and I1 OTUs
# Output: False positive proportion, True positive proportion, False discovery proportion
eval<-function(alpha,n.taxa,pval,de.ind){
  results<-c()
  for (j in (1:n.taxa)){
    if(is.na(pval[j])){results[j]<-'NA'
    next}   
    if (pval[j]<=alpha & isDA[j]=='TRUE'){results[j]='TP'}
    if (pval[j]<=alpha & isDA[j]=='FALSE'){results[j]='FP'}
    if (pval[j]>alpha &  isDA[j]=='TRUE'){results[j]='FN'}
    if (pval[j]>alpha &  isDA[j]=='FALSE'){results[j]='TN'}
  }
  TP<-sum(results=='TP')
  FP<-sum(results=='FP')
  TN<-sum(results=='TN')
  FN<-sum(results=='FN')
  
  sensitivity<-TP/(TP + FN) # recall, True positive rate
  specificity<-TN/(TN + FP)
  FDR<- if ((TP + FP) == 0) 
    0 else FP/(TP + FP)
  type1error <- FP/(TN+FP)
  precision<- if ((TP + FP) == 0) 
    0 else TP/(TP + FP)
  FPR <- FP/(TN+FP)
  F1 <- 2/((1/sensitivity)+(1/precision))
  
  
  return(cbind(sensitivity,FDR,type1error))
}
