## Data simulation using MIDAsim

# Load packages

library(MIDASim)
library(nleqslv)
library(survival)
library(microbiome)

# General setup 
data("dietswap")
dietswap2<-subset_samples(dietswap, timepoint==2)
count.dietswap <- t(as.matrix(dietswap2@otu_table))
count.dietswap.setup <- MIDASim.setup(count.dietswap, mode = 'parametric')

# Which settings to explore

alln <- c(50, 200) # total sample size
allb <- c(1, 1.5, 3, 5)  # effect size
rnt <- c(20,50) # number of causal taxa

eg <- expand.grid(n = alln, beta = allb, rnt = rnt, 
                  stringsAsFactors = F)

for(iscen in 1:nrow(eg)){
  
  condition <- eg[iscen,]
  n.data <- condition$n
  beta1 <- condition$beta # effect of the variable of interest
  rnt <- condition$rnt
  
  n.taxa = count.dietswap.setup$n.taxa
  
  #fn <- sprintf("Setting%d(n_%d_b_%s_rnt_%d)",iscen,n.data,beta1,rnt)
  fn <- sprintf("Setting%d(highcomp)",iscen)
  
  dir.create(file.path("./Sim_datasets/Dietswap", fn), showWarnings = FALSE)
  
  
  # Top 60 of most abundant taxa are selected for possible modification
  m = 60
  top.index = order(count.dietswap.setup$mean.rel.abund, decreasing =T)[1:m]
  x.index = sample(top.index, rnt, replace =F) # get causal taxa set: randomly select rnt=10/20 taxa from the most abundant 50 taxa
  
  # Group indicator
  x.true = c(rep(0, n.data/2), rep(1, n.data/2))
  x.true = x.true - mean(x.true)
  
  # Generate 100 datasets per scenario under consideration
  
  data.all = ra.all = matrix(NA, n.data, n.taxa)
  
  # effect1 <- (x.true * beta1)
  # new.mean.rel.abund = old.mean.rel.abund= count.dietswap.setup$mean.rel.abund
  # new.mean.rel.abund[x.index[1:(rnt/2)]] =   new.mean.rel.abund[x.index[1:(rnt/2)]] *exp(0.5*beta1)
  # new.mean.rel.abund[x.index[(rnt/2+1):(rnt)]] =  new.mean.rel.abund[x.index[(rnt/2+1):(rnt)]] *exp(-0.5*beta1)
  # 
  # new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  
  
  new.mean.rel.abund = old.mean.rel.abund= count.dietswap.setup$mean.rel.abund
  
  ## Following Hawinkel et al. (2017)
  
  # fraction1 = 1/(beta1+1)
  # fraction1_selected = x.index[sample(1:length(x.index),round(fraction1*length(x.index)),replace=FALSE)]
  # fraction1_a = sum(old.mean.rel.abund[fraction1_selected])
  # fraction2_selected = x.index[!(x.index %in%fraction1_selected)]
  # fraction2_b = sum(old.mean.rel.abund[fraction2_selected])
  # 
  # new.mean.rel.abund[fraction1_selected] =   new.mean.rel.abund[fraction1_selected] * beta1
  # new.mean.rel.abund[fraction2_selected] =  new.mean.rel.abund[fraction2_selected] * ((fraction1_a/fraction2_b)*(1-beta1)+1)
  # 
  # new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  # table(new.mean.rel.abund/old.mean.rel.abund)
  
  ## Following Thas (personal communication 7/2/2025)
  fraction1 = 0.9
  fraction1_selected = x.index[sample(1:length(x.index),round(fraction1*length(x.index)),replace=FALSE)]
  fraction2_selected = x.index[!(x.index %in%fraction1_selected)]
  
  new.mean.rel.abund[fraction1_selected] =   old.mean.rel.abund[fraction1_selected] * 1/beta1
  new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] * 1.2
  
  sum0 = sum(old.mean.rel.abund[-x.index])
  sum1 = sum(new.mean.rel.abund[fraction1_selected])
  sum2 = sum(new.mean.rel.abund[fraction2_selected])
  table(new.mean.rel.abund/old.mean.rel.abund)
  gam_d = (1-sum0-sum1)/sum2
  new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] *  1.2 *gam_d
  
  
  
## High compositionality setting (no adjustment, so all LFC of all taxa change due to compositionality, setting 16 used)
  
  #new.mean.rel.abund[x.index] =   new.mean.rel.abund[x.index] * beta1
  #new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  
  new.lib.size = sample(min(count.dietswap.setup$lib.size):max(count.dietswap.setup$lib.size), size=n.data/2, replace=TRUE)
  
  count.dietswap.modified.1 <- try(MIDASim.modify(count.dietswap.setup,
                                             mean.rel.abund = old.mean.rel.abund,
                                             lib.size = new.lib.size))
  
  count.dietswap.modified.2 <- try(MIDASim.modify(count.dietswap.setup,
                                             mean.rel.abund = new.mean.rel.abund,
                                             lib.size = new.lib.size))
  
  
  if ("try-error" %in% class(count.dietswap.modified.2)){
    break()
  }else{
    for(sim in 1:100){
      temp1 <- MIDASim(count.dietswap.modified.1)
      temp2 <- MIDASim(count.dietswap.modified.2)
      
      data.all <- rbind(temp1$sim_count,temp2$sim_count)
      ra.all<- rbind(temp1$sim_rel,temp2$sim_rel)
      
      if (any(is.na(rowSums(data.all)))){
        next()
      }
      
      cens.prob <- colMeans(data.all == 0) #the proportion of samples where the count of that taxon is zero.
      
      trueDA <- rep(0, n.taxa)
      trueDA[x.index] <- 1
      colnames(data.all) <- paste0("otu", 1:ncol(data.all))
      midas.data <- list(trueDA, data.all, cens.prob, x.true)
      fn2 <- file.path(file.path("./Sim_datasets/Dietswap", fn),sprintf("sim%d.Rdata",sim))
      save(midas.data, file = fn2)
    }}}

#### Create phyloseq object ####
create_phyloseq <- function(midas_data) {  
  # otu table
  midas_otutable <- midas_data[[2]]  # Extract the OTU table (counts)
  rownames(midas_otutable) <- paste0("sample", 1:nrow(midas_otutable))  # Set row names for samples
  
  # taxa table
  n.taxa <- ncol(midas_otutable)  # Number of taxa is the number of columns in the OTU table
  midas_taxtable <- matrix(nrow = n.taxa, ncol = 2)  # Create a 2-column matrix for taxa info
  midas_taxtable[, 1:2] <- c(midas_data[[1]], midas_data[[3]])  # First column: DA_ind, Second column: cens.prob
  colnames(midas_taxtable) <- c("DA_ind", "cens.prob")  # Column names for the taxa table
  midas_taxtable <- as.data.frame(midas_taxtable)  # Convert to data frame
  midas_taxtable$isDA <- ifelse(midas_taxtable$DA_ind == 1, "TRUE", "FALSE")  # Add isDA column
  rownames(midas_taxtable) <- colnames(midas_otutable)  # Set taxa (OTU) names as rownames for taxa table
  midas_taxtable <- as.matrix(midas_taxtable)  # Convert to matrix (as required by phyloseq)
  
  # sample data
  midas_samdata <- as.data.frame(midas_data[[4]])  # Extract sample data (x.index values)
  rownames(midas_samdata) <- rownames(midas_otutable)  # Set row names for sample data (same as OTU table)
  colnames(midas_samdata) <- "x.index"  # Column name for the variable of interest
  midas_samdata$group <- ifelse(midas_samdata$x.index == 0.5, 1, 0)  # Create a 'group' column for the groupings
  
  # Create phyloseq object
  midas_phyloseq <- phyloseq(otu_table(midas_otutable, taxa_are_rows = FALSE),
                             sample_data(midas_samdata),
                             tax_table(midas_taxtable))
  
  return(midas_phyloseq)  # Return the phyloseq object
}
