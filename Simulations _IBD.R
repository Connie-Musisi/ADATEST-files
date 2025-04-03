## Data simulation using MIDAsim
# Load packages

library(MIDASim)
library(nleqslv)
library(survival)
library(microbiome)

# General setup using count.ibd data
data("count.ibd")
count.ibd.setup <- MIDASim.setup(count.ibd, mode = 'parametric')

# Which settings to explore

alln <- c(50, 200) # total sample size
allb <- c(1, 1.5, 3, 5)  # effect size
rnt <- c(20,50) # number of causal taxa

eg <- expand.grid(n = alln, beta = allb, rnt = rnt, 
                  stringsAsFactors = F)
iscen = 16
for(iscen in 1:nrow(eg)){

  condition <- eg[iscen,]
  n.data <- condition$n
  beta1 <- condition$beta # effect of the variable of interest
  rnt <- condition$rnt
  
  n.taxa = count.ibd.setup$n.taxa
  
  fn <- sprintf("Setting%d(highcomp)",iscen)
  dir.create(file.path("./Sim_datasets", fn), showWarnings = FALSE)
  
  
  # Top 60 of most abundant taxa are selected for possible modification
  m = 60
  top.index = order(count.ibd.setup$mean.rel.abund, decreasing =T)[1:m]
  x.index = sample(top.index, rnt, replace =F) # get causal taxa set: randomly select rnt=10/20 taxa from the most abundant 50 taxa
  
  # Group indicator
  x.true = c(rep(0, n.data/2), rep(1, n.data/2))
  x.true = x.true - mean(x.true)
  
  # Generate 100 datasets per scenario under consideration
  
  data.all = ra.all = matrix(NA, n.data, n.taxa)
  
  effect1 <- (x.true * beta1)
  new.mean.rel.abund = old.mean.rel.abund= count.ibd.setup$mean.rel.abund
  new.mean.rel.abund[x.index[1:(rnt/2)]] =   new.mean.rel.abund[x.index[1:(rnt/2)]] *exp(0.5*beta1)
  new.mean.rel.abund[x.index[(rnt/2+1):(rnt)]] =  new.mean.rel.abund[x.index[(rnt/2+1):(rnt)]] *exp(-0.5*beta1)
  
  new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  
  
  new.mean.rel.abund = old.mean.rel.abund= count.ibd.setup$mean.rel.abund
  
  ## Following Hawinkel et al. (2017)
  
   fraction1 = 1/(beta1+1)
   fraction1_selected = x.index[sample(1:length(x.index),round(fraction1*length(x.index)),replace=FALSE)]
   fraction1_a = sum(old.mean.rel.abund[fraction1_selected])
   fraction2_selected = x.index[!(x.index %in%fraction1_selected)]
   fraction2_b = sum(old.mean.rel.abund[fraction2_selected])
   
   new.mean.rel.abund[fraction1_selected] =   new.mean.rel.abund[fraction1_selected] * beta1
   new.mean.rel.abund[fraction2_selected] =  new.mean.rel.abund[fraction2_selected] * ((fraction1_a/fraction2_b)*(1-beta1)+1)
   
   new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
   table(new.mean.rel.abund/old.mean.rel.abund)

  # Following Thas (personal communication 7/2/2025) -> used for all settings to make sure only the selected taxa differ
   fraction1 = 0.9
   fraction1_selected = x.index[sample(1:length(x.index),round(fraction1*length(x.index)),replace=FALSE)]
   fraction2_selected = x.index[!(x.index %in%fraction1_selected)]
   
   new.mean.rel.abund[fraction1_selected] =   old.mean.rel.abund[fraction1_selected] * 1/beta1
   new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] * 1.2
   
   sum0 = sum(old.mean.rel.abund[-x.index])
   sum1 = sum(new.mean.rel.abund[fraction1_selected])
   sum2 = sum(new.mean.rel.abund[fraction2_selected])
   
   gam_d = (1-sum0-sum1)/sum2
   new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] *  1.2 *gam_d
   
  
  
## High compositionality setting (no adjustment, so all LFC of all taxa change due to compositionality, setting 16 used)
  
 # new.mean.rel.abund = old.mean.rel.abund= count.ibd.setup$mean.rel.abund
 # new.mean.rel.abund[x.index] =   new.mean.rel.abund[x.index] * beta1
 # new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  
  new.lib.size = sample(min(count.ibd.setup$lib.size):max(count.ibd.setup$lib.size), size=n.data/2, replace=TRUE)
  
  count.ibd.modified.1 <- try(MIDASim.modify(count.ibd.setup,
                                             mean.rel.abund = old.mean.rel.abund,
                                             lib.size = new.lib.size))
  
  count.ibd.modified.2 <- try(MIDASim.modify(count.ibd.setup,
                                             mean.rel.abund = new.mean.rel.abund,
                                             lib.size = new.lib.size))
  
  
  if ("try-error" %in% class(count.ibd.modified.2)){
    break()
  }else{
    for(sim in 1:100){
    temp1 <- MIDASim(count.ibd.modified.1)
    temp2 <- MIDASim(count.ibd.modified.2)
    
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
  fn2 <- file.path(file.path("./Sim_datasets", fn),sprintf("sim%d.Rdata",sim))
  save(midas.data, file = fn2)
  }}}


min.without.zero<-apply(ra.all,2,FUN=
                          function(x) {
                            if(length(x[x!=0])>0){min(x[x!=0])}
                            else{0.00001}
                          })

hist(log2((apply(ra.all[1:100,],2,mean)+min.without.zero)/(apply(ra.all[101:200,],2,mean)+min.without.zero))[which(trueDA==1)], nclass = 20,main=paste0("Before trimming, median LFC = ",median(log2((apply(ra.all[1:100,],2,mean)+min.without.zero)/(apply(ra.all[101:200,],2,mean)+min.without.zero))[which(trueDA==1)],na.rm=TRUE)))
hist(log2((apply(ra.all[1:100,],2,mean)+min.without.zero)/(apply(ra.all[101:200,],2,mean)+min.without.zero))[which(trueDA==0)], nclass = 20,main=paste0("Before trimming, median LFC = ",median(log2((apply(ra.all[1:100,],2,mean)+min.without.zero)/(apply(ra.all[101:200,],2,mean)+min.without.zero))[which(trueDA==0)],na.rm=TRUE)))




lfc2<-log2(colMeans(tmp[(n1+1):n,])/colMeans(tmp[1:n1,]))
lfc2_check = log2((apply(ra.all[101:200,],2,mean))/(apply(ra.all[1:100,],2,mean)))

tmp<-tax_table(simdata_filter)
hist(lfc2[as.data.frame(tmp)$isDA=="FALSE"], nclass = 20,main=paste0("Before trimming, median LFC = ",median(lfc2[as.data.frame(tmp)$isDA=="FALSE"],na.rm=TRUE)))
hist(lfc2_check[which(trueDA==0)], nclass = 20,main=paste0("Before trimming, median LFC = ",median(lfc2_check[which(trueDA==0)],na.rm=TRUE)))
hist(log2((apply(ra.all[101:200,],2,mean)+min.without.zero)/(apply(ra.all[1:100,],2,mean)+min.without.zero))[which(trueDA==0)], nclass = 20,main=paste0("Before trimming, median LFC = ",median(log2((apply(ra.all[101:200,],2,mean)+min.without.zero)/(apply(ra.all[1:100,],2,mean)+min.without.zero))[which(trueDA==0)],na.rm=TRUE)))
hist(log2((apply(ra.all[101:200,],2,mean))/(apply(ra.all[1:100,],2,mean)))[which(trueDA==0)], nclass = 20,main=paste0("Before trimming, median LFC = ",median(log2((apply(ra.all[101:200,],2,mean))/(apply(ra.all[1:100,],2,mean)))[which(trueDA==0)],na.rm=TRUE)))

tmp

max(min.without.zero)

