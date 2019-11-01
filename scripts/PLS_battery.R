#### this needs some cleaning & testing (on a baby example)
  ### specifically check the BSRs are correct and also 

rm(list=ls())
library(oro.nifti)
library(TInPosition)
library(neuroim)
library(ours)
library(abind)
library(pheatmap)
library(corrplot)
library(gplots)

load('../../rdata/aggregate.design_2mm.rda')
load('../../rdata/spm.betas_2mm.rda')

aggregate.design.df <- as.data.frame(aggregate.design)
rownames(aggregate.design.df) <- NULL


# on.global.mask <- readNIfTI('./masks/HarvardOxford-cort_subcort_cereb_binary_mask-thr0_3mm.nii')
global.mask <- loadVolume('./masks/HarvardOxford-cort_subcort_cereb_binary_mask-thr0_2mm.nii')
global.mask.ind <- which(global.mask@.Data==1)



empty.brain <- array(0,dim=global.mask@space@Dim)

acc.mask <- loadVolume('./masks/HO_accumbens_2mm.nii')
acc.mask.ind <- which(acc.mask@.Data==1 & global.mask@.Data==1)
acc.columns <- which(global.mask.ind %in% acc.mask.ind)

amyg.mask <- loadVolume('./masks/HO_amygdala_2mm.nii')
amyg.mask.ind <- which(amyg.mask@.Data==1 & global.mask@.Data==1)
amyg.columns <- which(global.mask.ind %in% amyg.mask.ind)

mf.mask <- loadVolume('./masks/HO_medFront_2mm.nii')
mf.mask.ind <- which(mf.mask@.Data==1 & global.mask@.Data==1)
mf.columns <- which(global.mask.ind %in% mf.mask.ind)

interaction.only <- model.matrix(~ GROUP:CONDITION, data=aggregate.design.df, contrasts.arg = lapply(aggregate.design.df, contrasts, contrasts=FALSE))[,-1]
# me.inter <- model.matrix(~ GROUP*CONDITION, data=aggregate.design.df, contrasts.arg = lapply(aggregate.design.df, contrasts, contrasts=FALSE))[,-1]

iterations <- 250
pls1.res <- tepPLS(interaction.only, spm.betas, center1 = T, center2 = T, scale1 = F, scale2 = F, graphs=F, DESIGN = interaction.only, make_design_nominal = F)

acc.mean.fjs <- colMeans(pls1.res$TExPosition.Data$fj[acc.columns,])
amyg.mean.fjs <- colMeans(pls1.res$TExPosition.Data$fj[amyg.columns,])
mf.mean.fjs <- colMeans(pls1.res$TExPosition.Data$fj[mf.columns,])

mean.fjs <- rbind(acc.mean.fjs, amyg.mean.fjs, mf.mean.fjs)
rownames(mean.fjs) <- c("Acc","Amyg","MF")

acc.cors <- cor(rowMeans(spm.betas[,acc.columns]), pls1.res$TExPosition.Data$lx, method = "spearman")
amyg.cors <- cor(rowMeans(spm.betas[,amyg.columns]), pls1.res$TExPosition.Data$lx, method = "spearman")
mf.cors <- cor(rowMeans(spm.betas[,mf.columns]), pls1.res$TExPosition.Data$lx, method = "spearman")


  ## grp1 is increase + equal range
grp1.to.sample.from <- which(aggregate.design[,"CONDITION"]=="gain_increase" & aggregate.design[,"GROUP"]=="equalRange")
  ## grp2 is increase + equal indifference
grp2.to.sample.from <- which(aggregate.design[,"CONDITION"]=="gain_increase" & aggregate.design[,"GROUP"]=="equalIndifference")
  ## to get loss group, just +1 to each of these.
  
### these things all need some double checks...

fi.bsrs.m2 <- fi.bsrs.mean <- matrix(0, nrow(pls1.res$TExPosition.Data$fi), ncol(pls1.res$TExPosition.Data$fi)) ->
    fi.bsrs.m2_gain.range -> fi.bsrs.mean_gain.range -> fi.bsrs.m2_gain.indiff -> fi.bsrs.mean_gain.indiff ->
    fi.bsrs.m2_loss.range -> fi.bsrs.mean_loss.range -> fi.bsrs.m2_loss.indiff -> fi.bsrs.mean_loss.indiff

fj.bsrs.m2 <- fj.bsrs.mean <- matrix(0, nrow(pls1.res$TExPosition.Data$fj), ncol(pls1.res$TExPosition.Data$fj)) ->
    fj.bsrs.m2_gain.range -> fj.bsrs.mean_gain.range -> fj.bsrs.m2_gain.indiff -> fj.bsrs.mean_gain.indiff ->
    fj.bsrs.m2_loss.range -> fj.bsrs.mean_loss.range -> fj.bsrs.m2_loss.indiff -> fj.bsrs.mean_loss.indiff ->
    fj.bsrs.m2_loss -> fj.bsrs.mean_loss ->
    fj.bsrs.m2_loss_v2 -> fj.bsrs.mean_loss_v2


fi.top.5 <- array(-Inf, dim=c(nrow(pls1.res$TExPosition.Data$fi),ncol(pls1.res$TExPosition.Data$fi), floor(round(iterations * .025))+3 )) # ->
  # fi.top.5_gain.range -> fi.top.5_gain.indiff ->
  # fi.top.5_loss.range -> fi.top.5_loss.indiff

fi.bottom.5 <- array(Inf, dim=c(nrow(pls1.res$TExPosition.Data$fi),ncol(pls1.res$TExPosition.Data$fi), floor(round(iterations * .025))+3 )) # ->
  # fi.bottom.5_gain.range -> fi.bottom.5_gain.indiff ->
  # fi.bottom.5_loss.range -> fi.bottom.5_loss.indiff

# fj.top.5 <- array(-Inf, dim=c(nrow(pls1.res$TExPosition.Data$fj),ncol(pls1.res$TExPosition.Data$fj), floor(round(iterations * .025))+3 )) #matrix(-Inf, nrow(pls1.res$TExPosition.Data$fi), floor(round(iterations * .025)))
# fj.bottom.5 <- array(Inf, dim=c(nrow(pls1.res$TExPosition.Data$fj),ncol(pls1.res$TExPosition.Data$fj), floor(round(iterations * .025))+3 ))

# fi.arr <- array(NA,dim=c(nrow(pls1.res$TExPosition.Data$fi),ncol(pls1.res$TExPosition.Data$fi),iterations)) # -> fi.split.arr

# lx.arr <- ly.arr <- array(NA,dim=c(nrow(pls1.res$TExPosition.Data$lx),ncol(pls1.res$TExPosition.Data$lx),iterations))
lx.cors <- ly.cors <- array(NA,dim=c(ncol(pls1.res$TExPosition.Data$lx),ncol(pls1.res$TExPosition.Data$lx),iterations))
# fj.split.cors <- matrix(NA,iterations,ncol(pls1.res$TExPosition.Data$lx))

sh1.indices <- matrix(NA, iterations, nrow(interaction.only)/2)

for(i in 1:iterations){
  
  ## resampling: break the design into gain/loss (which automatically links to participants)
      ## then resample within for the group?
      ## or just resample from within gain constrained by group; then duplicate that sample for loss
  
  ## conditional bootstrap
  grp1.boot.indices_init <- sample(grp1.to.sample.from,replace = T)
    grp1.boot.indices <- c(grp1.boot.indices_init, grp1.boot.indices_init+1)
    
  grp2.boot.indices_init <- sample(grp2.to.sample.from,replace = T)
    grp2.boot.indices <- c(grp2.boot.indices_init, grp2.boot.indices_init+1)
    
  boot.indices <- sort(c(grp1.boot.indices, grp2.boot.indices))

  boot.interaction.only <- interaction.only[boot.indices, ]
  
  boot.R <- t(expo.scale(boot.interaction.only, center = pls1.res$TExPosition.Data$data1.norm$center, scale = pls1.res$TExPosition.Data$data1.norm$scale) ) %*% expo.scale( spm.betas[boot.indices,], center = pls1.res$TExPosition.Data$data2.norm$center, scale = pls1.res$TExPosition.Data$data2.norm$scale)
  boot.fi <- boot.R %*% pls1.res$TExPosition.Data$pdq$q
    boot.fi[is.nan(boot.fi)] <- 0
  boot.fj <- t(boot.R) %*% pls1.res$TExPosition.Data$pdq$p
    boot.fj[is.nan(boot.fj)] <- 0
    
  ## ok so here I need to compute these boot.fjs; then probably trim down the subsequent code.
  boot.R_gain.range <- t(expo.scale(interaction.only[grp1.boot.indices_init,],center = pls1.res$TExPosition.Data$data1.norm$center, scale = pls1.res$TExPosition.Data$data1.norm$scale)) %*% expo.scale( spm.betas[grp1.boot.indices_init,], center = pls1.res$TExPosition.Data$data2.norm$center, scale = pls1.res$TExPosition.Data$data2.norm$scale)
  boot.R_loss.range <- t(expo.scale(interaction.only[grp1.boot.indices_init+1,],center = pls1.res$TExPosition.Data$data1.norm$center, scale = pls1.res$TExPosition.Data$data1.norm$scale)) %*% expo.scale( spm.betas[grp1.boot.indices_init+1,], center = pls1.res$TExPosition.Data$data2.norm$center, scale = pls1.res$TExPosition.Data$data2.norm$scale)
  boot.R_gain.indiff <- t(expo.scale(interaction.only[grp2.boot.indices_init,],center = pls1.res$TExPosition.Data$data1.norm$center, scale = pls1.res$TExPosition.Data$data1.norm$scale)) %*% expo.scale( spm.betas[grp2.boot.indices_init,], center = pls1.res$TExPosition.Data$data2.norm$center, scale = pls1.res$TExPosition.Data$data2.norm$scale)
  boot.R_loss.indiff <- t(expo.scale(interaction.only[grp2.boot.indices_init+1,],center = pls1.res$TExPosition.Data$data1.norm$center, scale = pls1.res$TExPosition.Data$data1.norm$scale)) %*% expo.scale( spm.betas[grp2.boot.indices_init+1,], center = pls1.res$TExPosition.Data$data2.norm$center, scale = pls1.res$TExPosition.Data$data2.norm$scale)

    ## revisit this 
      ## also model it better...
  boot.R_loss <- t((matrix(c(0,0,-sqrt(4),sqrt(4)),nrow(interaction.only),ncol(interaction.only),byrow=T) * expo.scale(interaction.only,scale=F))[boot.indices,]) %*% expo.scale( spm.betas[boot.indices,], center = pls1.res$TExPosition.Data$data2.norm$center, scale = pls1.res$TExPosition.Data$data2.norm$scale)
  # boot.R_loss <- t( (matrix(c(0,0,-1,1),nrow(interaction.only),ncol(interaction.only),byrow=T) * interaction.only[boot.indices,]) ) %*% expo.scale( spm.betas[boot.indices,], center = pls1.res$TExPosition.Data$data2.norm$center, scale = pls1.res$TExPosition.Data$data2.norm$scale)
  # 
    ## specifically to get brains per group w.r.t the omnibus discriminant/interaction model.
    boot.fj_gain.range <- t(boot.R_gain.range) %*% pls1.res$TExPosition.Data$pdq$p
      boot.fj_gain.range[is.nan(boot.fj_gain.range)] <- 0
    boot.fj_loss.range <- t(boot.R_loss.range) %*% pls1.res$TExPosition.Data$pdq$p
      boot.fj_loss.range[is.nan(boot.fj_loss.range)] <- 0
    boot.fj_gain.indiff <- t(boot.R_gain.indiff) %*% pls1.res$TExPosition.Data$pdq$p
      boot.fj_gain.indiff[is.nan(boot.fj_gain.indiff)] <- 0
    boot.fj_loss.indiff <- t(boot.R_loss.indiff) %*% pls1.res$TExPosition.Data$pdq$p
      boot.fj_loss.indiff[is.nan(boot.fj_loss.indiff)] <- 0

      ## revisit this      
    boot.fj_loss <- t(boot.R_loss) %*% pls1.res$TExPosition.Data$pdq$p
      boot.fj_loss[is.nan(boot.fj_loss)] <- 0
  
      
      
  fi.delta <- boot.fi - fi.bsrs.mean
  fi.bsrs.mean <- fi.bsrs.mean + (fi.delta / i)
  fi.bsrs.m2 <- fi.bsrs.m2 + fi.delta * (boot.fi - fi.bsrs.mean)

  fj.delta <- boot.fj - fj.bsrs.mean
  fj.bsrs.mean <- fj.bsrs.mean + (fj.delta / i)
  fj.bsrs.m2 <- fj.bsrs.m2 + fj.delta * (boot.fj - fj.bsrs.mean)
  
    fj.delta_gain.range <- boot.fj_gain.range - fj.bsrs.mean_gain.range
    fj.bsrs.mean_gain.range <- fj.bsrs.mean_gain.range + (fj.delta_gain.range / i)
    fj.bsrs.m2_gain.range <- fj.bsrs.m2_gain.range + fj.delta_gain.range * (boot.fj_gain.range - fj.bsrs.mean_gain.range)
    
    fj.delta_loss.range <- boot.fj_loss.range - fj.bsrs.mean_loss.range
    fj.bsrs.mean_loss.range <- fj.bsrs.mean_loss.range + (fj.delta_loss.range / i)
    fj.bsrs.m2_loss.range <- fj.bsrs.m2_loss.range + fj.delta_loss.range * (boot.fj_loss.range - fj.bsrs.mean_loss.range)
    
    fj.delta_gain.indiff <- boot.fj_gain.indiff - fj.bsrs.mean_gain.indiff
    fj.bsrs.mean_gain.indiff <- fj.bsrs.mean_gain.indiff + (fj.delta_gain.indiff / i)
    fj.bsrs.m2_gain.indiff <- fj.bsrs.m2_gain.indiff + fj.delta_gain.indiff * (boot.fj_gain.indiff - fj.bsrs.mean_gain.indiff)
    
    fj.delta_loss.indiff <- boot.fj_loss.indiff - fj.bsrs.mean_loss.indiff
    fj.bsrs.mean_loss.indiff <- fj.bsrs.mean_loss.indiff + (fj.delta_loss.indiff / i)
    fj.bsrs.m2_loss.indiff <- fj.bsrs.m2_loss.indiff + fj.delta_loss.indiff * (boot.fj_loss.indiff - fj.bsrs.mean_loss.indiff)
    
    
    fj.delta_loss <- boot.fj_loss - fj.bsrs.mean_loss
    fj.bsrs.mean_loss <- fj.bsrs.mean_loss + (fj.delta_loss / i)
    fj.bsrs.m2_loss <- fj.bsrs.m2_loss + fj.delta_loss * (boot.fj_loss - fj.bsrs.mean_loss)
    
  
  ## will do the binding and sorting later...
  fi.bottom.5[,,dim(fi.bottom.5)[3]] <- fi.top.5[,,dim(fi.top.5)[3]] <- boot.fi
  fi.bottom.5 <- aperm(apply (fi.bottom.5, 1:2, sort), c(2,3,1))
  fi.top.5 <- aperm(apply (fi.top.5, 1:2, sort, decreasing=T), c(2,3,1))
  
  
  ## conditional split-half
  grp1.split.1.indices <- sample(grp1.to.sample.from, replace = F, floor(length(grp1.to.sample.from)/2))
    grp1.split.1.indices <- c(grp1.split.1.indices, grp1.split.1.indices+1)
    
  grp2.split.1.indices <- sample(grp2.to.sample.from, replace = F, floor(length(grp2.to.sample.from)/2))
    grp2.split.1.indices <- c(grp2.split.1.indices, grp2.split.1.indices+1)
    
  sh1.indices[i,] <- split.1.indices <- sort(c(grp1.split.1.indices, grp2.split.1.indices))
  split.2.indices <- sort(setdiff(1:nrow(aggregate.design), split.1.indices))

  split.1.interaction.only <- interaction.only[split.1.indices,]
  split.2.interaction.only <- interaction.only[split.2.indices,]
  
  split.1.pls <- tepPLS(split.1.interaction.only, spm.betas[split.1.indices,], center1 = T, center2 = T, scale1 = F, scale2 = F, graphs=F)
  split.2.pls <- tepPLS(split.2.interaction.only, spm.betas[split.2.indices,], center1 = T, center2 = T, scale1 = F, scale2 = F, graphs=F)

  ## hang on to all of these...  
  pred.lx1.from.2 <- expo.scale(split.1.interaction.only, center=split.2.pls$TExPosition.Data$data1.norm$center, scale=split.2.pls$TExPosition.Data$data1.norm$scale) %*% split.2.pls$TExPosition.Data$pdq$p
  pred.ly1.from.2 <- expo.scale(spm.betas[split.1.indices,], center=split.2.pls$TExPosition.Data$data2.norm$center, scale=split.2.pls$TExPosition.Data$data2.norm$scale) %*% split.2.pls$TExPosition.Data$pdq$q
  
  pred.lx2.from.1 <- expo.scale(split.2.interaction.only, center=split.1.pls$TExPosition.Data$data1.norm$center, scale=split.1.pls$TExPosition.Data$data1.norm$scale) %*% split.1.pls$TExPosition.Data$pdq$p
  pred.ly2.from.1 <- expo.scale(spm.betas[split.2.indices,], center=split.1.pls$TExPosition.Data$data2.norm$center, scale=split.1.pls$TExPosition.Data$data2.norm$scale) %*% split.1.pls$TExPosition.Data$pdq$q
  
  lx.cors[,,i] <- ((cor(pred.lx1.from.2, split.1.pls$TExPosition.Data$lx)^2) + (cor(pred.lx2.from.1, split.2.pls$TExPosition.Data$lx)^2)) / 2
  ly.cors[,,i] <- ((cor(pred.ly1.from.2, split.1.pls$TExPosition.Data$ly)^2) + (cor(pred.ly2.from.1, split.2.pls$TExPosition.Data$ly)^2)) / 2
  
  print(i)
  
}

fi.variance = fi.bsrs.m2/iterations
fi.bsrs = fi.bsrs.mean/sqrt(fi.variance)

fj.variance = fj.bsrs.m2/iterations
fj.bsrs = fj.bsrs.mean/sqrt(fj.variance)

  fj.variance_gain.range = fj.bsrs.m2_gain.range/iterations
  fj.bsrs_gain.range = fj.bsrs.mean_gain.range/sqrt(fj.variance_gain.range)
  
  fj.variance_loss.range = fj.bsrs.m2_loss.range/iterations
  fj.bsrs_loss.range = fj.bsrs.mean_loss.range/sqrt(fj.variance_loss.range)
  
  fj.variance_gain.indiff = fj.bsrs.m2_gain.indiff/iterations
  fj.bsrs_gain.indiff = fj.bsrs.mean_gain.indiff/sqrt(fj.variance_gain.indiff)
  
  fj.variance_loss.indiff = fj.bsrs.m2_loss.indiff/iterations
  fj.bsrs_loss.indiff = fj.bsrs.mean_loss.indiff/sqrt(fj.variance_loss.indiff)
  
  fj.variance_loss = fj.bsrs.m2_loss/iterations
  fj.bsrs_loss = fj.bsrs.mean_loss/sqrt(fj.variance_loss)

  
  ## split-half tells us about stability of components.
    ## bootstrap tells us about stability of variables for entire model
    ## bootstrap also tells us about stability of voxels for each subgroup w.r.t. the moel.

  colnames(interaction.only) <- rownames(pls1.res$TExPosition.Data$fi) <- c("GAIN.INDIFF","GAIN.RANGE","LOSS.INDIFF","LOSS.RANGE")
  corrplot(apply(ly.cors, c(1,2), mean), method = "number")
  hist(ly.cors[1,1,],breaks=50)
  barplot2(pls1.res$TExPosition.Data$fi[,1], plot.ci = T, ci.u = fi.bottom.5[,1,(dim(fi.bottom.5)[3]-1)], ci.l = fi.top.5[,1,(dim(fi.top.5)[3]-1)])
  
  interact.cols <- createColorVectorsByDesign(interaction.only)
  prettyPlot(cbind(pls1.res$TExPosition.Data$ly[,1],pls1.res$TExPosition.Data$lx[,1]), col=interact.cols$oc, asp=NA, ylab="Design scores",xlab="Brain scores", dev.new = F)
  legend("topleft", legend=c("GAIN.INDIFF","GAIN.RANGE","LOSS.INDIFF","LOSS.RANGE"), col=interact.cols$gc, pch=20, bty = 0)
  
  ## here get the intersection of global BSRS and groups (for hypotheses)
    ## then check against the columns (e.g,. acc.columns) and then place those values back into the indices. 
   
  ## omnibus
  omnibus.thresholded <- omnibus.unthresholded <- array(NA, dim=c(global.mask@space@Dim))
  omnibus.unthresholded[global.mask.ind] <- fj.bsrs[,1]
  omnibus.thresholded[global.mask.ind[which(abs(fj.bsrs[,1])>3)] ] <- fj.bsrs[which(abs(fj.bsrs[,1])>3),1]

  writeVolume(
    BrainVolume(omnibus.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
  , fileName = "./Images/brains/omnibus_unthresholded_2mm.nii.gz")

  writeVolume(
    BrainVolume(omnibus.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/omnibus_thresholded_2mm.nii.gz")

  
  writeVolume(
    BrainVolume(omnibus.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/omnibus_unthresh.nii.gz")
  
  writeVolume(
    BrainVolume(omnibus.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/omnibus_thresh.nii.gz")

  ## specific groups (for hypotheses)

  ## post-hoc gain indiff range
  fj.bsrs_gain.indiff.intersect.thresholded <- fj.bsrs_gain.indiff.thresholded <- fj.bsrs_gain.indiff.unthresholded <- array(NA, dim=c(global.mask@space@Dim))
  fj.bsrs_gain.indiff.unthresholded[global.mask.ind] <- fj.bsrs_gain.indiff[,1]
  fj.bsrs_gain.indiff.thresholded[global.mask.ind[which(abs(fj.bsrs_gain.indiff[,1])>3)] ] <- fj.bsrs_gain.indiff[which(abs(fj.bsrs_gain.indiff[,1])>3),1]
  fj.bsrs_gain.indiff.intersect.thresholded[global.mask.ind[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_gain.indiff[,1]) > 3 )] ] <- fj.bsrs_gain.indiff[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_gain.indiff[,1]) > 3 ),1]

  writeVolume(
    BrainVolume(fj.bsrs_gain.indiff.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_gain.indiff.unthresholded_2mm.nii.gz")

  writeVolume(
    BrainVolume(fj.bsrs_gain.indiff.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_gain.indiff.thresholded_2mm.nii.gz")


  writeVolume(
    BrainVolume(fj.bsrs_gain.indiff.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_gain.indiff.intersect.thresholded_2mm.nii.gz")
  
  
  
  writeVolume(
    BrainVolume(fj.bsrs_gain.indiff.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo1_unthresh.nii.gz")
  writeVolume(
    BrainVolume(fj.bsrs_gain.indiff.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo3_unthresh.nii.gz")
  
  writeVolume(
    BrainVolume(fj.bsrs_gain.indiff.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo1_thresh.nii.gz")
  writeVolume(
    BrainVolume(fj.bsrs_gain.indiff.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo3_thresh.nii.gz")


  ## post-hoc gain equal range
  fj.bsrs_gain.range.intersect.thresholded <- fj.bsrs_gain.range.thresholded <- fj.bsrs_gain.range.unthresholded <- array(NA, dim=c(global.mask@space@Dim))
  fj.bsrs_gain.range.unthresholded[global.mask.ind] <- fj.bsrs_gain.range[,1]
  fj.bsrs_gain.range.thresholded[global.mask.ind[which(abs(fj.bsrs_gain.range[,1])>3)] ] <- fj.bsrs_gain.range[which(abs(fj.bsrs_gain.range[,1])>3),1]
  fj.bsrs_gain.range.intersect.thresholded[global.mask.ind[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_gain.range[,1]) > 3 )] ] <- fj.bsrs_gain.range[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_gain.range[,1]) > 3 ),1]

  writeVolume(
    BrainVolume(fj.bsrs_gain.range.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_gain.range.unthresholded_2mm.nii.gz")

  writeVolume(
    BrainVolume(fj.bsrs_gain.range.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_gain.range.thresholded_2mm.nii.gz")


  writeVolume(
    BrainVolume(fj.bsrs_gain.range.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_gain.range.intersect.thresholded_2mm.nii.gz")
  
  
  
  writeVolume(
    BrainVolume(fj.bsrs_gain.range.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo2_unthresh.nii.gz")
  writeVolume(
    BrainVolume(fj.bsrs_gain.range.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo4_unthresh.nii.gz")
  
  writeVolume(
    BrainVolume(fj.bsrs_gain.range.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo2_thresh.nii.gz")
  writeVolume(
    BrainVolume(fj.bsrs_gain.range.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo4_thresh.nii.gz")


  ## post-hoc loss indiff range
  fj.bsrs_loss.indiff.intersect.thresholded <- fj.bsrs_loss.indiff.thresholded <- fj.bsrs_loss.indiff.unthresholded <- array(NA, dim=c(global.mask@space@Dim))
  fj.bsrs_loss.indiff.unthresholded[global.mask.ind] <- fj.bsrs_loss.indiff[,1]
  fj.bsrs_loss.indiff.thresholded[global.mask.ind[which(abs(fj.bsrs_loss.indiff[,1])>3)] ] <- fj.bsrs_loss.indiff[which(abs(fj.bsrs_loss.indiff[,1])>3),1]
  fj.bsrs_loss.indiff.intersect.thresholded[global.mask.ind[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss.indiff[,1]) > 3 )] ] <- fj.bsrs_loss.indiff[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss.indiff[,1]) > 3 ),1]

  writeVolume(
    BrainVolume(fj.bsrs_loss.indiff.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_loss.indiff.unthresholded_2mm.nii.gz")

  writeVolume(
    BrainVolume(fj.bsrs_loss.indiff.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_loss.indiff.thresholded_2mm.nii.gz")


  writeVolume(
    BrainVolume(fj.bsrs_loss.indiff.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_loss.indiff.intersect.thresholded_2mm.nii.gz")
  
  
  
  writeVolume(
    BrainVolume(fj.bsrs_loss.indiff.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo5_unthresh.nii.gz")
  writeVolume(
    BrainVolume(fj.bsrs_loss.indiff.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo7_unthresh.nii.gz")
  
  writeVolume(
    BrainVolume(fj.bsrs_loss.indiff.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo5_thresh.nii.gz")
  writeVolume(
    BrainVolume(fj.bsrs_loss.indiff.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo7_thresh.nii.gz")


  ## post-hoc loss equal range
  fj.bsrs_loss.range.intersect.thresholded <- fj.bsrs_loss.range.thresholded <- fj.bsrs_loss.range.unthresholded <- array(NA, dim=c(global.mask@space@Dim))
  fj.bsrs_loss.range.unthresholded[global.mask.ind] <- fj.bsrs_loss.range[,1]
  fj.bsrs_loss.range.thresholded[global.mask.ind[which(abs(fj.bsrs_loss.range[,1])>3)] ] <- fj.bsrs_loss.range[which(abs(fj.bsrs_loss.range[,1])>3),1]
  fj.bsrs_loss.range.intersect.thresholded[global.mask.ind[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss.range[,1]) > 3 )] ] <- fj.bsrs_loss.range[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss.range[,1]) > 3 ),1]

  writeVolume(
    BrainVolume(fj.bsrs_loss.range.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_loss.range.unthresholded_2mm.nii.gz")

  writeVolume(
    BrainVolume(fj.bsrs_loss.range.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_loss.range.thresholded_2mm.nii.gz")


  writeVolume(
    BrainVolume(fj.bsrs_loss.range.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_loss.range.intersect.thresholded_2mm.nii.gz")

  
  
  writeVolume(
    BrainVolume(fj.bsrs_loss.range.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo6_unthresh.nii.gz")
  writeVolume(
    BrainVolume(fj.bsrs_loss.range.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo8_unthresh.nii.gz")
  
  
  writeVolume(
    BrainVolume(fj.bsrs_loss.range.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo6_thresh.nii.gz")
  
  writeVolume(
    BrainVolume(fj.bsrs_loss.range.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo8_thresh.nii.gz")

  ## post-hoc loss equal range
  fj.bsrs_loss.intersect.thresholded <- fj.bsrs_loss.thresholded <- fj.bsrs_loss.unthresholded <- array(NA, dim=c(global.mask@space@Dim))
  fj.bsrs_loss.unthresholded[global.mask.ind] <- fj.bsrs_loss[,1]
  fj.bsrs_loss.thresholded[global.mask.ind[which(abs(fj.bsrs_loss[,1])>3)] ] <- fj.bsrs_loss[which(abs(fj.bsrs_loss[,1])>3),1]
  fj.bsrs_loss.intersect.thresholded[global.mask.ind[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss[,1]) > 3 )] ] <- fj.bsrs_loss[which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss[,1]) > 3 ),1]

  writeVolume(
    BrainVolume(fj.bsrs_loss.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_loss.unthresholded.nii_2mm.gz")

  writeVolume(
    BrainVolume(fj.bsrs_loss.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_loss.thresholded.nii_2mm.gz")


  writeVolume(
    BrainVolume(fj.bsrs_loss.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/fj.bsrs_loss.intersect.thresholded.nii_2mm.gz")
  
  
  
  writeVolume(
    BrainVolume(fj.bsrs_loss.unthresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo9_unthresh.nii.gz")
  
  writeVolume(
    BrainVolume(fj.bsrs_loss.intersect.thresholded,
                BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
    , fileName = "./Images/brains/forNeuroVault/hypo9_thresh.nii.gz")
  
  
  # ## (1); yes?
  # hyp.1.voxels <- intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_gain.indiff[,1]) > 3 ), mf.columns)
  # # intersect(which( abs(fj.bsrs_gain.indiff[,1]) > 3 ), mf.columns)
  # pheatmap(t(interaction.only) %*% spm.betas[,hyp.1.voxels], cluster_rows = F, cluster_cols = F)
    
  # hyp.1.unthresolded <- hyp.1.thresholded <- array(NA, dim=c(global.mask@space@Dim))
  #   ### this makes no sense...; I have to just export the bsrs but we have to look at the betas to answer the question.
  # # quick diagnostic.  
  # plot(cbind(colMeans((spm.betas[which(interaction.only[,"GAIN.INDIFF"] == 1),hyp.1.voxels])), fj.bsrs_gain.indiff[hyp.1.voxels,1]))
  # # abs(fj.bsrs_gain.indiff[hyp.1.voxels,1]) * sign(colMeans(spm.betas[which(interaction.only[,"GAIN.INDIFF"] == 1),hyp.1.voxels]))
  # # fj.bsrs_gain.indiff[hyp.1.voxels,1]
  # hyp.1.unthresolded[global.mask.ind] <- fj.bsrs_gain.indiff[,1]
  # 
  # writeVolume(
  #   BrainVolume(hyp.1.unthresolded, 
  #               BrainSpace(Dim=c(global.mask@space@Dim), spacing = global.mask@space@spacing, origin = global.mask@space@origin, axes = global.mask@space@axes, trans = global.mask@space@trans))
  #   , fileName = "./Images/brains/hyp_1_unthresholded.nii.gz")
  
  
  ## (2); no
  intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_gain.range[,1]) > 3 ), mf.columns)
  intersect(which( abs(fj.bsrs_gain.range[,1]) > 3 ), mf.columns)
  
  ## (3); no
  intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_gain.indiff[,1]) > 3 ), acc.columns)
  intersect(which( abs(fj.bsrs_gain.indiff[,1]) > 3 ), acc.columns)
  
  ## (4); no
  intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_gain.range[,1]) > 3 ), acc.columns)
  intersect(which( abs(fj.bsrs_gain.range[,1]) > 3 ), acc.columns)
  
  
  ## (5); yes
  intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss.indiff[,1]) > 3 ), mf.columns)
  intersect(which( abs(fj.bsrs_loss.indiff[,1]) > 3 ), mf.columns)
  pheatmap(t(interaction.only) %*% spm.betas[,intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss.indiff[,1]) > 3 ), mf.columns)], cluster_rows = F, cluster_cols = F)
  
  ## (6); no
  intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss.range[,1]) > 3 ), mf.columns)
  intersect(which( abs(fj.bsrs_loss.range[,1]) > 3 ), mf.columns)
  
  ## (7); no
  intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss.indiff[,1]) > 3 ), amyg.columns)
  intersect(which( abs(fj.bsrs_loss.indiff[,1]) > 3 ), amyg.columns)
  
  ## (8); no
  intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss.range[,1]) > 3 ), amyg.columns)
  intersect(which( abs(fj.bsrs_loss.range[,1]) > 3 ), amyg.columns)
  
  ## (9); no
    ## this one is effectively a comparison based on # 7 & 8.
  intersect(which( abs(fj.bsrs[,1]) > 3 & abs(fj.bsrs_loss[,1]) > 3 ), amyg.columns)
  intersect(which( abs(fj.bsrs_loss[,1]) > 3 ), amyg.columns)
  
  
  
