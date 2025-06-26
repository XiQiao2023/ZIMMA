#' Summarize ZIMMA Posterior Output
#'
#' This function summarizes the posterior results from ZIMMA.
#'
#' @param ZIMMA.pos The output from ZIMMA.
#' @return The parameter's posterior inclusion probability, mean, median, highest posterior density.
#' @examples
#' ZIMMA.summary(posterior_list)
#' @import coda
#' @export
ZIMMA.summary = function(ZIMMA.pos = posterior_list){
  nP = dim(ZIMMA.pos$pos_beta)[3]
  nC.med = dim(ZIMMA.pos$pos_gamma)[2]
  nC.out = dim(ZIMMA.pos$pos_alpha)[2]-nP
  taxa.names = ZIMMA.pos$object@tax_table[,ZIMMA.pos$M.level]

  ###########################
  ########## betaT ##########
  ###########################
  betaT = cbind(
    apply(ZIMMA.pos$pos_kbeta, 2, mean),
    apply(ZIMMA.pos$pos_beta[,2,], 2, mean),
    apply(ZIMMA.pos$pos_beta[,2,], 2, median),
    HPDinterval(mcmc(ZIMMA.pos$pos_beta[,2,])))
  colnames(betaT) = c("PIP","mean","median","HPD.Lower","HPD.Upper")
  rownames(betaT) = taxa.names

  ###########################
  ######### gammaT ##########
  ###########################
  gammaT = cbind(
    apply(ZIMMA.pos$pos_kgamma, 2, mean),
    apply(ZIMMA.pos$pos_gamma[,2,], 2, mean),
    apply(ZIMMA.pos$pos_gamma[,2,], 2, median),
    HPDinterval(mcmc(ZIMMA.pos$pos_gamma[,2,])))
  colnames(gammaT) = c("PIP","mean","median","HPD.Lower","HPD.Upper")
  rownames(gammaT) = taxa.names

  ###########################
  ######### alphaM ##########
  ###########################
  alphaM = cbind(
    apply(ZIMMA.pos$pos_kalpha, 2, mean),
    apply(ZIMMA.pos$pos_alpha[,3:(2+nP)], 2, mean),
    apply(ZIMMA.pos$pos_alpha[,3:(2+nP)], 2, median),
    HPDinterval(mcmc(ZIMMA.pos$pos_alpha[,3:(2+nP)])))
  colnames(alphaM) = c("PIP","Mean","Median","HPD.Lower","HPD.Upper")
  rownames(alphaM) = taxa.names

  ###########################
  ########## alphaT #########
  ###########################
  alphaT = cbind(
    mean(ZIMMA.pos$pos_alpha[,2]),
    median(ZIMMA.pos$pos_alpha[,2]),
    HPDinterval(mcmc(ZIMMA.pos$pos_alpha[,2])))
  colnames(alphaT) = c("Mean","Median","HPD.Lower","HPD.Upper")

  ###########################
  ########## Tau ############
  ###########################
  tau = cbind(
    apply(ZIMMA.pos$pos_tau[,1,], 2, mean),
    apply(ZIMMA.pos$pos_tau[,1,], 2, median),
    HPDinterval(mcmc(ZIMMA.pos$pos_tau[,1,])))
  colnames(tau) = c("mean","median","HPD.Lower","HPD.Upper")
  rownames(tau) = taxa.names

  #######################################
  ########## Acceptance Rate ############
  #######################################
  AR = cbind(
    apply(ZIMMA.pos$pos_beta[,nC.med+1,], 2, mean),
    apply(ZIMMA.pos$pos_beta[,nC.med+1,], 2, sd),
    apply(ZIMMA.pos$pos_tau[,2,], 2, mean),
    apply(ZIMMA.pos$pos_tau[,2,], 2, sd)
  )
  colnames(AR) = c("Beta_AR_Mean","Beta_AR_SD","Tau_AR_Mean","Tau_AR_SD")
  rownames(AR) = taxa.names

  #######################################
  ########## Structure Zeros ############
  #######################################
  S_Zero <- apply(ZIMMA.pos$structure_Zero, c(2, 3), mean)
  colnames(S_Zero) <- taxa.names
  rownames(S_Zero) <- colnames(ZIMMA.pos$object@otu_table)

  ###############################
  ############ Output ###########
  ###############################
  out = list(betaT,gammaT,alphaM,alphaT,tau,AR,S_Zero)
  names(out) = c("betaT","gammaT","alphaM","alphaT","tau",
                 "Acceptance_Rate","Structural_Zero")
  out$M.level = ZIMMA.pos$M.level
  out$object = ZIMMA.pos$object

  return(out)
}
