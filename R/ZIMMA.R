#' ZIMMA: Zero-Inflated Microbiome Mediation Analysis
#'
#' This function performs Bayesian dual mediation analysis for microbiome data, accounting for zero inflation in taxa counts.
#' It models both prevalence and abundance components of each microbial feature to identify mediators of treatment effects on outcomes.
#'
#' @param object A phyloseq object containing microbiome count data and associated metadata.
#' @param Outcome Character. The name of the outcome variable in the sample metadata.
#' @param Treat Character. The name of the treatment/exposure variable in the sample metadata.
#' @param C.med Optional. A vector of covariate names used as confounders in the mediator model.
#' @param C.out Optional. A vector of covariate names used as confounders in the outcome model.
#' @param M.level Character. Taxonomic level to use (e.g., "Genus", "Family").
#' @param n_iter Integer. Total number of MCMC iterations.
#' @param burn_in Integer. Number of burn-in iterations to discard.
#' @param Size Character. Library size estimation method (e.g., "GMPR", "TSS", "RLEpseudo").
#' @param pseudo Numeric. Pseudo-count added during library size estimation to avoid division by zero.
#'
#' @return A list containing parameter posterior samples, and the phyloseq object.
#' @examples
#' ZIMMA(object = physeq, Outcome = Outcome, Treat = Treat, C.med = NULL, C.out = NULL, M.level = M.level,n_iter = 100, burn_in = 20,Size = "RLEpseudo", pseudo = 0.)
#' @export
ZIMMA = function(object = physeq, Outcome = Outcome, Treat = Treat,
                 C.med = NULL, C.out = NULL, M.level = M.level,
                 n_iter = 100, burn_in = 20,Size = "RLEpseudo", pseudo = 0.5){

  physeq_level = tax_glom(object,taxrank = M.level)
  metadata = sample_data(physeq_level)

  Treat <- metadata[[Treat]]
  Y <- metadata[[Outcome]]
  C.med_df <- metadata[, C.med, drop = FALSE]
  C.med_df <- data.frame(C.med_df)
  C.out_df <- metadata[, C.out, drop = FALSE]
  C.out_df <- data.frame(C.out_df)
  M = as.matrix(as.data.frame(t(physeq_level@otu_table)))


  # Design matrix
  if (ncol(C.med_df) == 0) {
    x <- model.matrix(~ Treat)
  } else {
    x <- model.matrix(~ Treat + ., data = C.med_df)
  }

  if (is.null(C.out)){
    X <- model.matrix(~Treat+apply(M,2,scale))
  } else {
    X <- model.matrix(~Treat+apply(M,2,scale)+., data = C.out_df)
  }


  nP = ncol(M)
  nO = nrow(M)
  nC.med = ncol(x)
  nC.out = ncol(X)- nP

  ## Size Factor estimation from real data
  if (Size == "RLE"){
    mu_gm = apply(M,2,function(x) exp(mean(log(x[x > 0]))))
    S = apply(M, 1, function(x) median((x / mu_gm)[x>0]))
  }
  if (Size == "RLEpseudo"){
    mu_gm = apply(M,2,function(x) exp(mean(log(x+ pseudo))))
    S = apply(M, 1, function(x) median((x + pseudo) / mu_gm))
  }
  if (Size == "GMPR"){
    S = GMPR(t(M), intersect.no = 1)
  }
  if (Size == "CSS"){
    S = apply(M, 1, function(x) {
      q = quantile(x[x>0], probs = 0.5)
      sum(x[x <= q])
    })
  }

  ## Dispersion parameter estimation from real data
  tau_prior = apply(M, 2, function(variable){
    non_zero_counts = variable[variable !=0]
    result <- tryCatch({
      non_zero_bn = glm.nb(non_zero_counts ~ x[variable != 0,], control = glm.control(maxit = 5000))
      tau_mean = non_zero_bn[["theta"]]
      tau_var = 0.1
      tau_rate = tau_mean / tau_var
      tau_shape = tau_mean * tau_rate
      c(tau_shape, tau_rate)
    }, error = function(e) {
      tau_mean = 1e5
      tau_var = 0.1
      tau_rate = tau_mean / tau_var
      tau_shape = tau_mean * tau_rate
      c(tau_shape, tau_rate)
    })

    return(result)
  })

  # Prevalence model
  pos_gamma <- array(rep(0, n_iter * nC.med * nP), dim = c(n_iter,nC.med,nP))
  pos_kgamma <- array(rep(0, n_iter * 1 * nP), dim = c(n_iter, 1, nP))

  gamma = pos_gamma[1,,]
  kgamma = pos_kgamma[1,,]

  # Abundance model
  pos_beta <- array(rep(0, n_iter * (nC.med+1) * nP), dim = c(n_iter,nC.med+1,nP))
  pos_kbeta <- array(rep(0, n_iter * 1 * nP), dim = c(n_iter, 1, nP))
  pos_tau <- array(rep(1, n_iter * 2 * nP), dim = c(n_iter, 2, nP))

  beta_current = pos_beta[1,1:nC.med,]
  beta_proposed = beta_current
  kbeta = pos_kbeta[1,,]
  tau_current = rinvgamma(nP,shape = tau_prior[1,],rate = tau_prior[2,])
  tau_proposed = tau_current

  # Outcome model
  pos_alpha <- array(rep(0, n_iter * (nC.out + nP)), dim = c(n_iter, (nC.out + nP), 1))
  pos_kalpha <- array(rep(0, n_iter * 1 * nP), dim = c(n_iter, 1, nP))
  pos_sigmaY <- array(rep(1, n_iter * 1 * 1), dim = c(n_iter, 1, 1))

  alpha = pos_alpha[1,,]
  kalpha = pos_kalpha[1,,]
  sigmaY = pos_sigmaY[1,,]

  ## Structure_Zero_Indicator
  Structure_Zero = array(rep(0,n_iter * nP * nO),dim = c(n_iter, nO, nP))

  ## Variable Selection Parameter
  delta_gamma = rep(1, nP)
  theta_gamma = rep(0.5, nP)

  delta_beta = rep(1, nP)
  theta_beta = rep(0.5, nP)

  delta_alpha = rep(1, nP)
  theta_alpha = rep(0.5, nP)

  for (i in 2:n_iter) {


    for (k in 1:nP) {

      ## Variable selection -------------------------------------------------------------------------

      ## kalpha - inclusion indicator
      nu_alpha = rinvgamma(1, 2 + 0.5, 2 + alpha[2+k]^2 / (2 * delta_alpha[k]))
      a = dnorm(alpha[2+k], 0, sqrt(1 * nu_alpha)) * theta_alpha[k]
      b = dnorm(alpha[2+k], 0, sqrt(0.0001 * nu_alpha)) * (1 - theta_alpha[k])
      kalpha = rbinom(1, 1, a / (a + b))
      pos_kalpha[i,,k] =  kalpha #Update

      ## va0 - prior variance
      delta_alpha[k] = kalpha * 1 + (1 - kalpha) * 0.0001
      va0 = delta_alpha[k] * nu_alpha

      vaM = 1 / (1/va0 + sum(X[,2+k]^2)/ sigmaY)
      saM = X[,2+k] %*% (Y - X %*% alpha + X[,2+k] *alpha[2+k])
      eaM = vaM * saM / sigmaY
      alpha[2+k] = rnorm(1, eaM, sqrt(vaM))

      theta_alpha[k] = rbeta(1, 0.5 +  kalpha, 0.5 + 1 - kalpha)

      ## kgamma - inclusion indicator
      nu_gamma = rinvgamma(1, 2 + 0.5, 2 + gamma[2,k]^2 / (2 * delta_gamma[k]))
      a = dnorm(gamma[2,k], 0, sqrt(1 * nu_gamma)) * theta_gamma[k]
      b = dnorm(gamma[2,k], 0, sqrt(0.1 * nu_gamma)) * (1 - theta_gamma[k])
      kgamma = rbinom(1, 1, a / (a + b))
      pos_kgamma[i,,k] =  kgamma #Update

      ## vg0 - prior variance
      delta_gamma[k] = kgamma * 1 + (1 - kgamma) * 0.1
      vg0 = delta_gamma[k] * nu_gamma

      theta_gamma[k] = rbeta(1, 0.5 +  kgamma, 0.5 + 1 - kgamma)

      ## kbeta - inclusion indicator
      nu_beta = rinvgamma(1, 2 + 0.5, 2 + beta_current[2,k]^2 / (2 * delta_beta[k]))
      a = dnorm(beta_current[2,k], 0, sqrt(1 * nu_beta)) * theta_beta[k]
      b = dnorm(beta_current[2,k], 0, sqrt(0.1 * nu_beta)) * (1 - theta_beta[k])
      kbeta = rbinom(1, 1, a / (a + b))
      pos_kbeta[i,,k] =  kbeta #Update

      ## vb0 - prior variance
      delta_beta[k] = kbeta * 1 + (1 - kbeta) * 0.1
      vb0 =  delta_beta[k] * nu_beta

      theta_beta[k] = rbeta(1, 0.5 +  kbeta, 0.5 + 1 - kbeta)

      #vg0  = 2 ## without variable selection
      #vb0 = 10E4 ## without variable selection

      ## Identify Structural Zeros ----------------------------------------------------------------
      zeros <- (M[,k] == 0)
      c <- 1/(1+exp(- (x[zeros,] %*% gamma[,k])))
      d <- dnbinom(M[zeros,k], size = tau_current[k], mu = S[zeros] * exp(x[zeros,] %*% beta_current[,k]), log = FALSE) * (1-c)
      Omega <- rep(0,nO)
      Omega[zeros] <- rbinom(sum(zeros), 1, c / (c + d))
      Structure_Zero[i,,k] = Omega

      ## Prevalence Model: Polya-gamma Augmentation ----------------------------------------------
      ## gamma --------
      phi = rpg(nO,1,x %*% gamma[,k])
      Sigma_g = solve(diag(1/c(1, vg0,rep(10E4,nC.med-2))) + t(x) %*% diag(phi) %*% x)
      mu_g = Sigma_g %*% t(x) %*% (Omega-1/2)
      gamma[,k] <- mvrnorm(1, mu_g, Sigma_g)
      pos_gamma[i,,k] = gamma[,k] #Update

      ## Abundance Model: random walk MH --------------------------------------------------------
      M_NB = M[Omega == 0,k]
      x_NB = x[Omega == 0,]
      S_NB = S[Omega == 0]

      ## beta ---------
      beta_proposed[,k] <- mvrnorm(1, beta_current[,k], diag(0.01,nrow = nC.med))

      # Calculate log-likelihoods
      loglik_current <-  sum(dnbinom(M_NB, size = tau_current[k], mu = S_NB * exp(x_NB %*% beta_current[,k]), log = TRUE)) +
        sum(dnorm(beta_current[-2,k], 0, 10E4, log = TRUE)) +
        sum(dnorm(beta_current[2,k], 0, sqrt(vb0), log = TRUE)) +
        dgamma(tau_current[k], shape = tau_prior[1,k], rate = tau_prior[2,k], log = TRUE)

      loglik_proposed <-  sum(dnbinom(M_NB, size = tau_current[k], mu = S_NB * exp(x_NB %*% beta_proposed[,k]), log = TRUE)) +
        sum(dnorm(beta_proposed[-2,k], 0, 10E4, log = TRUE)) +
        sum(dnorm(beta_proposed[2,k], 0, sqrt(vb0), log = TRUE)) +
        dgamma(tau_current[k], shape = tau_prior[1,k], rate = tau_prior[2,k], log = TRUE)

      AR <- exp(loglik_proposed - loglik_current) # acceptance probability

      # Accept or reject
      if (runif(1) < AR) {
        beta_current[,k] <- beta_proposed[,k]
      }
      pos_beta[i, 1:ncol(x),k] <- beta_current[,k] #Update
      pos_beta[i, ncol(x)+1,k] = min(1,AR) #Update

      ## tau --------
      repeat {
        tau_proposed[k] = rnorm(1, tau_current, sd = 0.1)
        if (tau_proposed[k] >= 0) break
      }

      # Calculate log-likelihoods
      loglik_current <- sum(dnbinom(M_NB, size = tau_current[k], mu = S_NB * exp(x_NB %*% beta_current[,k]), log = TRUE)) +
        sum(dnorm(beta_current[-2,k], 0, 10E4, log = TRUE)) +
        sum(dnorm(beta_current[2,k], 0, sqrt(vb0), log = TRUE)) +
        dgamma(tau_current[k], shape = tau_prior[1,k], rate = tau_prior[2,k], log = TRUE)

      loglik_proposed <- sum(dnbinom(M_NB, size = tau_proposed[k], mu = S_NB * exp(x_NB %*% beta_current[,k]), log = TRUE)) +
        sum(dnorm(beta_current[-2,k], 0, 10E4, log = TRUE)) +
        sum(dnorm(beta_current[2,k], 0, sqrt(vb0), log = TRUE)) +
        dgamma(tau_proposed[k], shape = tau_prior[1,k], rate = tau_prior[2,k], log = TRUE)

      AR <- exp(loglik_proposed - loglik_current) # acceptance probability

      # Accept or reject
      if (runif(1) < AR) {
        tau_current[k] <- tau_proposed[k]
      }
      pos_tau[i, 1,k] <- tau_current[k]
      pos_tau[i, 2,k] = min(1,AR)
    }

    ## Outcome model: regular Gibbs -----------------------------------------------------------------
    alpha[1] = rnorm(1, sum(Y - X %*% alpha + alpha[1] * X[,1])/nO, sqrt(sigmaY/nO))
    alpha[2] = rnorm(1, X[,2] %*% (Y - X %*% alpha +  alpha[2] * X[,2])/sum(X[,2]^2), sqrt(sigmaY/sum(X[,2]^2)))
    if (nC.out > 2) {
      for (c in 1:(nC.out-2)) {
        alpha[2+nP+c] = rnorm(1, X[,(2+nP+c)] %*% (Y - X %*% alpha  + alpha[2+nP+c] * X[,(2+nP+c)])/sum(X[,(2+nP+c)]^2), sqrt(sigmaY/sum(X[,(2+nP+c)]^2)))
      }
    }
    YSSR = sum((Y - X %*% alpha)^2)
    sigmaY = rinvgamma(1, nO/2 + 10E-4,  YSSR/2 + 10E-4)

    pos_sigmaY[i,,] =  sigmaY #Update
    pos_alpha[i,,] = alpha #Update

    gc()

  }
  posterior_list = list(pos_beta = pos_beta, pos_gamma = pos_gamma,pos_alpha = pos_alpha,
                        pos_kbeta = pos_kbeta,pos_kgamma = pos_kgamma,pos_kalpha = pos_kalpha,
                        pos_tau = pos_tau, pos_sigmaY = pos_sigmaY,structure_Zero = Structure_Zero)

  posterior_list <- lapply(posterior_list, function(array) {
    array[-(1:burn_in),,]
  })

  posterior_list$prior_tau = tau_prior
  posterior_list$M.level = M.level
  posterior_list$object = physeq_level
  return(posterior_list)
}




#' Geometric Mean of Pairwise Ratios (GMPR) Normalization
#'
#' This function calculates normalizing factors for microbiome sequencing data (or more generally, zero-inflated sequencing data)
#' using the GMPR method. The resulting size factors can be used as offsets in count-based regression models or to produce normalized data.
#'
#' @param comm A matrix or data frame of raw count data (features x samples), where rows are taxa/features and columns are samples.
#' @param intersect.no Integer. Minimum number of shared features required for calculating pairwise ratios (default = 1).
#' @param ct.min Integer. Minimum count threshold to consider a feature as present (default = 1).
#' @param trace Logical. Whether to print progress messages during computation (default = TRUE).
#'
#' @return A numeric vector of sample-specific size factors, with names corresponding to the sample IDs.
#'
#' @import matrixStats
#' @export
GMPR <- function (comm, intersect.no = 1, ct.min = 1, trace = TRUE) {
  # Computes the GMPR size factor
  #
  # Args:
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
  #   ct.min: the minimum number of counts required to calculate ratios

  #
  # Returns:
  #   a vector of the size factors with attribute 'NSS'. Samples with distinct sets of features will be output as NA.
  #         NSS:   number of samples with significant sharing (> intersect.no) including itself

  # mask counts < ct.min
  comm[comm < ct.min] <- 0

  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }

  if (trace) cat('Begin GMPR size factor calculation ...\n')

  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {
    if (i %% 50 == 0) {
      cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )

  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'),
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n',
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }

  if (trace) cat('Completed!\n')
  if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  names(gmpr) <- names(comm.no) <- colnames(comm)

  attr(gmpr, 'NSS') <- comm.no

  return(gmpr)
}

