# DEFINE EACH OF 5 MODELS FOR KNOCKOUT EXPERIMENT
# 12/13/2024: removed covars on pc detection for now
# 12/16/2024: toggling on PPCs

# FULL MODEL [H-A-S] --------------------------------------------------------------
# includes point count data, binomial ARU data from thresholded ML scores, ML scores

cat("
model {

  # Priors
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) 
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  
#  alpha0 ~ dunif(-5, 5)
#  alpha1 ~ dunif(-5, 5)
#  alpha2 ~ dunif(-5, 5)

  beta0 ~ dunif(-5, 5) # Occupancy intercept on prob. scale
  beta1 ~ dunif(-5, 5)
#  beta2 ~ dunif(-5, 5)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, .1)
  mu[2] ~ dnorm(-2, .1)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    
    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*burn[i] )/(1+exp(beta0 + beta1*burn[i] ))
    zSim[i] ~ dbern(psiSim[i])


    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }

    # GOF Point Count - Tukey-Freeman Discrepancy
      T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
      y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
      T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
  
 # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dbern(p_aru[i])  # n_surveys.aru = Total files processed
      }

    # GOF ARU Count - Tukey-Freeman Discrepancy
    T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
    y_aru_Sim[i] ~ dbin(p_aru[i], n_s[i])
    T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
  }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    #GOF for ARU scores
      LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }

  # GOF assessment
  T_pc_obs <- sum(T_pc_obs0)
  T_pc_sim <- sum(T_pc_sim0)
  T_aru_obs <- sum(T_aru_obs0)
  T_aru_sim <- sum(T_aru_sim0)

  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)

  #GOF - Regression: Difference in mean vegetation
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TvegObs <- sum(burn*z) / ifelse(cz1>0,cz1, 1)  - sum(burn*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(burn*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(burn*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])

  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta1 * Xburn[k] # psi predictions for burn over the range of burn
  }

}
    ", file = "COMB_minimal/models/jags_files/model_HAS.txt")

# HUMANS ONLY [H] -----------------------------------------------------------------

cat("
model {

  # Priors
  p11 ~ dbeta(2, 2)  # comment out when modeling pc detection w/ covariates
#  alpha0 ~ dunif(-5, 5)
#  alpha1 ~ dunif(-5, 5)
#  alpha2 ~ dunif(-5, 5)

  beta0 ~ dunif(-5, 5) # Occupancy intercept on prob. scale
  beta1 ~ dunif(-5, 5)
#  beta2 ~ dunif(-5, 5)

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*burn[i]
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    
    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*burn[i] )/(1+exp(beta0 + beta1*burn[i] ))
    zSim[i] ~ dbern(psiSim[i])

        # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
   
  #  for(j in 1:nsurveys.pc) {
  #    y.ind[i,j] ~ dbern(p11[i,j]*z[i]) # Observed occ. data (if available)
  #    logit(p11[i,j]) <- alpha0 + alpha1*Time[i,j] + alpha2*Date[i,j]
    
    # GOF Point Count - Tukey-Freeman Discrepancy (for when p11 is parameterized)
  #  expected_pc[i,j] = p11[i,j]*z[i] # in j loop when y.ind is parameterized
  #  y_pc_Sim[i,j] ~ dbern(p11[i,j] * z[i])
    
 # }
        # GOF Point Count - Tukey-Freeman Discrepancy
      T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
      y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
      T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
      
  #  expected_sum_pc[i] = sum(expected_pc[i, ])
  #  T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(expected_sum_pc[i]))^2 
  #  y_pc_Sim_sum[i] = sum(y_pc_Sim[i,])
  #  T_pc_sim0[i] <- (sqrt(y_pc_Sim_sum[i]) - sqrt(expected_sum_pc[i]))^2  # ...and for simulated data
    

  } # end of site loop i


  # GOF assessment
  T_pc_obs <- sum(T_pc_obs0)
  T_pc_sim <- sum(T_pc_sim0)

  #GOF - Regression: Difference in mean vegetation
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TvegObs <- sum(burn*z) / ifelse(cz1>0,cz1, 1)  - sum(burn*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(burn*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(burn*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])

  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta1 * Xburn[k] # psi predictions for burn with canopy cover held at mean
  }

} # end of model loop
    ", file = "COMB_minimal/models/jags_files/model_H.txt")


# ARU DETECTIONS + SCORES [A-S]----------------------------------------------------------------
cat("
model {
  
  # Priors
  #p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) 
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  
  #  alpha0 ~ dunif(-5, 5)
  #  alpha1 ~ dunif(-5, 5)
  #  alpha2 ~ dunif(-5, 5)

  beta0 ~ dunif(-5, 5) # Occupancy intercept on prob. scale
  beta1 ~ dunif(-5, 5)
  #  beta2 ~ dunif(-5, 5)
  
  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, 5)
  mu[2] ~ dnorm(-2, 5)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])
  
  # Likelihood part 1 & 2: PC and ARU detections
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    
    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*burn[i] )/(1+exp(beta0 + beta1*burn[i] ))
    zSim[i] ~ dbern(psiSim[i])
    
    # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dbern(p_aru[i]) 
    }
    
    # GOF ARU Count - Tukey-Freeman Discrepancy
    T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
    y_aru_Sim[i] ~ dbin(p_aru[i], n_v[i])
    T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
  } # END of site loop
  
  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    #GOF for ARU scores
    LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }
  
  # GOF assessment
  T_aru_obs <- sum(T_aru_obs0)
  T_aru_sim <- sum(T_aru_sim0)
  
  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)
  
  #GOF - Regression: Difference in mean env. variable
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TvegObs <- sum(burn*z) / ifelse(cz1>0,cz1, 1)  - sum(burn*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(burn*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(burn*(1-zSim))/ifelse(czSim0>0,czSim0, 1)
  
  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  
  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta1 * Xburn[k] # psi predictions for burn
  }
  
} # end of model loop
", file = "COMB_minimal/models/jags_files/model_AS.txt")


# HUMANS + ARU DETECTIONS [H-A]---------------------------------------------------------------

cat("
model {
  # Priors
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) # commenting out p11 bc we are providing parameters for p below 
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  
#  alpha0 ~ dunif(-5, 5)
#  alpha1 ~ dunif(-5, 5)
#  alpha2 ~ dunif(-5, 5)

  beta0 ~ dunif(-5, 5) # Occupancy intercept on prob. scale
  beta1 ~ dunif(-5, 5)
#  beta2 ~ dunif(-5, 5)


  # Likelihood part 1 & 2: PC and ARU detections
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    
    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*burn[i] )/(1+exp(beta0 + beta1*burn[i] ))
    zSim[i] ~ dbern(psiSim[i])
    
    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
   
  #  for(j in 1:nsurveys.pc) {
  #    y.ind[i,j] ~ dbern(p11[i,j]*z[i]) # Observed occ. data (if available)
  #    logit(p11[i,j]) <- alpha0 + alpha1*Time[i,j] + alpha2*Date[i,j]
    
    # GOF Point Count - Tukey-Freeman Discrepancy (for when p11 is parameterized)
  #  expected_pc[i,j] = p11[i,j]*z[i] # in j loop when y.ind is parameterized
  #  y_pc_Sim[i,j] ~ dbern(p11[i,j] * z[i])
    
 # }
        # GOF Point Count - Tukey-Freeman Discrepancy
      T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
      y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
      T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
      
  #  expected_sum_pc[i] = sum(expected_pc[i, ])
  #  T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(expected_sum_pc[i]))^2 
  #  y_pc_Sim_sum[i] = sum(y_pc_Sim[i,])
  #  T_pc_sim0[i] <- (sqrt(y_pc_Sim_sum[i]) - sqrt(expected_sum_pc[i]))^2  # ...and for simulated data


    # ARU - binomial
      p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
      for(j in 1:nsurveys.aru) {
        y.aru[i,j] ~ dbern(p_aru[i]) 
    }

    # GOF ARU Count - Tukey-Freeman Discrepancy
    T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
    y_aru_Sim[i] ~ dbin(p_aru[i], n_v[i])
    T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
  
  } # end of site loop i
  # GOF assessment
  T_pc_obs <- sum(T_pc_obs0)
  T_pc_sim <- sum(T_pc_sim0)
  T_aru_obs <- sum(T_aru_obs0)
  T_aru_sim <- sum(T_aru_sim0)

  #GOF - Regression: Difference in mean env. variable
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TvegObs <- sum(burn*z) / ifelse(cz1>0,cz1, 1)  - sum(burn*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(burn*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(burn*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])

  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta1 * Xburn[k] # psi predictions for burn with canopy cover held at mean
  }
} # end of model loop
", file = "COMB_minimal/models/jags_files/model_HA.txt")


# HUMANS + SCORES [H-S]-------------------------------------------------------------
cat("
model {

  # Priors
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) # commenting out p11 bc we are providing parameters for p below 
#  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
#  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  
#  alpha0 ~ dunif(-5, 5)
#  alpha1 ~ dunif(-5, 5)
#  alpha2 ~ dunif(-5, 5)

  beta0 ~ dunif(-5, 5) # Occupancy intercept on prob. scale
  beta1 ~ dunif(-5, 5)
#  beta2 ~ dunif(-5, 5)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, 5)
  mu[2] ~ dnorm(-2, 5)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    
    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*burn[i] )/(1+exp(beta0 + beta1*burn[i] ))
    zSim[i] ~ dbern(psiSim[i])

    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
   
  #  for(j in 1:nsurveys.pc) {
  #    y.ind[i,j] ~ dbern(p11[i,j]*z[i]) # Observed occ. data (if available)
  #    logit(p11[i,j]) <- alpha0 + alpha1*Time[i,j] + alpha2*Date[i,j]
    
    # GOF Point Count - Tukey-Freeman Discrepancy (for when p11 is parameterized)
  #  expected_pc[i,j] = p11[i,j]*z[i] # in j loop when y.ind is parameterized
  #  y_pc_Sim[i,j] ~ dbern(p11[i,j] * z[i])
    
 # }
        # GOF Point Count - Tukey-Freeman Discrepancy
      T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
      y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
      T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
      
  #  expected_sum_pc[i] = sum(expected_pc[i, ])
  #  T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(expected_sum_pc[i]))^2 
  #  y_pc_Sim_sum[i] = sum(y_pc_Sim[i,])
  #  T_pc_sim0[i] <- (sqrt(y_pc_Sim_sum[i]) - sqrt(expected_sum_pc[i]))^2  # ...and for simulated data

  }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    #GOF for ARU scores
      LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }

  # GOF assessment
  T_pc_obs <- sum(T_pc_obs0)
  T_pc_sim <- sum(T_pc_sim0)

  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)

  #GOF - Regression: Difference in mean vegetation
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TvegObs <- sum(burn*z) / ifelse(cz1>0,cz1, 1)  - sum(burn*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(burn*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(burn*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])

  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta1 * Xburn[k] # psi predictions for burn with canopy cover held at mean
  }
} # end of model loop
", file = "COMB_minimal/models/jags_files/model_HS.txt")


# SCORES ONLY [S] ---------------------------------------------------------

cat("
model {

  # Priors
#  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) 
#  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
#  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  
#  alpha0 ~ dunif(-5, 5)
#  alpha1 ~ dunif(-5, 5)
#  alpha2 ~ dunif(-5, 5)

  beta0 ~ dunif(-5, 5) # Occupancy intercept on prob. scale
  beta1 ~ dunif(-5, 5)
#  beta2 ~ dunif(-5, 5)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, .1)
  mu[2] ~ dnorm(-2, .1)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    
    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*burn[i] )/(1+exp(beta0 + beta1*burn[i] ))
    zSim[i] ~ dbern(psiSim[i])
   }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    #GOF for ARU scores
      LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }

  # GOF assessment
  # T_pc_obs <- sum(T_pc_obs0)
  # T_pc_sim <- sum(T_pc_sim0)
  # T_aru_obs <- sum(T_aru_obs0)
  # T_aru_sim <- sum(T_aru_sim0)

  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)

  #GOF - Regression: Difference in mean vegetation
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TvegObs <- sum(burn*z) / ifelse(cz1>0,cz1, 1)  - sum(burn*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(burn*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(burn*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])

  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta1 * Xburn[k] # psi predictions for burn over the range of burn
  }
}
    ", file = "COMB_minimal/models/jags_files/model_S.txt")


# ARU DETECTIONS ONLY [A] ------------------------------------------------

cat("
model {

  # Priors
#  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) 
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  
#  alpha0 ~ dunif(-5, 5)
#  alpha1 ~ dunif(-5, 5)
#  alpha2 ~ dunif(-5, 5)

  beta0 ~ dunif(-5, 5) # Occupancy intercept on prob. scale
  beta1 ~ dunif(-5, 5)
#  beta2 ~ dunif(-5, 5)

  # Parameters of the observation model for the scores
  # mu[1] ~ dnorm(-2, .1)
  # mu[2] ~ dnorm(-2, .1)
  # sigma[1] ~ dunif(0.1, 5)
  # tau[1] <- 1 / (sigma[1] * sigma[1])
  # sigma[2] ~ dunif(0.1, 5)
  # tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    
    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*burn[i] )/(1+exp(beta0 + beta1*burn[i] ))
    zSim[i] ~ dbern(psiSim[i])

#     # Point count detection process
#     p[i] <- z[i]*p11 # Detection probability
#     for(j in 1:nsurveys.pc) {
#       y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
#     }
# 
#     # GOF Point Count - Tukey-Freeman Discrepancy
#       T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
#       y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
#       T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
  
 # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dbern(p_aru[i])  # n_surveys.aru = Total files processed
      }

    # GOF ARU Count - Tukey-Freeman Discrepancy
    T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
    y_aru_Sim[i] ~ dbin(p_aru[i], n_s[i])
    T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
  }

  # # Likelihood part 2: feature score data
  # for(k in 1:nsamples) {
  #   score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  #   #GOF for ARU scores
  #     LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  #     score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  #     LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  # }

  # GOF assessment
  # T_pc_obs <- sum(T_pc_obs0)
  # T_pc_sim <- sum(T_pc_sim0)
  T_aru_obs <- sum(T_aru_obs0)
  T_aru_sim <- sum(T_aru_sim0)

  # #GOF assessment for scores
  # Dobs <- -2 * sum(LLobs_score)
  # Dsim <- -2 * sum(LLsim_score)

  #GOF - Regression: Difference in mean vegetation
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TvegObs <- sum(burn*z) / ifelse(cz1>0,cz1, 1)  - sum(burn*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(burn*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(burn*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])

  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta1 * Xburn[k] # psi predictions for burn over the range of burn
  }

}
    ", file = "COMB_minimal/models/jags_files/model_A.txt")
