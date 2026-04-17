# DEFINE EACH OF 5 MODELS FOR KNOCKOUT EXPERIMENT
# 12/13/2024: removed covars on pc detection for now
# 12/16/2024: toggling on PPCs

# FULL MODEL [H-A-S] --------------------------------------------------------------
# includes point count data, binomial ARU data from thresholded ML scores, ML scores

base_dir <- "COMB_minimal/models/src/paper_experiments/jags_files/no_covars/"

if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE) 
} else {
  print("Directory already exists.")
}

cat("
model {

  # Priors
  psi ~ dunif(0, 1)
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) 
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, .1)
  mu[2] ~ dnorm(-2, .1)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # Latent occupancy states

    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:n_v[i]) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }

    # GOF Point Count - Tukey-Freeman Discrepancy
    T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
    y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
    T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data

  
    # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:n_s[i]) {
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

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites
}
    ", file = paste0(base_dir, "model_HAS.txt"))

# HUMANS ONLY [H] -----------------------------------------------------------------

cat("
model {

  # Priors
  psi ~ dunif(0, 1)
  p11 ~ dbeta(2, 2)  

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) 

    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:n_v[i]) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
   
    # GOF Point Count - Tukey-Freeman Discrepancy
    T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
    y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
    T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
      
  } # end of site loop i


  # GOF assessment
  T_pc_obs <- sum(T_pc_obs0)
  T_pc_sim <- sum(T_pc_sim0)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites

} # end of model loop
    ", file = paste0(base_dir, "model_H.txt"))


# ARU DETECTIONS + SCORES [A-S]----------------------------------------------------------------
cat("
model {
  
  # Priors
  psi ~ dunif(0, 1)
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, 5)
  mu[2] ~ dnorm(-2, 5)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])
  
  # Likelihood part 1 & 2: PC and ARU detections
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # Latent occupancy states
    
    # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:n_s[i]) {
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

  mean_psi <- mean(psi)
  NOcc <- sum(z[])
  mu_diff_2_1 <- mu[2] - mu[1]
} 
", file = paste0(base_dir, "model_AS.txt"))


# HUMANS + ARU DETECTIONS [H-A]---------------------------------------------------------------

cat("
model {
  # Priors
  psi ~ dunif(0,1)
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)  
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)

  # Likelihood part 1 & 2: PC and ARU detections
  for (i in 1:nsites) { # Loop over sites

    z[i] ~ dbern(psi) # Latent occupancy states

    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:n_v[i]) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
   
    # GOF Point Count - Tukey-Freeman Discrepancy
    T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
    y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
    T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data

    # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 
    for(j in 1:n_s[i]) {
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

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites
  
} # end of model loop
", file = paste0(base_dir, "model_HA.txt"))


# HUMANS + SCORES [H-S]-------------------------------------------------------------
cat("
model {

  # Priors
  psi ~ dunif(0, 1)
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) # commenting out p11 bc we are providing parameters for p below 

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, 5)
  mu[2] ~ dnorm(-2, 5)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites

    z[i] ~ dbern(psi) # Latent occupancy states

    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:n_v[i]) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
   
    # GOF Point Count - Tukey-Freeman Discrepancy
    T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
    y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
    T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
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

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites
  mu_diff_2_1 <- mu[2] - mu[1]

} # end of model loop
", file = paste0(base_dir, "model_HS.txt"))


# SCORES ONLY [S] ---------------------------------------------------------

cat("
model {

  # Priors
  psi ~ dunif(0,1)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, .1)
  mu[2] ~ dnorm(-2, .1)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1: No ARU or PC
  
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # Latent occupancy states
  }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    #GOF for ARU scores
      LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }


  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites
  mu_diff_2_1 <- mu[2] - mu[1]
}
    ", file = paste0(base_dir, "model_S.txt"))


# ARU DETECTIONS ONLY [A] ------------------------------------------------

cat("
model {

  # Priors
  psi ~ dunif(0,1)
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # Latent occupancy states
    
    # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:n_s[i]) {
      y.aru[i,j] ~ dbern(p_aru[i])  
    }

    # GOF ARU Count - Tukey-Freeman Discrepancy
    T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
    y_aru_Sim[i] ~ dbin(p_aru[i], n_s[i])
    T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
  }

  T_aru_obs <- sum(T_aru_obs0)
  T_aru_sim <- sum(T_aru_sim0)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites
}
    ", file = paste0(base_dir, "model_A.txt"))

