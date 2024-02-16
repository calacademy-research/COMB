#library(R2jags)
library(runjags)
library(mcmcplots)
library(jagsUI)
library(tidyverse)

gen_data_text <- function(
    include_aru = TRUE,
    include_pc = TRUE,
    include_scores = TRUE,
    include_covar = TRUE) {
  string_init <-
    "
data {"
  
  if (include_covar) {
    occ_lik <- " 
  # Occupancy
  for (i in 1:nsites) { 
    logit(psi[i]) <- beta0 + beta1*covar[i]
    z[i] ~ dbern(psi[i])"
  } else {
    occ_lik <- "
  # Occupancy
    for (i in 1:nsites) { 
      z[i] ~ dbern(psi)" 
  }
  
  if (include_pc) {
    pc_lik <- "
    # Point count 
    p[i] <- z[i]*p11 
    for(j in 1:nsurveys.pc) {
      y.pc[i,j] ~ dbern(p[i]) 
    }"
  } else {
    pc_lik <- ""
  }
  
  if (include_aru & include_scores) {
    aru_lik <- "
    # ARU: detection + scores 
    p_aru[i] <- z[i]*p_aru11 + p_aru01  
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dbern(p_aru[i]) 
      score[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1]) T(ifelse(y.aru[i,j] == 1, threshold, -10), ifelse(y.aru[i,j] == 1, 10, threshold))
    }
  }"
  } else if (include_aru) {
    aru_lik <- "    # ARU - binomial
       p_aru[i] <- z[i]*p_aru11 + p_aru01
       for(j in 1:nsurveys.aru) {
         y.aru[i,j] ~ dbern(p_aru[i])
       }
     }"
  } else if (include_scores) {
    aru_lik <- "
       for(j in 1:nsurveys.aru) {
         score[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1])
       }
     }"
  } else {
    aru_lik <- "
   	 }"
  }
  
  close_string <- "
}
model{
	fake <- 0
}
"
data_mod_string <- paste(string_init, occ_lik, pc_lik, aru_lik, close_string, sep = "\n")
return(data_mod_string)
}




sim_input_data <- function(
    nsites = 100,
    p11 = 0.1,
    p_aru11 = 0.1,
    p_aru01 = 0.025,
    beta0 = 0,
    beta1 = 1,
    mu = c(-1, 1),
    sigma = c(1, 1),
    nsurveys.aru = 24,
    nsurveys.pc = 3,
    covar_prob = 0.5,
    threshold = 0,
    covar_continuous = FALSE) {
  tau <- 1 / (sigma * sigma)
  siteid <- rep(1:nsites, each = 1)
  nsamples <- length(siteid)
  if (covar_continuous) {
    covar <- rnorm(n = nsites)
  } else {
    covar <- rbinom(n = nsites, size = 1, prob = covar_prob)
  }
  
  input_data <- list(
    nsites = nsites,
    p11 = p11,
    p_aru11 = p_aru11,
    p_aru01 = p_aru01,
    beta0 = beta0,
    beta1 = beta1,
    mu = mu,
    tau = tau,
    siteid = siteid,
    nsamples = nsamples,
    nsurveys.aru = nsurveys.aru,
    nsurveys.pc = nsurveys.pc,
    covar = covar,
    threshold = threshold
  )
  return(input_data)
}


gen_simulated_data <- function(
    include_aru = TRUE,
    include_pc = TRUE,
    include_scores = TRUE,
    include_covar = TRUE,
    nsites = 100,
    p11 = 0.1,
    p_aru11 = 0.1,
    p_aru01 = 0.025,
    beta0 = 0,
    beta1 = 1,
    mu = c(-1, 1),
    sigma = c(1, 1),
    nsurveys.aru = 24,
    nsurveys.pc = 3,
    covar_prob = 0.5,
    threshold = 0,
    covar_continuous = FALSE,
    return_inputs = TRUE) {
  data_gen_string <- gen_data_text(
    include_aru,
    include_pc,
    include_scores,
    include_covar
  )
  
  input_data <- sim_input_data(
    nsites = nsites,
    p11 = p11,
    p_aru11 = p_aru11,
    p_aru01 = p_aru01,
    beta0 = beta0,
    beta1 = beta1,
    mu = mu,
    sigma = sigma,
    nsurveys.aru = nsurveys.aru,
    nsurveys.pc = nsurveys.pc,
    covar_prob = covar_prob,
    threshold = threshold,
    covar_continuous = covar_continuous
  )
  
  out <- run.jags(data_gen_string, data = input_data, monitor = c("y.aru", "y.pc", "score"), sample = 1, n.chains = 1, summarise = FALSE)
  sim_data <- out$mcmc[[1]]
  y.aru.sim <- matrix(sim_data[1:(nsurveys.aru * nsites)], nrow = nsites)
  y.pc.sim <- matrix(sim_data[(1 + (nsurveys.aru * nsites)):((nsurveys.aru * nsites) + (nsurveys.pc * nsites))], nrow = nsites)
  score.sim <- matrix(sim_data[(1 + (nsurveys.aru * nsites) + (nsurveys.pc * nsites)):length(sim_data)], nrow = nsites)
  
  if (return_inputs) {
    sim_data_out <- list(
      nsites = nsites,
      siteid = input_data$siteid,
      nsamples = input_data$nsamples,
      nsurveys.aru = nsurveys.aru,
      nsurveys.pc = nsurveys.pc,
      covar = input_data$covar,
      threshold = threshold,
      y.aru = y.aru.sim,
      y.pc = y.pc.sim,
      score = score.sim
    )
  } else {
    sim_data_out <- list(
      y.aru = y.aru.sim,
      y.pc = y.pc.sim,
      score = score.sim
    )
  }
  return(sim_data_out)
}


gen_jags_model_text <- function(
    include_aru_model = TRUE,
    include_pc_model = TRUE,
    include_scores_model = TRUE,
    include_covar_model = TRUE,
    psi_prior = "dbeta(2,2)",
    beta0_prior = "dnorm(0, 1)",
    beta1_prior = "dnorm(0, 1)",
    p11_prior = "dbeta(2, 2)",
    mu1_prior = "dnorm(-2, 0.2)",
    mu2_prior = "dnorm(-2, 0.2)",
    sigma1_prior = "dunif(0.1, 5)",
    sigma2_prior = "dunif(0.1, 5)",
    p_aru11_prior = "dbeta(2, 2)",
    p_aru01_prior = "dbeta(1, 3)I(0, 1 - p_aru11)") {
  string_init <-
    "
  model {
  "
  priors_string <- "# Priors
  "
  
  if (sum(include_pc_model, include_aru_model, include_scores_model) == 0) stop("Must include at least one source of observed data")
  
  if (include_covar_model) {
    reg_priors <- paste0("beta0 ~ ", beta0_prior, "\n beta1 ~ ", beta1_prior)
    priors_string <- paste(priors_string, reg_priors, sep = "\n")
    
    occ_lik <- " 
  	# Likelihood for occupancy
    for (i in 1:nsites) {
      logit(psi[i]) <- beta0 + beta1*covar[i]
      z[i] ~ dbern(psi[i])
    "
  } else {
    priors_string <- paste(priors_string, "psi ~ ", psi_prior, sep = "\n")
    occ_lik <- "
    # Likelihood for occupancy
    for (i in 1:nsites) { # Loop over sites
      z[i] ~ dbern(psi) # Latent occupancy states
    "
  }
  
  if (include_pc_model) {
    pc_lik <- "
      # Point count
      p[i] <- z[i]*p11
      for(j in 1:nsurveys.pc) {
        y.pc[i,j] ~ dbern(p[i]) 
      }"
    p11_prior_string <- paste0("p11 ~ ", p11_prior)
    priors_string <- paste(priors_string, p11_prior_string, sep = "\n")
  } else {
    pc_lik <- ""
  }
  
  if (include_aru_model & include_scores_model) {
    aru_lik <- " 
      # ARU 
      p_aru[i] <- z[i]*p_aru11 + p_aru01 
      for(j in 1:nsurveys.aru) {
        y.aru[i,j] ~ dbern(p_aru[i])  
        score[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1]) T(ifelse(y.aru[i,j] == 1, threshold, -10), ifelse(y.aru[i,j] == 1, 10, threshold))
      }
    }"
    p_aru11_prior_string <- paste0("p_aru11 ~ ", p_aru11_prior)
    p_aru01_prior_string <- paste0("p_aru01 ~ ", p_aru01_prior)
    scores_prior <- paste("# Priors of the observation model for the scores",
                          paste0("mu[1] ~ ", mu1_prior),
                          paste0("mu[2] ~", mu2_prior),
                          paste0("sigma[1] ~ ", sigma1_prior),
                          "tau[1] <- 1 / (sigma[1] * sigma[1])",
                          paste0("sigma[2] ~ ", sigma2_prior),
                          "tau[2] <- 1 / (sigma[2] * sigma[2])",
                          sep = "\n"
    )
    priors_string <- paste(priors_string, p_aru11_prior_string, p_aru01_prior_string, scores_prior, sep = "\n")
  } else if (include_aru_model) {
    aru_lik <- "    # ARU - binomial
       p_aru[i] <- z[i]*p_aru11 + p_aru01
       for(j in 1:nsurveys.aru) {
         y.aru[i,j] ~ dbern(p_aru[i])
        }
     }"
    p_aru11_prior_string <- paste0("p_aru11 ~ ", p_aru11_prior)
    p_aru01_prior_string <- paste0("p_aru01 ~ ", p_aru01_prior)
    priors_string <- paste(priors_string, p_aru11_prior_string, p_aru01_prior_string, sep = "\n")
  } else if (include_scores_model) {
    aru_lik <- "     
       for(j in 1:nsurveys.aru) {
         score[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1])
       }
     }"
    scores_prior <- paste("# Priors of the observation model for the scores",
                          paste0("mu[1] ~ ", mu1_prior),
                          paste0("mu[2] ~", mu2_prior),
                          paste0("sigma[1] ~ ", sigma1_prior),
                          "tau[1] <- 1 / (sigma[1] * sigma[1])",
                          paste0("sigma[2] ~ ", sigma2_prior),
                          "tau[2] <- 1 / (sigma[2] * sigma[2])",
                          sep = "\n"
    )
    priors_string <- paste(priors_string, scores_prior, sep = "\n")
  } else {
    aru_lik <- "
     }"
  }
  
  close_string <- "
  mean_psi <- mean(psi)
  NOcc <- sum(z[]) 
  PropOcc <- NOcc/nsites
}
"
data_mod_string <- paste(string_init, priors_string, occ_lik, pc_lik, aru_lik, close_string, sep = "\n")
return(data_mod_string)
}

gen_monitored_param_list <- function(
    include_aru_model = TRUE,
    include_pc_model = TRUE,
    include_scores_model = TRUE,
    include_covar_model = TRUE) {
  monitored_list <- c("deviance", "PropOcc")
  
  if (include_covar_model) {
    monitored_list <- append(monitored_list, c("beta0", "beta1", "beta2", "mean_psi"))
  } else {
    monitored_list <- append(monitored_list, "psi")
  }
  
  if (include_pc_model) {
    monitored_list <- append(monitored_list, "p11")
  }
  
  if (include_aru_model) {
    monitored_list <- append(monitored_list, c("p_aru11", "p_aru01"))
  }
  
  if (include_scores_model) {
    monitored_list <- append(monitored_list, c("mu", "sigma"))
  }
  return(monitored_list)
}

fit_model_to_sim_data <- function(
    jagsInputData = NULL,
    jags_model_fit_text = NULL,
    monitored_list = NULL,
    nsites = 100,
    na = 1000,
    ni = 5000,
    nt = 1,
    nb = 1000,
    nc = 6) {
  # if(is.na(jagsInputData)) stop('Input dataset required')
  # if(is.na(jags_model_fit_text)) stop('No model text specified')
  # if(is.na(monitored_list)) stop('Need to include list of monitored parameters')
  
  modelFile <- tempfile()
  cat(file = modelFile, jags_model_fit_text)
  
  inits <- function() {
    list(
      mu = sort(rnorm(2)),
      sigma = runif(2, 0.5, 2.5),
      z = rep(1, nsites),
      beta0 = 0,
      beta1 = 0
    )
  }
  
  simjagsResult <- jagsUI::jags(
    jagsInputData,
    inits,
    monitored_list,
    modelFile,
    n.adapt = na,
    n.chains = nc,
    n.thin = nt,
    n.iter = ni,
    n.burnin = nb,
    parallel = TRUE
  )
  return(simjagsResult)
}


expand_param_list <- function(param_list) {
  param_names <- names(param_list)
  list_out <- list()
  for (param_name in param_names) {
    param <- param_list[[param_name]]
    if (length(param) == 1) {
      list_out[param_name] <- param
    } else {
      for (i in 1:length(param)) {
        name <- paste0(param_name, "[", i, "]")
        list_out[name] <- param[i]
      }
    }
  }
  return(list_out)
}


run_n_simulation_and_fit <- function(
    condition_name = "",
    n_sims = 1,
    include_aru = TRUE,
    include_pc = TRUE,
    include_scores = TRUE,
    include_covar = TRUE,
    nsites = 100,
    p11 = 0.1,
    p_aru11 = 0.1,
    p_aru01 = 0.025,
    beta0 = 0,
    beta1 = 1,
    mu = c(-1, 1),
    sigma = c(1, 1),
    nsurveys.aru = 24,
    nsurveys.pc = 3,
    covar_prob = 0.5,
    threshold = 0,
    covar_continuous = FALSE,
    include_aru_model = TRUE,
    include_pc_model = TRUE,
    include_scores_model = TRUE,
    include_covar_model = TRUE,
    psi_prior = "dbeta(2,2)",
    beta0_prior = "dnorm(0, 0.5)",
    beta1_prior = "dnorm(0, 0.5)",
    p11_prior = "dbeta(2, 2)",
    mu1_prior = "dnorm(-2, 0.2)",
    mu2_prior = "dnorm(-2, 0.2)",
    sigma1_prior = "dunif(0.1, 5)",
    sigma2_prior = "dunif(0.1, 5)",
    p_aru11_prior = "dbeta(2, 2)",
    p_aru01_prior = "dbeta(1, 3)I(0, 1 - p_aru11)",
    na = 1000,
    ni = 8000,
    nt = 1,
    nb = 1000,
    nc = 6) {
  jagsInputDataInit <- gen_simulated_data(
    include_aru = include_aru,
    include_pc = include_pc,
    include_scores = include_scores,
    include_covar = include_covar,
    nsites = nsites,
    p11 = p11,
    p_aru11 = p_aru11,
    p_aru01 = p_aru01,
    beta0 = beta0,
    beta1 = beta1,
    mu = mu,
    sigma = sigma,
    nsurveys.aru = nsurveys.aru,
    nsurveys.pc = nsurveys.pc,
    covar_prob = covar_prob,
    threshold = threshold,
    covar_continuous = covar_continuous
  )
  
  jags_model_fit_text <- gen_jags_model_text(
    include_aru_model = include_aru_model,
    include_pc_model = include_pc_model,
    include_scores_model = include_scores_model,
    include_covar_model = include_covar_model,
    psi_prior = psi_prior,
    beta0_prior = beta0_prior,
    beta1_prior = beta1_prior,
    p11_prior = p11_prior,
    mu1_prior = mu1_prior,
    mu2_prior = mu2_prior,
    sigma1_prior = sigma1_prior,
    sigma2_prior = sigma2_prior,
    p_aru11_prior = p_aru11_prior,
    p_aru01_prior = p_aru01_prior
  )
  
  monitored <- gen_monitored_param_list(
    include_aru_model = include_aru_model,
    include_pc_model = include_pc_model,
    include_scores_model = include_scores_model,
    include_covar_model = include_covar_model
  )
  
  jagsResultInit <- fit_model_to_sim_data(
    jagsInputData = jagsInputDataInit,
    jags_model_fit_text = jags_model_fit_text,
    monitored_list = monitored,
    nsites = nsites,
    na = na,
    ni = ni,
    nt = nt,
    nb = nb,
    nc = nc
  )
  
  all_param_list <- rownames(jagsResultInit$summary)
  sim_results_list <- vector("list", length(all_param_list))
  names(sim_results_list) <- all_param_list
  init_df <- data.frame(matrix(ncol = 11, nrow = n_sims))
  names(init_df) <- colnames(jagsResultInit$summary)
  
  for (param_name in all_param_list) {
    sim_results_list[[param_name]] <- init_df
  }
  
  convergence_fail_counter <- 0
  trials_per_sim <- rep(1, n_sims)
  max_niter <- rep(ni, n_sims)
  sim_result_dic <- rep(NA, n_sims)
  # Iterate n times
  for (i in 1:n_sims) {
    jagsInputData <- gen_simulated_data(
      include_aru = include_aru,
      include_pc = include_pc,
      include_scores = include_scores,
      include_covar = include_covar,
      nsites = nsites,
      p11 = p11,
      p_aru11 = p_aru11,
      p_aru01 = p_aru01,
      beta0 = beta0,
      beta1 = beta1,
      mu = mu,
      sigma = sigma,
      nsurveys.aru = nsurveys.aru,
      nsurveys.pc = nsurveys.pc,
      covar_prob = covar_prob,
      threshold = threshold
    )
    
    jagsResult <- fit_model_to_sim_data(
      jagsInputData = jagsInputData,
      jags_model_fit_text = jags_model_fit_text,
      monitored_list = monitored
    )
    
    f1 <- function(x) abs(x) > 1.1
    while (TRUE %in% lapply(jagsResult$Rhat, f1)) {
      convergence_fail_counter <- convergence_fail_counter + 1
      trials_per_sim[i] <- trials_per_sim[i] + 1
      max_niter[i] <- max_niter[i] + 2000
      jagsResult <- fit_model_to_sim_data(
        jagsInputData = jagsInputData,
        jags_model_fit_text = jags_model_fit_text,
        monitored_list = monitored,
        ni = max_niter[i]
      )
      if (max_niter[i] > 15000) {
        break
      }
    }
    
    for (param_name in rownames(jagsResult$summary)) {
      sim_results_list[[param_name]][i, ] <- jagsResult$summary[param_name, ]
    }
    sim_result_dic[i] <- jagsResult$DIC
  }
  
  
  sim_input_params <- expand_param_list(list(
    include_aru = include_aru,
    include_pc = include_pc,
    include_scores = include_scores,
    include_covar = include_covar,
    nsites = nsites,
    p11 = p11,
    p_aru11 = p_aru11,
    p_aru01 = p_aru01,
    beta0 = beta0,
    beta1 = beta1,
    mu = mu,
    sigma = sigma,
    nsurveys.aru = nsurveys.aru,
    nsurveys.pc = nsurveys.pc,
    covar_prob = covar_prob,
    threshold = threshold
  ))
  
  
  model_input_params <- list(
    include_aru_model = include_aru_model,
    include_pc_model = include_pc_model,
    include_scores_model = include_scores_model,
    include_covar_model = include_covar_model,
    covar_continous = covar_continuous,
    psi_prior = psi_prior,
    beta0_prior = beta0_prior,
    beta1_prior = beta1_prior,
    p11_prior = p11_prior,
    mu1_prior = mu1_prior,
    mu2_prior = mu2_prior,
    sigma1_prior = sigma1_prior,
    sigma2_prior = sigma2_prior,
    p_aru11_prior = p_aru11_prior,
    p_aru01_prior = p_aru01_prior
  )
  
  return(list(
    condition_name = condition_name,
    sim_results_list = sim_results_list,
    sim_input_params = sim_input_params,
    model_input_params = model_input_params,
    convergence_fail_counter = convergence_fail_counter,
    max_niter = max_niter,
    trials_per_sim = trials_per_sim,
    sim_result_dic = sim_result_dic
  ))
}




get_all_metric_df <- function(stacked_list, metric_name, param_list) {
  df <- data.frame(
    metric_vec = numeric(0),
    condition_name_vec = character(0),
    param_name_vec = character(0),
    true_value = numeric(0)
  )
  
  for (i in 1:length(stacked_list)) {
    for (parameter_name in param_list) {
      metric_vec <- stacked_list[[i]]$sim_results_list[[parameter_name]][[metric_name]]
      condition_name_vec <- rep(stacked_list[[i]]$condition_name, length(metric_vec))
      param_name_vec <- rep(parameter_name, length(metric_vec))
      true_value_scalar <- stacked_list[[i]]$sim_input_params[[parameter_name]]
      true_value <- rep(
        ifelse(is.null(true_value_scalar), NA, true_value_scalar),
        length(metric_vec)
      )
      df_temp <- data.frame(metric_vec, condition_name_vec, param_name_vec, true_value)
      df <- rbind(df, df_temp)
    }
  }
  return(df)
}


results_test <- run_n_simulation_and_fit(
  n_sims = 5,
  condition_name = "Base model - test"
)

results_list <- run_n_simulation_and_fit(
  n_sims = 100,
  condition_name = "Base model"
)
results_list2 <- run_n_simulation_and_fit(
  n_sims = 100,
  include_scores_model = FALSE,
  condition_name = "Base model - no scores in model"
)
results_list3 <- run_n_simulation_and_fit(
  n_sims = 100,
  include_scores_model = FALSE,
  include_aru_model = FALSE,
  condition_name = "Base model - PC Only"
)
results_list4 <- run_n_simulation_and_fit(
  n_sims = 100,
  include_pc_model = FALSE,
  condition_name = "Base model - No PC"
)
results_list5 <- run_n_simulation_and_fit(
  n_sims = 100,
  include_aru_model = FALSE,
  condition_name = "Base model - Scores + PC"
)
results_stacked <- list(
  results_list, results_list2, results_list3,
  results_list4, results_list5
)

results_df_mean <- get_all_metric_df(results_stacked, "mean", param_names)

ggplot(
  results_df_mean,
  aes(x = condition_name_vec, y = metric_vec, fill = condition_name_vec)
) +
  geom_boxplot() +
  geom_hline(aes(yintercept = true_value), size = 1) +
  facet_wrap(~param_name_vec, scale = "free") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_discrete(
    labels =
      c("PC + ARU + Scores", "No PC", "PC + ARU (no Scores)", "PC Only", "Scores + PC")
  ) +
  labs(
    title = "Posterior Mean",
    y = "Posterior Mean Parameter Value",
    fill = "Fitted Model"
  ) +
  theme(
    title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12)
  )
