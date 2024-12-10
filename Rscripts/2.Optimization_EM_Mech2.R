# This R script only optimizes the loglikelihood over the incidence and prevalence parts (glm and CoxPH)


# ---------------------------
# ----------------------------------------------------------------------
# Derive the function for the baseline survival function
# ----------------------------------------------------------------------
# ---------------------------
base_surv = function(time, delta, w, design_mat, params){
  
  event_point = time[which(w == 1 & delta == 1)] |> unique() |> sort()
  
  coxexp = (design_mat %*% params) |> exp() |> drop()
  
  lambda = numeric()
  
  event = numeric()
  
  # Estimate Λ_0(t)
  for(i in 1:length(time)){
    
    if (any(delta[which(time == time[i])] == 1 & w[which(time == time[i])] != 0, na.rm = T)){
      
      event[i] = sum(delta * as.numeric(w == 1) * as.numeric(time == time[i]), na.rm = T) # Number of events at the timepoint
      
      temp = sum(as.numeric(time >= time[i]) * w * coxexp, na.rm=T) # Denominator (Cases at risk - weighted)
      
      lambda[i] = event[i]/temp # Λ_0(t)
      
    } else {
      
      lambda[i] = 0
      
    }
    
  }
  
  HHazard = numeric()
  
  for(i in 1:length(time)){
    
    HHazard[i] = sum(as.numeric(time <= time[i]) * lambda)
    
    if (time[i] > max(event_point)){
      HHazard[i] = Inf # No survival probability after the last death timepoint (for the susceptible)
    }
    
    if (time[i] < min(event_point)){
      HHazard[i] = 0 # No risk before the first death timepoint
    }
    
  }
  
  survival = exp(-HHazard)
  
  return(survival)
  
  
}


# ---------------------------
# ----------------------------------------------------------------------
# Derive the function for the EM algorithm (known cured information)
# ----------------------------------------------------------------------
# ---------------------------

em = function(mydf, incidence_vars, latency_vars, cure_id_vars, tol = 10^-4){
  
  incidence_mat = mydf[, incidence_vars] |> as.matrix()
  latency_mat = mydf[, latency_vars] |> as.matrix()
  cure_id_mat = mydf[, cure_id_vars] |> as.matrix()
  
  difference = 100
  iter = 1
  
  while(difference > tol & iter < 100){
    
    if (iter == 1){ # Initialization of Weights and Parameters
      
      ## w_i = 1 if susceptible and uncensored; 0 otherwise
      mydf$w_i = ifelse(is.na(mydf$observed_cure_status), 0, mydf$observed_cure_status)
      
      ## Incidence part
      mod_inc = glm(w_i ~ incidence_mat, data = mydf, family = 'binomial')
      
      beta = mod_inc$coefficients
      names(beta) = gsub('incidence_mat', 'incidence_mat_', names(beta))
      
      mydf$predicted_inc = predict(mod_inc, newdata = mydf, type = 'response')
      
      ## Latency part
      mod_latency = coxph(Surv(observed_time, censor_status) ~ latency_mat + offset(log(w_i)), subset=w_i!=0, method="breslow",
                          data = mydf)
      
      beta_latency = mod_latency$coefficients
      
      names(beta_latency) = gsub('latency_mat', 'latency_mat_', names(beta_latency))
      
      mydf$base_lat = base_surv(time = mydf$observed_time,
                                delta = mydf$censor_status,
                                design_mat = latency_mat,
                                params = beta_latency,
                                w = mydf$w_i)
      
      mydf$predicted_lat = mydf$base_lat^(exp(latency_mat %*% beta_latency))
      
      ## Cure identification part
      mod_cure_id = coxph(Surv(observed_time, censor_status) ~ cure_id_mat + offset(log(1-w_i)),
                          subset=w_i!=1, method="breslow",
                          data = mydf)
      
      beta_cure_id = mod_cure_id$coefficients
      
      names(beta_cure_id) = gsub('cure_id_mat', 'cure_id_mat_', names(beta_cure_id))
      
      mydf$base_cure_id = base_surv(time = mydf$observed_time,
                                    delta = mydf$censor_status,
                                    design_mat = cure_id_mat,
                                    params = beta_cure_id,
                                    w = 1-mydf$w_i)
      
      
      mydf$predicted_cure_id = mydf$base_cure_id^(exp(cure_id_mat %*% beta_cure_id))
      
      params = c(beta, beta_latency, beta_cure_id)
      
      mydf$w_i = ifelse(!is.na(mydf$observed_cure_status),
                        mydf$observed_cure_status,
                        mydf$predicted_inc * mydf$predicted_lat / ((1-mydf$predicted_inc)*mydf$predicted_cure_id + mydf$predicted_inc * mydf$predicted_lat))
      
      
    }
    
    ## Incidence part
    mod_inc = glm(w_i ~ incidence_mat, data = mydf, family = 'binomial')
    
    beta = mod_inc$coefficients
    names(beta) = gsub('incidence_mat', 'incidence_mat_', names(beta))
    
    mydf$predicted_inc = predict(mod_inc, newdata = mydf, type = 'response')
    
    ## Latency part
    mod_latency = coxph(Surv(observed_time, censor_status) ~ latency_mat + offset(log(w_i)), subset=w_i!=0, method="breslow",
                        data = mydf)
    
    beta_latency = mod_latency$coefficients
    
    names(beta_latency) = gsub('latency_mat', 'latency_mat_', names(beta_latency))
    
    mydf$base_lat = base_surv(time = mydf$observed_time,
                              delta = mydf$censor_status,
                              design_mat = latency_mat,
                              params = beta_latency,
                              w = mydf$w_i)
    
    mydf$predicted_lat = mydf$base_lat^(exp(latency_mat %*% beta_latency))
    
    ## Cure identification part
    mod_cure_id = coxph(Surv(observed_time, censor_status) ~ cure_id_mat + offset(log(1-w_i)),
                        subset=w_i!=1, method="breslow",
                        data = mydf)
    
    beta_cure_id = mod_cure_id$coefficients
    
    names(beta_cure_id) = gsub('cure_id_mat', 'cure_id_mat_', names(beta_cure_id))
    
    mydf$base_cure_id = base_surv(time = mydf$observed_time,
                                  delta = mydf$censor_status,
                                  design_mat = cure_id_mat,
                                  params = beta_cure_id,
                                  w = 1-mydf$w_i)
    
    
    mydf$predicted_cure_id = mydf$base_cure_id^(exp(cure_id_mat %*% beta_cure_id))
    
    params_updated = c(beta, beta_latency, beta_cure_id)
    
    mydf$w_i = ifelse(!is.na(mydf$observed_cure_status),
                      mydf$observed_cure_status,
                      mydf$predicted_inc * mydf$predicted_lat / ((1-mydf$predicted_inc)*mydf$predicted_cure_id + mydf$predicted_inc * mydf$predicted_lat))
    
    
    difference = (params_updated - params)^2 |> sum() |> sqrt()
    
    params = params_updated
    
    iter = iter + 1
    # print(iter)
  }
  
  result = cbind('Coef' = names(params), 'Estimate' = params) |> as.data.frame() %>%
    
    mutate(Estimate = as.numeric(Estimate))
  
  
  return(result)
  
}



# -------------------------------------------------------------
# -------------------------------------------------------------

#                 Apply EM algorithm

# -------------------------------------------------------------
# -------------------------------------------------------------

tmp = lapply(df, function(mydf){
  
  ## Observed cure status: if uncensored -> true_cure_status, else unknown
  mydf$observed_cure_status = ifelse(mydf$censor_status == 1,
                                     mydf$true_cure_status, NA)

  result = em(mydf, incidence_vars, latency_vars, cure_id_vars)
  
  
  # Bootstrap ------------------------------------------

  
  ## Derive the subsets
  mydf_uncured = subset(mydf, observed_cure_status == 1)
  mydf_cured = subset(mydf, observed_cure_status == 0)
  mydf_censored = subset(mydf, censor_status == 0)
  
  ## Subset sample sizes
  n_uncured = nrow(mydf_uncured)
  n_cured = nrow(mydf_cured)
  n_censored = nrow(mydf_censored)
  
  i = 1
  
  bootparams = list()
  
  while (i <= nboots) {
    
    id_cured = sample(1:n_cured, n_cured, replace = TRUE)
    id_uncured = sample(1:n_uncured, n_uncured, replace = TRUE)
    id_censored = sample(1:n_censored, n_censored, replace = TRUE)
    
    
    
    bootdata = rbind(
      mydf_uncured[id_uncured, ],
      mydf_cured[id_cured, ],
      mydf_censored[id_censored, ]
      )
    
    
    bootparams[[i]] = em(bootdata, incidence_vars, latency_vars, cure_id_vars)
      
    i = i+1
  
  }
  
  
  total_results = result
  
  if (nboots > 0){
    
    boot_results = data.table::rbindlist(bootparams) %>%
      
      group_by(Coef) %>%
      
      summarise(n = n(),
                Var = var(Estimate),
                SD = sd(Estimate),
                lower_boot = quantile(Estimate, probs = 0.025),
                upper_boot = quantile(Estimate, probs = 0.975))
    
    total_results = merge(result, boot_results)
    
  }
  
  
  return(total_results)
  
  
}) 


