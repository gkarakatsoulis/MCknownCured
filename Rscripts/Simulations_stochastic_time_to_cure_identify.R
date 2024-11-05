library(dplyr)


# To simulate the (randomly right) censored observations, we need to first simulate a lifetime
# vector and, independently, a termination time (=censoring time) vector. Then we observe
# whichever comes first (ztimes) and related indicator (status)
#
# For the cured individuals, we simulate a time to cure identification 
# Then, we follow the same process as above


myfunction = function(N, exp_rate, seed, prob_known_cured, incidence_par, latency_par, cure_id_par){
  
  # We start with determining the population ------------------------------
  set.seed(seed)
  
  # Total sample size
  # N = 500
  
  # -------------------------------------------------------------------
  #                       1. Model Covariates generation
  # -------------------------------------------------------------------
  
  
  # a. Generate the incidence part risk factors (one continuous, one binary)
  z1 = rnorm(N, 0, 1)
  z2 = rbinom(n = N, size = 1, prob = 0.5)
  
  # b. Generate the both parts risk factors (one continuous, one binary)
  q1 = rnorm(N, 0, 1)
  q2 = rbinom(N, 1, 0.7)
  
  # c. Generate the latency part risk factors (one continuous, one binary)
  x1 = rnorm(N, 0, 1)
  x2 = rbinom(N, 1, 0.6)
  
  mydata = as.data.frame(cbind('x1' = x1, 'x2' = x2, 'q1' = q1, 'q2' = q2))
  
  
  # -------------------------------------------------------------------
  #                       2. Logistic regression part
  # -------------------------------------------------------------------
  
  # a. Generate the logistic regression
  logit_incidence = 2 + cbind(z1, z2, q1, q2) %*% incidence_par
  
  # b. Calculate the probabilities
  prob_incidence = plogis(logit_incidence)
  
  # c. Generate the vector of the cure status (Y_i: 1 if uncured, 0 if cured)
  mydata$y_i = rbinom(n = N, size = 1, prob = prob_incidence)
  
  
  # -------------------------------------------------------------------
  #                       3. Censoring times Mechanism
  # -------------------------------------------------------------------
  
  
  # a. Generate the censoring times (coming from exponential distribution)
  mydata$cens_time = rexp(n = N, rate = exp_rate)
  
  
  # -------------------------------------------------------------------
  #                       4. Time-to-event for susceptible
  # -------------------------------------------------------------------
  
  shape_par = 2
  scale_par = 2 
  X = mydata[, c('x1', 'x2', 'q1', 'q2')] |> as.matrix()
  scale_par_indiv = scale_par * exp(- X %*% (latency_par / shape_par))
  
  mydata$susc_time = rweibull(N, shape = shape_par, scale = scale_par_indiv)
  
  
  # -------------------------------------------------------------------
  #                       5. Time-to-event for uncured identification
  # -------------------------------------------------------------------
  
  shape_par_cure_id = 2
  scale_par_cure_id = 1 
  X = mydata[, c('x1')] |> as.matrix()
  
  scale_par_indiv_cure_id = scale_par_cure_id * exp(- X %*% (cure_id_par / shape_par_cure_id))
  
  mydata$cure_id_time = rweibull(N, shape = shape_par_cure_id, scale = scale_par_indiv_cure_id)
  
  # -------------------------------------------------------------------
  #                       6. Final Data Generation
  # -------------------------------------------------------------------
  
  mydata = mydata %>%
    
    transmute(
      
      true_cure_status = y_i,
      true_time = ifelse(true_cure_status == 1,
                         susc_time,
                         cure_id_time),
      cens_time,
      censor_status = ifelse(true_time <= cens_time, 1, 0),
      observed_time = ifelse(censor_status == 1,
                             true_time,
                             cens_time),
      z1, z2, q1, q2, x1, x2,     
      censor_status0 = ifelse(censor_status == 1 & true_cure_status == 1, 1, 0),
      observed_cure_status0 = ifelse(censor_status0 == 1, 1, NA),
      observed_cure_status100 = true_cure_status,
      censor_status100 = censor_status
      
      )
  
  
  
  results = list()
  
  results[['Data']] = df
  
  results[['Parameters']] = cbind(
    
    c('Sample_size', 'Censor', 'Beta', 'Gamma'),
    c(N, paste0('Exp(', exp_rate, ')'),
      paste0(incidence_par, collapse = '_'),
      paste0(latency_par, collapse = '_'))
    
    
  ) %>%
    
    as.data.frame()
  
  colnames(results$Parameters) = c('Parameter', 'Value')
  
  prob_known_cured = sum(mydata$censor_status == 1 & mydata$true_cure_status == 1) / sum(mydata$true_cure_status == 1)*100
  
  file_name = paste0('Mech2_', 'N_', N,
                     '_Censor_Exp(', exp_rate, ')',
                     '_Prob_known_', prob_known_cured |> round(),
                     '_Beta_(', paste0(incidence_par, collapse = '_'),
                     ')_Gamma_(', paste0(latency_par, collapse = '_'), ')')
  
  myfile_save = gsub('Rscripts/Simulations', paste0('Simulated_Datasets/', file_name,'//'), getwd())
  
  ifelse(!dir.exists(myfile_save), dir.create(myfile_save), FALSE)
  
  
  openxlsx::write.xlsx(results, paste0(myfile_save, 'Dataset_', seed, '.xlsx'))
  
  return(results)
  
  
}

myseed = 1:100

for (j in myseed){
    
    myfunction(N = 500,
               exp_rate = 1/5,
               seed = j,
               incidence_par = c(1.5, -2, 4, 0.5),
               latency_par = c(0.9, 4, -2, 2),
               cure_id_par = 0.5)
    
}
  





