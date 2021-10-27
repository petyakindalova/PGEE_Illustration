### Function to run GEE
#'Inputs
#'@param y (n_subj x n_visits)-vector of binary values for the response.
#'@param X ((n_subj x n_visits) x P) design matrix containing all covariates P included in the model (incl intercept)
#'@param n_subj Number of subjects
#'@param n_visits Number of time points per subject. Assume a balanced data set.
#'@param covariance Covariance structure, options implemented are "Independence" and "Exchangeable". 
#'@param tol Tolerance for the iterative algorithm. 
#'@param max_iter Maximum number of iterations for the iterative algorithm if convergence not achieved (tolerence).
#'@param phi_est True if dispersion parameter is to be be estimates, False if phi set to 1.
#'
#'Outputs
#'@return List of results containing ..
#'
gee_run = function(y, X, n_subj, n_visits, covariance = "Independence", tol = 1e-4, max_iter = 10, phi_est = T){
  
  #if(sum(y==0)==length(y)) {return(NaN); stop("All zero response")}
  
  print_all = F
  conv_flag = F
  
  # helper fn
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  ar1_rev <- function(n, res) {
    print(sign(res))
    exponent <- 1/abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                        (1:n - 1))
    diag(exponent) = rep(1, n)
    #print(exponent)
    (abs(res)^exponent)*sign(res)
  }
  
  P = ncol(X)
  n_obs = n_subj * n_visits
  
  # Set starting values
  beta_old = c(log(mean(y)+0.0001), rep(0, P-1))
  iter_mu = exp(X%*%beta_old)
  iter_W = diag(iter_mu[,1], ncol = n_obs)
  iter_I = t(X)%*%iter_W%*%X
  iter_U = t(X)%*%(y - iter_mu)
  iter_alpha = 0
  iter_phi = 1 #added 8 Apr '21
  iter_R = matrix(iter_alpha, ncol=n_visits, nrow=n_visits)
  diag(iter_R)= rep(1, n_visits)
  iter_res = (y - iter_mu) / (sqrt(iter_mu))
  
  # Set temporary values
  # beta_old = beta_0
  # iter_mu = mu_0
  # iter_W = W_0
  # iter_I = I_0
  # iter_U = U_0
  # iter_alpha = alpha_0
  # iter_R = R_0
  # iter_res = res_0
  fisher_scoring_iterations = 0
  
  se_trace= matrix(ncol=P, nrow=max_iter)
  se_model_trace= matrix(ncol=P, nrow=max_iter)
  
  # IRLS
  for (i in 1:max_iter) {
    #print(i)
    if(i==1) {iter_I_inv = solve(iter_I)}
    
    beta_new    = beta_old + iter_I_inv%*%iter_U
    if(print_all) {print("beta"); print(beta_new)}
    iter_mu     = (exp(X%*%beta_new))
    #print(dim(iter_mu))
    
    # Covariance choice
    switch(covariance,
           # Independence={
           #   iter_alpha  = 0
           #   iter_R      = matrix(iter_alpha, ncol=n_visits, nrow=n_visits)
           #   diag(iter_R)= rep(1, n_visits)
           # },
           # AR1={
           #   iter_res = (y - iter_mu) / (sqrt(iter_mu)/sqrt(iter_phi))
           #   sum_res = 0
           #   res_temp = matrix(rep(0, n_visits*n_visits), ncol=n_visits)
           #   
           #   for (i in 1:n_subj) {
           #     idx = ((i-1)*n_visits+1):(i*n_visits)
           #     #res_temp = ar1_rev(n_visits, iter_res[idx] %*% t(iter_res[idx]))
           #     #sum_res = sum_res + sum(res_temp) - sum(diag(res_temp))
           #     res_temp = res_temp + iter_res[idx] %*% t(iter_res[idx])
           #   }
           #   
           #   for (i in 1:(n_visits - 1)){
           #     sum_res = sum_res + res_temp[i, i+1]
           #   }
           #   
           #   #iter_alpha  = (sum_res / (n_visits*(n_visits-1)))/(n_subj)
           #   iter_alpha  = iter_phi * (sum_res / ((n_visits-1)))/(n_subj)
           #   iter_R      = ar1_cor(n_visits, iter_alpha)
           #   
           # },
           Exchangeable={ 
             
             iter_res = (y - iter_mu) / (sqrt(iter_mu))
             #print(dim(iter_res^2))
             sum_res = 0
             #
             if(phi_est) {iter_phi = 1/(sum(iter_res^2)/(n_obs - P))}
             
             if(print_all) {print("phi"); print(iter_phi)}
             #
             for (i in 1:n_subj) {
               idx = ((i-1)*n_visits+1):(i*n_visits)
               res_temp = iter_res[idx] %*% t(iter_res[idx])
               sum_res = sum_res + sum(res_temp) - sum(diag(res_temp))
             }
             
             iter_alpha  = iter_phi * (sum_res / (n_visits*(n_visits-1)))/(n_subj)
             if(print_all) {print("alpha"); print(iter_alpha)}
             iter_R      = matrix(iter_alpha, ncol=n_visits, nrow=n_visits)
             diag(iter_R)= rep(1, n_visits)
           },
           {
             stop('Unknown covariance structure!')
           }
    )
    
    #iter_W_ind  = diag(iter_mu[,1], ncol=n_obs)
    
    #tmp1        = diag(sqrt(iter_mu[,1]), ncol=n_obs)
    #iter_W      = as.matrix(tmp1 %*% bdiag(rep(list(iter_R), n_subj)) %*% tmp1)
    tmp1        = diag(sqrt(1/iter_mu[,1]), ncol=n_obs)
    
    tmp2        = diag(iter_mu[,1], ncol=n_obs)%*%X
    #tmp2        = mat.mult(diag(iter_mu[,1], ncol=n_obs), X)
    tmp3        = t(tmp2)
    
    #tmp4        = spdinv(iter_W)
    if(print_all) {print("R inverse"); print(solve(iter_R))}
    iter_R_inv = solve(iter_R)
    
    if(n_subj<200) {
      tmp4      = as.matrix(tmp1 %*% bdiag(rep(list(solve(iter_R)), n_subj)) %*% tmp1) * iter_phi
    } else {
      tmp4 = matrix(rep(0, n_obs*n_obs), ncol=n_obs, nrow=n_obs)
      for (i in 1:n_subj) {
        idx_vec = ((i-1)*n_visits+1):(i*n_visits)
        #print(idx_vec)
        tmp4[idx_vec, idx_vec]  = (tmp1[idx_vec, idx_vec]%*%iter_R_inv%*%tmp1[idx_vec, idx_vec]) * iter_phi
      }
      
    }
    
    tmp5        = tmp3 %*% tmp4
    iter_I      = tmp5 %*% tmp2
    iter_U      = tmp5 %*% (y - iter_mu)
    rm(tmp1); rm(tmp3); rm(tmp4);
    gc()
    
    # Estimate variance of estimates
    # `model variance'
    if(print_all) {print("I inverse"); print(solve(iter_I))}
    iter_I_inv = solve(iter_I)
    se_model      = sqrt(diag(iter_I_inv))
    se_model_trace[fisher_scoring_iterations+1, ] = se_model
    
    # # sandwich estimator
    # temp_v1 = tmp5 %*% tmp2
    # temp_v2 = as.matrix((0^!bdiag(lapply(1:n_subj, matrix, data=1,ncol=n_visits, nrow=n_visits)))*((y-iter_mu) %*% t(y-iter_mu)))
    # temp_v3 = tmp5 %*% temp_v2 %*% t(tmp5)
    # 
    # temp_v4 = iter_I_inv %*% temp_v3 %*% iter_I_inv
    # se_sandwich = sqrt(diag(temp_v4))
    # if(print_all) {print("SE sandwich"); print(se_sandwich)}
    # se_trace[fisher_scoring_iterations+1, ] = se_sandwich
    
    if(all(abs(beta_old - beta_new)<tol)) {
      fisher_scoring_iterations = fisher_scoring_iterations + 1
      conv_flag=T
      break
    } else {
      beta_old  = beta_new
      fisher_scoring_iterations = fisher_scoring_iterations + 1
    }
  }
  
  #se_trace = se_trace[1:fisher_scoring_iterations, ]
  se_model_trace = se_model_trace[1:fisher_scoring_iterations, ]
  
  # sandwich estimator
  #temp_v1 = tmp5 %*% tmp2
  rm(tmp2)
  temp_v2 = as.matrix((0^!bdiag(lapply(1:n_subj, matrix, data=1,ncol=n_visits, nrow=n_visits)))*((y-iter_mu) %*% t(y-iter_mu)))
  temp_v3 = tmp5 %*% temp_v2 %*% t(tmp5)
  rm(tmp5)
  
  temp_v4 = iter_I_inv %*% temp_v3 %*% iter_I_inv
  se_sandwich = sqrt(diag(temp_v4))
  if(print_all) {print("SE sandwich"); print(se_sandwich)}

  
  #print(beta_new)
  
  out           = list(beta = beta_new,
                       beta_se_model = se_model,
                       beta_se_model_trace = se_model_trace,
                       beta_se_sandwich = se_sandwich,
                       #beta_se_sandwich_trace = se_trace,
                       alpha = iter_alpha,
                       phi = iter_phi,
                       iterations = fisher_scoring_iterations,
                       conv = conv_flag)
  
  out
  
}
