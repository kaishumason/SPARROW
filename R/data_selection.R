
#X is the data matrix(p by n where p is the number of covariates and n is the number of data points)
#data size is how many points we wish to select

#expand covariate matrix based on X, lambda, family,and min expresson
expand_covariate_matrix = function(X,lambda,family,keep_coef = matrix(1,ncol(X),ncol(lambda)),
                                   lib_size = rep(1,nrow(X)), min_reads_per_1000 = 1,min_freq = 500){
  #get scaled library sizes
  if(family == "poisson" | family == "negative binomial"){
    scaled_spot_size = sqrt((lib_size/1000)*min_reads_per_1000)
    X = sweep(X,1,scaled_spot_size,"*")
  } else if(family == "binomial"){
    scaled_spot_size = sqrt((lib_size/1000)*min_reads_per_1000)
    X = sweep(X,1,scaled_spot_size,"*")
  }

  #initialize expanded covariate matrix
  big_X  = list()
  #iterate over each cell type
  for(j in c(1:ncol(lambda))){
    #get which index to remove
    bad_ind = which(keep_coef[,j] == 0)
    #get scaled X
    scaled_X = sweep(X,1,lambda[,j],"*")
    colnames(scaled_X) = paste0(colnames(lambda)[j],":",colnames(X))
    #remove index
    if(length(bad_ind)>0){
      scaled_X = scaled_X[,-bad_ind]
    }
    #save value
    big_X[[j]] = scaled_X
  }
  big_X = do.call(cbind, big_X)
  #remove those columns that don't appear enough
  freq = apply(big_X,2,function(x){sum(x != 0)})
  big_X = big_X[,(freq > min_freq)]

  return(big_X)
}

#compute target standard erros
compute_target_standard_error = function(X,min_effect,target_power_approx,alpha,acc = 1e-3){
  if(length(min_effect) == 1){
    min_effect = rep(min_effect,ncol(X))
  }
  #step 1: get global min standard error (make signal/noise at least 3 for smallest effects)
  global_min_SE = SE_optimizer(min_effect,pow = 0.999,alpha = 0.05, acc = 1e-3)

  #step 2: Get standard errors on full dataset
  A_full = solve(t(X)%*%X)
  SE_full = sqrt(diag(A_full))

  #get number of covariates
  nSE = length(SE_full)
  #initialize list of minimum SE
  ct_min_SE = rep(NA,nSE)
  #iterate over covariates
  for(j in c(1:nSE)){
    #compute power on min effect
    ct_min_power = power_calc(SE_full[j],min_effect[j],0.05)
    #if we have high for smallest effect
    if(ct_min_power > (0.999)){
      ct_min_SE[j] = global_min_SE[j]
    }else{
      #get max effect size for our standard error
      max_effect = min_effect_optimizer(SE_full[j],0.999,alpha = 0.05)
      #compute power of SE full
      ct_min_SE[j] = SE_power_optimizer(SE_full[j],target_power_approx,min_effect[j],max_effect,alpha = 0.05)
    }
  }

  #get cell type specific global SE
  ct_min_SE[ct_min_SE < global_min_SE] = global_min_SE[ct_min_SE < global_min_SE]

  return(ct_min_SE)
}

#function to get subset of data of size data size
#X is of dimension #cov times #observations
data_selection = function(X,data_size,log = FALSE,min_SE = NULL){
  #get number of points and number of covariates
  n = ncol(X)
  p = nrow(X)
  #get minimum batch size such that we sacrifice less than 1/1000 of efficiency
  subsample_size = ceiling(log(1000)*n/data_size)
  #get initial varcov matrix
  initial_sigma = diag(x = 100, p)
  #initialize list of points we choose
  points_added = rep(NA,data_size)
  #keep track of how mnay points we have left to pick from
  points_left = c(1:n)
  t1 = Sys.time()
  for(j in c(1:data_size)){
    if(log & j%%1000 == 0){
      print(paste0("On iteration ",j," out of ",data_size))
      print(Sys.time() - t1)
    }
    #subsample points
    m = length(points_left)
    #sample points
    samp = sample.int(m, min(subsample_size, m), replace = FALSE)
    #get sampled indices
    sampled_inds = points_left[samp]
    #get the corresponding x points
    big_X_temp = X[,sampled_inds]
    #take transpose
    big_X_temp = t(big_X_temp)


    #check gain for each data point
    A = big_X_temp%*%initial_sigma
    #numerator and denominator of gain
    num = rowSums(A^2)
    denom = 1 + rowSums(A*big_X_temp)
    #identify best index
    ind = which.max(num/denom)
    if(length(ind) > 1){
      ind = sample(ind,1)
    }
    #get best index
    best_cand = sampled_inds[ind]

    #add point to points added
    points_added[j] = best_cand
    #compute new covaraince matrix via woodbury's lemma
    x = t(big_X_temp[ind,,drop = FALSE])
    A = initial_sigma%*%x
    initial_sigma = initial_sigma - A%*%t(A)/(1 + t(x)%*%A)[1,1]
    #remove point from points left '
    points_left = points_left[points_left != best_cand]
  }
  return(points_added)
}

#function to get subset of data of size data size
#with stopping criterion if all SE is less than min_SE
#X is of dimension #cov times #observations
data_selection_terminal = function(X,data_size,log = FALSE,min_SE = NULL){
  #scale by min_SE
  X = sweep(X,1,min_SE,"*")
  #get number of points and number of covariates
  n = ncol(X)
  p = nrow(X)
  #get frequencies of how many covaraites are non null
  freqs = colSums(X > 0)
  #keep track of how many covariates converged
  conv_cov = c()
  #get minimum batch size such that we sacrifice less than 1/1000 of efficiency
  subsample_size = ceiling(log(1000)*n/data_size)
  #get initial varcov matrix
  initial_sigma = diag(x = 100, p)
  #initialize list of points we choose
  points_added = rep(NA,data_size)
  #keep track of how mnay points we have left to pick from
  valid_points_left = c(1:n)
  invalid_points_left = c()

  print("starting")
  t1 = Sys.time()
  for(j in c(1:data_size)){
    if(length(conv_cov) == p){
      print("ended early")
      break
    }
    if(log & j%%1000 == 0){
      print(paste0("On iteration ",j," out of ",data_size))
      if(is.null(min_SE) == FALSE){
        print("Ratio of Current SE to Target SE")
        stand_errs = sqrt(diag(initial_sigma))
        print(stand_errs)
      }
      print(Sys.time() - t1)
    }

    #subsample points
    m_eff <- length(valid_points_left)
    sampled_inds <- sample(valid_points_left, min(c(subsample_size, m_eff)), replace = FALSE)

    #get the corresponding x points
    big_X_temp = X[,sampled_inds,drop = F]
    #take transpose
    big_X_temp = t(big_X_temp)

    #check gain for each data point
    A = big_X_temp%*%initial_sigma
    #numerator and denominator of gain
    num = rowSums(A^2)
    denom = 1 + rowSums(A*big_X_temp)

    #identify best index
    ind = which.max(num/denom)

    if(length(ind) > 1){
      ind = sample(ind,1)
    }
    #get best index
    best_cand = sampled_inds[ind]

    #add point to points added
    points_added[j] = best_cand
    #compute new covaraince matrix via woodbury's lemma
    x = t(big_X_temp[ind,,drop = FALSE])
    A = initial_sigma%*%x
    initial_sigma = initial_sigma - A%*%t(A)/(1 + t(x)%*%A)[1,1]
    if(is.null(min_SE) == FALSE){
      stand_errs = sqrt(diag(initial_sigma))
      conv_cov_curr = which(stand_errs <=1)
      new_covs = setdiff(conv_cov_curr, conv_cov)
      #update freq accordingly
      if(length(new_covs) > 0){
        for(cov in new_covs){
          freqs = freqs - as.integer(X[cov,] != 0)
        }
      }
      #get covaraites that are removes
      removed_covs = setdiff(conv_cov,conv_cov_curr)
      if(length(removed_covs) > 0){
        for(cov in removed_covs){
          freqs = freqs + as.integer(X[cov,] != 0)
        }
      }
    }
    #remove point from points left '
    valid_points_left <- valid_points_left[valid_points_left != best_cand]
    #update converged convs
    conv_cov <- conv_cov_curr
    #update weights if needed
    if(length(removed_covs) > 0 | length(new_covs) > 0){
      points_left = c(valid_points_left,invalid_points_left)
      valid_points_left = points_left[freqs[points_left] != 0]
      invalid_points_left = points_left[freqs[points_left] == 0]
    }
  }
  return(points_added)
}



#for a given effect size, desired power, and type 1 error rate, compute min SE
SE_optimizer = function(min_effect,pow,alpha,acc = 1e-3){
  #assume effect
  lower = 0
  upper = 1000
  z = qnorm(1-alpha/2)
  while((upper - lower) > acc){
    mid = (lower + upper)/ 2
    val = pnorm(z - mid) - pnorm(-z - mid)
    if(val > (1-pow)){
      lower = mid
    }else{
      upper = mid
    }
  }
  return(min_effect/upper)
}


#for a given effect size and type 1 error rate, compute min effect that can be found at desired power
min_effect_optimizer = function(SE,pow,alpha,acc = 1e-3){
  #assume effect
  lower = 0
  upper = 1000
  z = qnorm(1-alpha/2)
  while((upper - lower) > acc){
    mid = (lower + upper)/ 2
    val = pnorm(z - mid) - pnorm(-z - mid)
    if(val > (1-pow)){
      lower = mid
    }else{
      upper = mid
    }
  }
  return(SE*upper)
}


#for a given SE and effect size, compute power
power_calc = function(SE,effect,alpha){
  z = qnorm(1-alpha/2)
  1 - pnorm(z - effect/SE) - pnorm(-z - effect/SE)

  return(1 - pnorm(z - effect/SE) - pnorm(-z - effect/SE))
}

#for a given SE and effect size, compute power  integrsl
power_integral = function(SE,effect_min,effect_max,alpha){
  #get what we're integrating over
  nbox = 1000
  x = seq(effect_min,effect_max,length.out = nbox)
  #evaluate integral
  powers = power_calc(SE,x,alpha)
  #return total power
  return(mean(powers))
}


#For a given target SE, find min SE that recovers target_power_approx of power
SE_power_optimizer = function(target_SE,target_power_approx,effect_min,effect_max,alpha,acc = 1e-3){
  upper = 100
  lower = 1
  #get target power
  target_power = power_integral(target_SE,effect_min,effect_max,alpha)
  #find min standard error that approximates power within some
  while((upper - lower) > acc){
    mid = (lower + upper)/ 2
    #get power on suggested new SE
    val = power_integral(target_SE*mid,effect_min,effect_max,alpha = 0.05)
    #evaluate ratio
    power_ratio = val/target_power
    #if we're too powerful, can have larger SE
    if(power_ratio > target_power_approx){
      lower = mid
    }else{
      upper = mid
    }
  }
  return(target_SE*upper)

}


#function to simulate data
simulate_data = function(n,nct,effect_scale,intercept_scale,library_size ,spot_ct,spot_size,p,num_null = 3,prob_ct = NULL){
  if(is.null(prob_ct)){
    prob_ct = rep(1,nct)
  }
  #covariate matrix
  X = matrix(rnorm(p*n)*rbinom(p*n,1,prob = 0.25),n,p)
  X = cbind(1,X)
  #true beta
  beta = matrix((rnorm(nct*ncol(X))),ncol(X),nct)*effect_scale
  null_beta = matrix(0, dim(beta)[1],dim(beta)[2])
  #set first row (intercept) to be large
  beta[1,] = -intercept_scale*abs(runif(nct,1,1.5))
  #make random covariates null
  null_cand = nct*ncol(X)
  #get intercept indices
  intercepts = 1 + nrow(beta)*c(0:(ncol(beta)-1))
  #get candidates for effects that can be null (i.e non intercepts)
  null_cand = c(1:null_cand)[-intercepts]
  #sample null effects
  random_null = sample(null_cand,num_null * nct, replace = F)
  #set these coefs to be 0
  null_beta[random_null] = 1
  beta[random_null] = 0
  #get lambda
  lambda = matrix(0,n,nct)
  #get cell type assignment
  CT = rep(NA,n)
  for(k in c(1:n)){
    #sample cell types
    L = sample(1:nct,size = spot_ct, replace = TRUE,prob = prob_ct)
    counter = 1
    #attribute cells to each cell type
    for(val in L){
      lambda[k,val] = lambda[k,val] + counter
      counter = counter + 1
    }
    #scale number of cells by spot size
    lambda[k,] = rdirichlet(1, 2.5*lambda[k,]/sum(lambda[k,]))
    CT[k] = which.max(lambda[k,])
  }

  #Calculate eta for all other cell types
  eta_rest = X%*%beta
  #Calculate mu for all other cell types
  mu_rest = spot_poisson[['marginal_mu']](eta_rest)
  #Calculate C for all spots
  C = rowSums(lambda*mu_rest)
  #simulate response
  y = rpois(n = length(C),lambda = C*library_size)

  #return data
  return(list(y = y, X = X,lambda = lambda ,beta = beta,null_beta = null_beta, CT = CT))
}











