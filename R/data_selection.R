#' Expand Covariate Matrix Based on Cell Type Proportions
#'
#' Expands the design matrix \code{X} using cell type proportions \code{lambda},
#' applying scaling for sequencing depth and removing low-frequency features.
#'
#' @param X Covariate matrix of dimension (n times p)
#' @param lambda Matrix of cell type proportions (n x k)
#' @param family Model family: "poisson", "negative binomial", or "gaussian"
#' @param keep_coef Matrix (p by k) indicating which coefficients to keep in case one knows what to remove a priori(default: all TRUEs)
#' @param lib_size Library sizes for each sample (default: 1 for all)
#' @param min_reads_per_1000 If family == "poisson" or "negative binomial" Minimum reads per 1000 for robust scaling (default: 1)
#' @param max_noise_sd If family == "gaussian" maximum standard deviation of the error term to be used for robust scaling (default: 1)
#'
#' @return A scaled and expanded covariate matrix
#' @export
expand_covariate_matrix = function(X,lambda,family = "gaussian",keep_coef = matrix(TRUE,ncol(X),ncol(lambda)),lib_size = rep(1,nrow(X)), min_reads_per_1000 = 1,max_noise_sd = 1){
  #get scaled library sizes
  if(family == "poisson" | family == "negative binomial"){
    scaled_spot_size = sqrt((lib_size/1000)*min_reads_per_1000)
    X = sweep(X,1,scaled_spot_size,"*")
  } else if(family == "gaussian"){
    X = X/max_noise_sd
  }else{
    stop("family should be one of poisson, negative binomial, or gaussian\n
        If you don't wish to apply any transformation to X, use gaussian(default)")
  }

  #initialize expanded covariate matrix
  big_X  = list()
  #iterate over each cell type
  for(j in c(1:ncol(lambda))){
    #get which index to remove
    bad_ind = which(keep_coef[,j] == FALSE)
    #get scaled X
    scaled_X = sweep(X,1,lambda[,j],"*")
    if(is.null(colnames(lambda)) == F & is.null(colnames(X)) == F){
      colnames(scaled_X) = paste0(colnames(lambda)[j],":",colnames(X))
    }
    #remove index
    if(length(bad_ind)>0){
      scaled_X = scaled_X[,-bad_ind]
    }
    #save value
    big_X[[j]] = scaled_X
  }
  big_X = do.call(cbind, big_X)

  return(big_X)
}



#' Compute Target Standard Errors for Covariates
#'
#' Estimates the minimum required standard errors for detecting given effect sizes with desired power.
#'
#' @param X Design matrix
#' @param min_effect Minimum effect size(s) to detect
#' @param target_power_approx Approximate desired statistical power (0-1)
#' @param alpha Type I error rate (e.g., 0.05)
#' @param acc Accuracy for numeric search (default: 1e-3)
#'
#' @return A vector of minimum standard errors for each covariate
#' @export
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







#' Select Data with Early Stopping Based on Standard Errors
#'
#' Performs greedy data selection with an early stopping criterion when standard errors converge.
#'
#' @param X Covariate matrix (p x n)
#' @param max_data_size maximum number of data points to select
#' @param min_standard_error Target standard errors for convergence
#' @param log Whether to print progress
#' @param period If log = TRUE, print progress every period
#'
#' @return A vector of selected data point indices
#' @export
data_selection = function(X,max_data_size,min_standard_error = NULL,log = FALSE,period = 1000){
  #scale by min_standard_error
  X = sweep(X,1,min_standard_error,"*")
  #get number of points and number of covariates
  n = ncol(X)
  p = nrow(X)
  #get frequencies of how many covaraites are non null
  freqs = colSums(X > 0)
  #keep track of how many covariates converged
  conv_cov = c()
  #get minimum batch size such that we sacrifice less than 1/1000 of efficiency
  subsample_size = ceiling(log(1000)*n/max_data_size)
  #get initial varcov matrix
  initial_sigma = diag(x = 100, p)
  #initialize list of points we choose
  points_added = rep(NA,max_data_size)
  #keep track of how mnay points we have left to pick from
  valid_points_left = c(1:n)
  invalid_points_left = c()

  print("starting")
  t1 = Sys.time()
  for(j in c(1:max_data_size)){
    if(length(conv_cov) == p){
      print("ended early")
      break
    }
    if(log & j%%period == 0){
      print(paste0("On iteration ",j," out of ",max_data_size))
      if(is.null(min_standard_error) == FALSE){
        print("Ratio of Current standard error to Target standard error")
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
    if(is.null(min_standard_error) == FALSE){
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
    if(is.null(min_standard_error) == FALSE){
      if(length(removed_covs) > 0 | length(new_covs) > 0){
        points_left = c(valid_points_left,invalid_points_left)
        valid_points_left = points_left[freqs[points_left] != 0]
        invalid_points_left = points_left[freqs[points_left] == 0]
      }
    }
  }
  return(points_added)
}

#' Select Informative Data Points for Modeling
#'
#' Selects a subset of data points to minimize standard error using greedy optimization.
#'
#' @param X Covariate matrix (p x n)
#' @param data_size Number of data points to select
#' @param log Whether to print progress every 1000 iterations
#' @param period If log = TRUE, print progress every period
#'
#' @return A vector of selected data point indices
#' @export
data_selection_fixed = function(X,data_size,log = FALSE,period = 1000){
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
    if(log & j%%period == 0){
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





#' Compute Minimum Standard Error for Given Power and Effect
#'
#' Calculates the minimum standard error needed to detect a given effect with a specified power.
#'
#' @param min_effect Minimum effect size
#' @param pow Desired power (e.g., 0.8 or 0.999)
#' @param alpha Type I error rate
#' @param acc Accuracy of the search
#'
#' @importFrom stats pnorm qnorm
#' @keywords internal
SE_optimizer = function(min_effect,pow,alpha = 0.05,acc = 1e-3){
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


#' Compute Minimum Detectable Effect Size
#'
#' Calculates the minimum detectable effect size given standard error and desired power.
#'
#' @param SE Standard error
#' @param pow Desired power
#' @param alpha Type I error rate
#' @param acc Accuracy of search
#'
#' @importFrom stats pnorm qnorm
#' @keywords internal
min_effect_optimizer = function(SE,pow,alpha = 0.05,acc = 1e-3){
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


#' Compute Statistical Power
#'
#' Calculates the power to detect an effect given a standard error and significance level.
#'
#' @param SE Standard error
#' @param effect Effect size
#' @param alpha Type I error rate
#'
#' @importFrom stats pnorm qnorm
#' @keywords internal
power_calc = function(SE,effect,alpha = 0.05){
  z = qnorm(1-alpha/2)
  1 - pnorm(z - effect/SE) - pnorm(-z - effect/SE)

  return(1 - pnorm(z - effect/SE) - pnorm(-z - effect/SE))
}

#' Average Power Over a Range of Effect Sizes
#'
#' Computes the average power over a uniform distribution of effect sizes.
#'
#' @param SE Standard error
#' @param effect_min Minimum effect size
#' @param effect_max Maximum effect size
#' @param alpha Type I error rate
#'
#' @keywords internal
power_integral = function(SE,effect_min,effect_max,alpha){
  #get what we're integrating over
  nbox = 1000
  x = seq(effect_min,effect_max,length.out = nbox)
  #evaluate integral
  powers = power_calc(SE,x,alpha)
  #return total power
  return(mean(powers))
}


#' Optimize Standard Error for Target Power Approximation
#'
#' Finds the smallest standard error that achieves a desired approximation to target power.
#'
#' @param target_SE Initial estimate of standard error
#' @param target_power_approx Desired power ratio (e.g., 0.95)
#' @param effect_min Minimum effect size
#' @param effect_max Maximum effect size
#' @param alpha Type I error rate
#' @param acc Accuracy for numeric optimization
#'
#' @keywords internal
SE_power_optimizer = function(target_SE,target_power_approx,effect_min,effect_max,alpha = 0.05,acc = 1e-3){
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


#' Simulate gene expression data
#'
#' This function simulates gene expression data using a spot-based model
#'
#' @param n Number of data points
#' @param nct Number of cell types
#' @param effect_range Range of effect sizes (effect sizes drawn uniformly on this interval)
#' @param intercept_range Range for intercepts (intercept sizes drawn uniformly on this interval)
#' @param min_effect Minimum absolute value of effect sizes.
#' @param library_size Number of transcripts per cell/spot
#' @param spot_ct Number of cell types per spot
#' @param p Number of covariates
#' @param num_null Number of null coefficients to include
#' @param prob_ct Optional probabilities for cell type sampling
#' @param family The data generating distribution. One of poisson, negative binomial,binomial, and gaussian
#' @param dispersion Dispersion parameter for NB or Gaussian family. Size parameter for NB and sd for Gaussian. Default: 1
#'
#' @return A list with simulated y, X, lambda, beta, null_beta, and CT
#' @importFrom gtools rdirichlet
#' @importFrom stats pnorm qnorm rbinom rnorm rpois runif
#' @export
simulate_data = function(n,nct,effect_range,intercept_range,min_effect = 0.05,library_size = 500,
                         spot_ct = min(2,nct),p = 6,num_null = 2,prob_ct = NULL,family = "poisson",
                         dispersion = 1){
  if(is.null(prob_ct)){
    prob_ct = rep(1,nct)
  }
  #reorder effect range and intercept range
  effect_range = sort(effect_range)
  intercept_range = sort(intercept_range)
  #function to sample from ranges
  sample_from_interval_diff <- function(a, b, c, n = 1) {
    interval_difference <- function(a, b, c) {
      if (a > b) stop("Error: a must be less than b.")
      intervals <- list()
      if (b <= -c || a >= c) {
        intervals <- list(c(a, b))
      } else {
        if (a < -c) intervals <- append(intervals, list(c(a, -c)))
        if (b > c) intervals <- append(intervals, list(c(c, b)))
      }
      intervals
    }

    intervals <- interval_difference(a, b, c)

    if (length(intervals) == 0) {
      stop("No intervals available for sampling. Make sure maximum of effect range is not less than min effect.")
    }

    # Compute lengths of intervals
    lengths <- sapply(intervals, function(x) x[2] - x[1])
    probs <- lengths / sum(lengths)

    # Sample intervals proportional to their lengths
    chosen_intervals <- sample(seq_along(intervals), size = n, replace = TRUE, prob = probs)

    # Sample uniformly within chosen intervals
    samples <- sapply(chosen_intervals, function(idx) {
      runif(1, intervals[[idx]][1], intervals[[idx]][2])
    })

    return(samples)
  }

  #covariate matrix
  X = matrix(rnorm(p*n)*rbinom(p*n,1,prob = 0.25),n,p)
  X = cbind(1,X)
  #true beta
  beta_samples = sample_from_interval_diff(a = effect_range[1],b = effect_range[2],c = min_effect, n = nct*ncol(X))
  beta = matrix(beta_samples,ncol(X),nct)
  null_beta = matrix(0, dim(beta)[1],dim(beta)[2])
  #set first row (intercept) to be large
  beta[1,] = (runif(nct,intercept_range[1],intercept_range[2]))
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
    lambda[k,] = gtools::rdirichlet(1, 2.5*lambda[k,]/sum(lambda[k,]))
    CT[k] = which.max(lambda[k,])
  }

  #Calculate eta for all other cell types
  eta_rest = X%*%beta
  #Calculate mu for all other cell types
  if(family == "poisson"){
    mu_rest = exp(eta_rest)
    #Calculate C for all spots
    C = rowSums(lambda*mu_rest)
    #simulate response
    y = rpois(n = length(C),lambda = C*library_size)
  }else if(family == "negative binomial"){
    mu_rest = exp(eta_rest)
    #Calculate C for all spots
    C = rowSums(lambda*mu_rest)
    #simulate response
    y = stats::rnbinom(n = length(C),mu = C*library_size, size = dispersion)
  }else if(family == "binomial"){
    mu_rest = gtools::inv.logit(eta_rest)
    #Calculate C for all spots
    C = rowSums(lambda*mu_rest)
    #simulate response
    y = rbinom(n = length(C),size = library_size, prob = C)
  }else if(family == "gaussian"){
    mu_rest = eta_rest
    #Calculate C for all spots
    C = rowSums(lambda*mu_rest)
    #simulate response
    y = rnorm(n = length(C),mean = C*library_size,sd = dispersion)
  }else{
    stop("Family must be one of poisson, binomiaa, negative binomial, or gaussian")
  }


  #return data
  return(list(y = y, X = X,lambda = lambda ,beta = beta,null_beta = null_beta, CT = CT))
}



#' Read MERFISH Example Data from GitHub
#'
#' Downloads and loads MERFISH example data from a public GitHub repository.
#' The function retrieves cell type annotations, spatial region labels, and
#' gene expression count matrices (in up to 10 chunks), then returns them as a list.
#'
#' @param num_chunks Integer between 0 and 10. Controls how many count matrix chunks
#' to download and load. Defaults to 10 (all chunks). Use 0 to skip loading counts.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{counts}{A matrix or data frame of gene expression counts (combined from multiple chunks).}
#'   \item{regions}{A vector or factor of region labels for each cell.}
#'   \item{CT}{A vector or factor of cell types for each cell.}
#' }
#'
#' @details This function requires an internet connection to download data from
#' \url{https://github.com/kaishumason/SpotGLM-Example-Data}. The function assumes the
#' repository structure and filenames are consistent with the expected format.
#'
#' @examples
#' \dontrun{
#' data_list <- read_merfish_data(num_chunks = 5)
#' head(data_list$counts)
#' table(data_list$CT)
#' }
#'
#' @export
read_merfish = function(num_chunks = 10){
  #get cell types
  # Define raw GitHub URL
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/merfish/cell_types.rds"

  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  cell_types <- readRDS(temp_file)


  #get regions
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/merfish/regions.rds"

  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  regions <- readRDS(temp_file)

  #get counts
  n = min(as.integer(num_chunks),10)
  if(n<0){
    stop("number of chunks should be an integer between 0(no counts) and 10 (all counts)")
  }
  data = vector("list",n)
  for(j in c(1:n)){
    url <- paste0("https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/merfish/counts_",j,".rds")
    # Temporary file to store the .rds
    temp_file <- tempfile(fileext = ".rds")

    # Download the file (use mode = "wb" for binary)
    download.file(url, destfile = temp_file, mode = "wb")

    data[[j]] = readRDS(temp_file)
  }
  data = do.call(rbind, data)

  return(list(counts = data,regions = regions, CT = cell_types))
}



#' Read Visium HD Example Data from GitHub
#'
#' Downloads and loads spatial transcriptomics data from the Visium HD example
#' dataset hosted on a public GitHub repository. This includes spatial coordinates,
#' gene expression counts, deconvolution results, and effective niche estimates.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{coords}{A matrix or data frame of spatial coordinates (x and y) for each spot.}
#'   \item{niche}{A matrix, data frame, or list representing effective niche composition per spot.}
#'   \item{deconv}{A matrix or data frame of cell type deconvolution proportions per spot.}
#'   \item{counts}{A gene expression count matrix (genes × spots).}
#' }
#'
#' @details This function requires an internet connection to download data from
#' \url{https://github.com/kaishumason/SpotGLM-Example-Data}. The repository must contain
#' the files \code{coords.rds}, \code{deconv_matrix.rds}, \code{count_matrix_subset.rds},
#' and \code{niche.rds} in the \code{visiumHD} folder.
#'
#' @examples
#' \dontrun{
#' data_list <- read_visiumHD()
#' head(data_list$coords)
#' dim(data_list$counts)
#' }
#'
#' @export
read_visiumHD = function(){
  #read coordinates
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/coords.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  coords <- as.matrix(readRDS(temp_file)[,4:5])

  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/deconv.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  deconv <- as.matrix(readRDS(temp_file))


  #read counts
  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/counts.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  counts <- readRDS(temp_file)


  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/niche.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  niche <- readRDS(temp_file)

  #read library size
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/visiumHD/library_size.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  library_size <- readRDS(temp_file)


  return(list(coords = coords,niche = niche, deconv = deconv, counts = counts,library_size = library_size))

}






#' Read Visium Example Data from GitHub
#'
#' Downloads and loads Visium spatial transcriptomics data from a GitHub repository.
#' Includes coordinates, deconvolution, effective niche covariates, library sizes, and gene counts.
#'
#' @return A list containing:
#' \describe{
#'   \item{coords}{Matrix of spatial coordinates.}
#'   \item{niche}{Effective niche covariate matrix.}
#'   \item{deconv}{Cell type deconvolution matrix.}
#'   \item{counts}{Gene expression count matrix.}
#'   \item{library_size}{Vector of library sizes per spot.}
#' }
#'
#' @examples
#' \dontrun{
#' visium_data <- read_visium()
#' str(visium_data)
#' }
#'
#' @export
read_visium = function(){
  #read coordinates
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/coords.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  coords <- as.matrix(readRDS(temp_file))


  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/deconv.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  deconv <- as.matrix(readRDS(temp_file))



  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/effective_niche_covariates.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  niche <- readRDS(temp_file)



  #read library size
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/library_size.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  library_size <- readRDS(temp_file)


  #get counts
  data = vector("list",4)
  for(j in c(1:4)){
    url <- paste0("https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Visium/counts_",j,".rds")
    # Temporary file to store the .rds
    temp_file <- tempfile(fileext = ".rds")

    # Download the file (use mode = "wb" for binary)
    download.file(url, destfile = temp_file, mode = "wb")

    data[[j]] = readRDS(temp_file)
  }
  data = do.call(cbind, data)



  return(list(coords = coords,niche = niche, deconv = deconv, counts = data,library_size = library_size))

}



#' Read Spatial ATAC-Seq Example Data
#'
#' Loads spatial ATAC-seq data from a GitHub repository, including coordinates, deconvolution,
#' region-level features, and motif score matrices.
#'
#' @return A list containing:
#' \describe{
#'   \item{coords}{Spatial coordinates matrix.}
#'   \item{regions}{Matrix or data frame of spatial regions per spot.}
#'   \item{deconv}{Cell type deconvolution matrix.}
#'   \item{motif_scores}{Matrix of motif activity scores per region.}
#' }
#'
#' @examples
#' \dontrun{
#' atac_data <- read_spatial_atac()
#' names(atac_data)
#' }
#'
#' @export
read_spatial_atac = function(){
  #read coordinates
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-ATAC/coord.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  coords <- as.matrix(readRDS(temp_file))


  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-ATAC/deconvolution.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  deconv <- as.matrix(readRDS(temp_file)$mat)



  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-ATAC/region_matrix.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  regions <- readRDS(temp_file)




  #get counts
  data = vector("list",3)
  for(j in c(1:3)){
    url <- paste0("https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-ATAC/scores_",j,".rds")
    # Temporary file to store the .rds
    temp_file <- tempfile(fileext = ".rds")

    # Download the file (use mode = "wb" for binary)
    download.file(url, destfile = temp_file, mode = "wb")

    data[[j]] = readRDS(temp_file)
  }
  data = do.call(cbind, data)



  return(list(coords = coords,regions = regions, deconv = deconv, motif_scores = data))

}



#' Read Spatial Long-Read RNA-Seq Data
#'
#' Loads spatial long-read RNA-seq data from a GitHub repository. Includes coordinates,
#' region annotations, deconvolution, library sizes, and expression matrices for genes and isoforms.
#'
#' @return A list containing:
#' \describe{
#'   \item{coords}{Matrix of spatial coordinates.}
#'   \item{regions}{Spatial region annotations.}
#'   \item{deconv}{Cell type deconvolution matrix.}
#'   \item{library_size}{Vector of library sizes per spot.}
#'   \item{total_gene_expression}{Matrix of total gene expression.}
#'   \item{isoform_expression}{Matrix of isoform-level expression.}
#' }
#'
#' @examples
#' \dontrun{
#' long_read_data <- read_spatial_long_read()
#' str(long_read_data)
#' }
#'
#' @export
read_spatial_long_read = function(){
  #read coordinates
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/coords.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  coords <- as.matrix(readRDS(temp_file))


  #read deconvolution
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/deconvolution.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  deconv <- as.matrix(readRDS(temp_file))



  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/region_matrix.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  regions <- readRDS(temp_file)


  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/library_size.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  library_size <- readRDS(temp_file)



  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/total_gene_expression.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  total_gene_expression <- readRDS(temp_file)




  #read effective niche
  url <- "https://raw.githubusercontent.com/kaishumason/SpotGLM-Example-Data/main/Spatial-Long-Read/isoform_expression.rds"
  # Temporary file to store the .rds
  temp_file <- tempfile(fileext = ".rds")

  # Download the file (use mode = "wb" for binary)
  download.file(url, destfile = temp_file, mode = "wb")

  # Read the RDS file
  isoform_expression <- readRDS(temp_file)






  return(list(coords = coords,regions = regions, deconv = deconv,library_size = library_size,
              total_gene_expression = total_gene_expression ,isoform_expression = isoform_expression))

}

