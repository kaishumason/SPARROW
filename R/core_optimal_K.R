#step 1: get global min standard error (make signal/noise at least 3 for smallest effects)
min_effect = abs(effect_scale/5)
global_min_SE = SE_optimizer(min_effect,pow = 0.999,alpha = 0.05, acc = 1e-3)

#step 2: Get min SE 
A_full = solve(t(big_X)%*%big_X)
SE_full = sqrt(diag(A_full))

nSE = length(SE_full)
ct_min_SE = rep(NA,nSE)
for(j in c(1:nSE)){
  ct_min_power = power_calc(SE_full[j],min_effect,0.05)
  #if we have full power for smallest effect 
  if(ct_min_power > (0.999)){
    ct_min_SE[j] = global_min_SE
  }else if(ct_min_power < 0.5){
    #if we have too little power, find new point 
    ct_min_effect = min_effect_optimizer(SE_full[j],0.80,0.05,acc = 1e-3)
    ct_min_SE[j] = SE_optimizer(ct_min_effect,pow = 0.79,alpha = 0.05, acc = 1e-3)
  }else{
    #otherwise make sure we have enough power
    ct_min_SE[j] = SE_optimizer(min_effect,pow = ct_min_power * 0.975,alpha = 0.05, acc = 1e-3)
  }
}

#get cell type specifici global SE
ct_min_SE[ct_min_SE < global_min_SE] = global_min_SE