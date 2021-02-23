syncon_vectR_inhomogen <- function(t,theta){
  dim_syst = 3
  
  
  g_l = theta[1]
  v_l = theta[2]
  v_e = theta[3]
  v_i = theta[4]
  val_I = theta[5]
  inv_tau_e = theta[6]
  inv_tau_i = theta[7]
  g_e = theta[8]
  g_i = theta[9]
  
  conduct = theta[10]
  
  res_rt  =matrix(0,dim_syst,1)
  
  res_rt[1,1] =  (g_l*v_l+val_I)/conduct
  res_rt[2,1] = inv_tau_e*g_e
  res_rt[3,1] = inv_tau_i*g_i

  
  return(res_rt)
}