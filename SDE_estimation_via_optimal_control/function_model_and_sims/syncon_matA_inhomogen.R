syncon_matA_inhomogen <- function(t,State,theta){
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
  
  res_at  =matrix(0,dim_syst,dim_syst)
  
  res_at[1,] = c(-g_l/conduct , (- State[1]+v_e)/conduct , (-State[1]+v_i)/conduct )
  res_at[2,] = c(0,-inv_tau_e,0)
  res_at[3,] = c(0,0,-inv_tau_i)

  
  return(res_at)
}