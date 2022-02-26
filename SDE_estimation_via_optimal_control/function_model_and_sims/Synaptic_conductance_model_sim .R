Synaptic_conductance_model_sim <- function(deb,h,fin,theta,y0,sigma){

Times_obs = seq(deb,fin, by=h) 
nb_obs = length(Times_obs)

dim_syst = 4
dim_u = 2

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

matA_t <-function (t,y){
  res_at  =matrix(rep(0,dim_syst^2),dim_syst,dim_syst) 
  res_at[1,] = c(-g_l/conduct , (- y[1]+v_e)/conduct , (-y[1]+v_i)/conduct , (g_l*v_l+val_I)/conduct )
  res_at[2,] = c(0,-inv_tau_e,0,inv_tau_e*g_e)
  res_at[3,] = c(0,0,-inv_tau_i,inv_tau_i*g_i)
  res_at[4,] = c(0,0,0,0) 
  
  return(res_at)
}

mat_b_f <-function (y,sig){
  res_bt  =matrix(rep(0,dim_syst*dim_u),dim_syst,dim_u)
  res_bt[2,] = c(sig[1]*sqrt(y[2]),0)
  res_bt[3,] = c(0,sig[2]*sqrt(y[3]))
 
 return(res_bt)
}

Res_sim = matrix(rep(0,dim_syst*nb_obs),dim_syst,nb_obs);
Res_sim[,1] = matrix(y0,dim_syst,1);

for (ni in 1 : (nb_obs-1)){
  delta_ni = Times_obs[ni+1]-Times_obs[ni]
  Res_sim[,ni+1]= (diag(dim_syst)+delta_ni*matA_t(Times_obs[ni],Res_sim[,ni]))%*%Res_sim[,ni]
  Res_sim[,ni+1] =Res_sim[,ni+1] + mat_b_f(Res_sim[,ni],sigma)%*%(sqrt(delta_ni)*rnorm(2))
}

return(list(Times_obs,Res_sim))
}
