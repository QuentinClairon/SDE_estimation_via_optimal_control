fitzHugh_sde_sim <- function(deb,h,fin,theta,y0,sigma){

  Times_obs = seq(deb,fin, by=h) 
  nb_obs = length(Times_obs)
  
  dim_syst =  2
  
  eps = theta[1]
  gamma = theta[2]
  beta = theta[3]
  
  Res_sim = matrix(rep(0,dim_syst*nb_obs),dim_syst,nb_obs)
  Stoch_part = matrix(rep(0,(dim_syst-1)*nb_obs),dim_syst-1,nb_obs)
  seq_u = matrix(rep(0,(dim_syst-1)*(nb_obs-1)),dim_syst-1,nb_obs-1)
  
  mat_b = matrix(c(0,sigma),dim_syst,1)

  
  for (i in 1 : (nb_obs-1)){
  delta_i = Times_obs[i+1]-Times_obs[i]
 
  Stoch_part[,i+1] = Stoch_part[,i]+sqrt(delta_i)*rnorm(dim_syst-1);
  seq_u[,i]  = Stoch_part[,i+1]-Stoch_part[,i]
  
  part_det_i = matrix(c(0,0),dim_syst,1)
  part_det_i[1] =  (Res_sim[1,i] - Res_sim[1,i]^3 - Res_sim[2,i])*(delta_i/eps)
  part_det_i[2] =  (gamma*Res_sim[1,i] - Res_sim[2,i] +beta)*delta_i
  Res_sim[,i+1] = Res_sim[,i] +part_det_i+ mat_b*seq_u[,i]
  }
  return(list(Times_obs,Res_sim))

}
