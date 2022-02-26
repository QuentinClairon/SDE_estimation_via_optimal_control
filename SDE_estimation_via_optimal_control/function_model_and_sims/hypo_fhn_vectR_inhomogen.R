hypo_fhn_vectR_inhomogen <- function(t,theta){
  dim_syst = 2
  vect_r= matrix(0,dim_syst,1)
  
  inv_eps = theta[1]
  gamma = theta[2]
  beta = theta[3]

  vect_r[2,1] = beta
  
  return(vect_r)
}
