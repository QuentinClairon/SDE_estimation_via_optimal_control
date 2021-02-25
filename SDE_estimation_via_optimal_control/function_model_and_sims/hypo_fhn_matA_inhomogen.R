hypo_fhn_matA_inhomogen <- function(t,State,theta){
  dim_syst = 2
  
  Mat_A = matrix(0,dim_syst,dim_syst)
  
  inv_eps = theta[1]
  gamma = theta[2]
  beta = theta[3]
  
  Mat_A[1,1] = (1-State[1]^2)*inv_eps
  Mat_A[1,2] = -inv_eps    
 # Mat_A[1,3] = 0
  
  Mat_A[2,1] = gamma
  Mat_A[2,2] = -1
  #  Mat_A[2,3] = beta
  
  return(Mat_A)
}