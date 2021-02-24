pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('deSolve')
library('optimx')

# Observation interval specification
  deb = 0
  fin = 10
  delta = 0.01
  
  # Parameter value specification
  eps =0.1
  gamma = 1.5
  beta = 0.8
  
  sigma = 0.3
  sigma_log = log(sigma)
  
  
  theta = c(eps,gamma,beta)
  theta_to_est = c(1/eps,gamma,beta)
  
  #Initial condition specification
  y0 = c(0,0)
  

  
  #Pseudo-linear representation with log transformation for volatility
  fun_mat_A = function(t,State,theta)hypo_fhn_matA_inhomogen(t,State,theta)
  fun_vectR = function(t,theta)hypo_fhn_vectR_inhomogen(t,theta)
  fun_mat_Gamma = function(t,State,volat){matrix(c(0,exp(volat)),2,1)}
  mat_C = matrix(c(1,0),1,2)
  
  #Initialisation of the optimization algorithm at true parameter value 
  param_ini = theta_to_est
  sigma_ini = sigma_log
  
  
  #Penalization hyperparameter selection
  weight_obs_trial = c(10^14,10^18,10^20)
  
  # Known initial condition
  x_0_known  = c()
  
  #mB value as in section 2
  mB_val =1
  
  
  ################## Data generation  ################## 
  out_hypofhn = fitzHugh_sde_sim(deb,delta,fin,theta,y0,sigma)
  
  Times_sim = out_hypofhn[[1]]
  Res_sim = out_hypofhn[[2]]
  
  Times_obs = Times_sim[seq(1,length(Times_sim),1)]
  Y_obs = Res_sim [1,seq(1,length(Times_sim),1)]
  Y_obs = matrix(Y_obs,1,length(Times_obs))
  
  #Initialisation for the state for the algorithm section 4.2.2
  State_ini = Res_sim [1:2,seq(1,length(Times_sim),1)]
  
  list_res_est = list()
  
  
  ################## Parameter estimation  ################## 
  # Estimation Results are embedded in list: list_res_est with the same length as  weight_obs_trial each entry is composed of
  # 1st entry: Estimation result
  # 2nd entry: Computation time
  # 3rd entry: log-transformed external criteria presented in section 4.3
  for (weight_obs in weight_obs_trial){
    
    parameter_estimation = c()
    computation_time = 10^20
    Crit_extern= 10^20
   
    
    out_try_catch <-tryCatch({
      
      
      
      T1 = Sys.time()
      out_sde_est = est_param_sde_oca(Times_obs,Y_obs,State_ini ,param_ini,sigma_ini,fun_mat_A, fun_vectR ,fun_mat_Gamma,mat_C,weight_obs,x_0_known,mB_val)
      T2 = Sys.time()
      
     
      parameter_estimation = c(1/out_sde_est$estimation[1],out_sde_est$estimation[2:3],exp(out_sde_est $estimation[4]))
      computation_time = difftime(T2, T1, units = "secs")
      Crit_extern = out_sde_est$Crit_extern

      
      print(parameter_estimation)

    },error=function(cond){
      return(1)}
    )
    
    list_res_est = append(list_res_est,list(list(parameter_estimation ,computation_time,Crit_extern)))
  
}

