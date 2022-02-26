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
  nu_k = 0.20
  c_k = 0.15
  
  theta = c(nu_k)
  sigma_log = log(c_k)
  
  #Initial condition specification
  x0 = c(0,0,0)

  #Pseudo-linear representation with log transformation for volatility
  fun_mat_A = function(t,State,theta)monotone_cyclic_neuron_matA(t,theta)
  fun_vect_R = function(t,theta){matrix(0,3,1)}
  fun_mat_Gamma = function(t,State,volat){matrix(c(0,0,exp(volat)),3,1)}
  mat_C = matrix(c(1,0,0),1,3)
  
  
  #Initialisation of the optimization algorithm at true parameter value 
  param_ini = theta
  sigma_ini = sigma_log
  
  
  #Penalization hyperparameter selection
  weight_obs_trial = c(10^15,10^20,10^25,10^30)
  
  # Known initial condition
  x_0_known = x0 
  
  #mB value as in Section 2
  mB_val =2
  
  
  ################## Data generation  ################## 
  out_cyclic= monotone_cyclic_neuron_sim(deb,delta,fin,theta,x0,c_k)
  
  Times_sim = out_cyclic[[1]]
  Res_sim = out_cyclic[[2]]
  
  Times_obs = Times_sim[seq(1,length(Times_sim),1)]
  Y_obs = Res_sim [1,seq(1,length(Times_sim),1)]
  Y_obs = matrix(Y_obs,1,length(Times_obs))
 
  list_res_est = list()

  ################## Parameter estimation  ################## 
  # Estimation Results are embedded in list: list_res_est with the same length as  weight_obs_trial each entry is composed of
  # 1st entry: Estimation result
  # 2nd entry: Computation time
  # 3rd entry: log-transformed external criteria presented in section 4.3
  for (weight_obs in weight_obs_trial){
    
    parameter_estimation = c()
    computation_time = 10^20
    Crit_extern = 10^20
 
    out_try_catch <-tryCatch({
      
      T1 = Sys.time()
      out_sde_est = est_param_sde_oca(Times_obs,Y_obs,State_ini=c()  ,param_ini,sigma_ini,fun_mat_A,fun_vect_R ,fun_mat_Gamma,mat_C,weight_obs,x_0_known,mB_val,
                                      linear_equation=1)
      T2 = Sys.time()
      
      parameter_estimation = c(out_sde_est$estimation[1],exp(out_sde_est $estimation[2]))
      computation_time = difftime(T2, T1, units = "secs")
      Crit_extern = out_sde_est$Crit_extern
      
      print(parameter_estimation)
    },error=function(cond){
     return(1)}
    )
    
     list_res_est = append(list_res_est,list(list(parameter_estimation ,computation_time,Crit_extern)))
    
}

