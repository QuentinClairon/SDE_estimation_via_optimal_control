pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('deSolve')
library('optimx')

# Observation interval specification
  deb =0;
  delta = 0.02
  fin = 20
  
  # Parameter value specification
  g_l = 50;
  v_l = -70;
  v_e = 0;
  v_i = -80;
  val_I = -60;
  inv_tau_e = 1/0.5;
  inv_tau_i = 1/1;
  g_e = 17.8;
  g_i = 9.4;
  conduct = 10;
  
  theta_gen = c(g_l,v_l,v_e,v_i,val_I,inv_tau_e,inv_tau_i,g_e,g_i,conduct)
  theta_to_est = log(c(inv_tau_e,inv_tau_i,g_i))
  
  sigma_e = 0.1;
  sigma_i = 0.1;
  sigma = c(sigma_e, sigma_i)
  sigma_log = log(sigma)
  
  #Initial condition specification
  
  V0 = -60;
  Ge0 = 10;
  Gi0 = 1;
  
  y0 = c(V0,Ge0,Gi0)
  
  #Pseudo-linear representation with log transformation for volatility
  fun_mat_A = function(t,State,theta)syncon_matA_inhomogen(t,State,c(theta_gen[1:5],exp(theta[1:2]),theta_gen[8],exp(theta[3]),theta_gen[10]))
  fun_vect_R = function(t,theta)syncon_vectR_inhomogen(t,c(theta_gen[1:5],exp(theta[1:2]),theta_gen[8],exp(theta[3]),theta_gen[10]))
  fun_mat_Gamma = function(t,State,volat)syncon_matGamma_inhomogen(t,State,exp(volat))
  mat_C = matrix(c(1,0,0),1,3)
  
  
  #Initialisation of the optimization algorithm at true parameter value 
  param_ini = theta_to_est
  sigma_ini = sigma_log
  
  #Penalization hyperparameter selection
  weight_obs_trial = c(10^8,5*10^8,10^9,5*10^9)
  
  # Known initial condition
  x_0_known = y0
  
  #mB value as in section 2
  mB_val =1
  
  ################## Data generation  ################## 
  out_syncon = Synaptic_conductance_model_sim(deb,delta,fin,theta_gen,c(y0,1),sigma)
  
  Times_sim = out_syncon[[1]]
  Res_sim = out_syncon[[2]]

  
  Times_obs = Times_sim[seq(1,length(Times_sim),1)]
  Y_obs = Res_sim [1,seq(1,length(Times_sim),1)]
  Y_obs = matrix(Y_obs,1,length(Times_obs))
  
  #Initialisation for the state for the algorithm section 4.2.2
  State_ini = matrix(rep( y0,length(Times_obs)),length(y0),length(Times_obs))
  
  
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
    out_sde_est = est_param_sde_oca(Times_obs,Y_obs,State_ini ,param_ini,sigma_ini,fun_mat_A,  fun_vect_R ,fun_mat_Gamma,mat_C,weight_obs,x_0_known,mB_val)
    T2 = Sys.time()
    
  
    parameter_estimation = c(1/exp(out_sde_est$estimation[1]),1/exp(out_sde_est$estimation[2]),exp(out_sde_est$estimation[3]),exp(out_sde_est $estimation[4:5]))
    computation_time = difftime(T2, T1, units = "secs")
    Crit_extern = out_sde_est$Crit_extern

    print(parameter_estimation)
    list_res_est = append(list_res_est,list(list(parameter_estimation ,computation_time,Crit_extern)))

          },error=function(cond){
           return(1)}
       )
    

  }
  


