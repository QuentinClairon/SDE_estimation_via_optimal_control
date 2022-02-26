est_param_sde_oca<- function(Times_obs,Yn,State_ini,param_ini,sigma_ini=c(),mat_A,vect_r,mat_Gamma,mat_C,weight_obs,x_0_known,mB_val,linear_equation=0){
  ############################################################
  # DESCRIPTION: compute the parameter estimation via the method presented for SDE models
  # 
  # FUNCTION SIGNATURE:  est_param_sde_oca<- function(Times_obs,Yn,State_ini,param_ini,sigma_ini=c(),mat_A,vect_r,mat_Gamma,mat_C,weight_obs,x_0_known,mB_val,linear_equation=0)
  # 
  # INPUT :
  # Times_obs                 :(vector) observation time points
  # Yn                        :(dim_obs x nb_obs sized matrix) observations matrix
  # State_ini                 :(dim_syst x nb_obs sized matrix) initial guess for the system state
  # param_ini                 :(nb_param sized vector)  initial guess for theta estimator
  # sigma_ini                 :(nb_volatility sized vector) initial guess for volatility estimator
  #                              should be written such that the corresponding states variables are written last
  # mat_A                     :(function handler) matrix valued function representing  A in the pseudolinear model dZ =  A(t,Z,theta)Zdt +r(t,theta)dt +Gamma(t,Z,sigma)dWt of signature: 
  #                                              mat_A<- function(t,State,param) 
  #                                                INPUT: t (real number) current time value
  #                                                     : State (dim_syst sized vector) state_variable value 
  #                                                     : param (vector) current parameter value
  #                                                OUTPUT: (dim_syst square matrix) the value of  A(t,State,param) 
  # vect_r                     :(function handler) vector valued function representing  r in the pseudolinear model dZ =  A(t,Z,theta)Zdt +r(t,theta)dt +Gamma(t,Z,sigma)dWt of signature: 
  #                                              vect_r- function(t,param) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param (vector) current parameter value
  #                                                OUTPUT: (dim_syst x 1 matrix) the value of  r(t,param) 
  # mat_Gamma                 :(function handler) matrix valued function representing  Gamma in the pseudolinear model dZ =  A(t,Z,theta)Zdt +r(t,theta)dt +Gamma(t,Z,sigma)dWt of signature: 
  #                                              mat_Gamma<- function(t,State,sigma) 
  #                                                INPUT: t (real number) current time value
  #                                                     : State (dim_syst sized vector) state_variable value 
  #                                                     : sigma (vector) current volatility value
  #                                                OUTPUT: (dim_syst x dim_rough sized matrix) the value of  Gamma(t,State,sigma) 
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the state variables
  # weight_obs                :(positive real) penalization parameter for Brownian increment
  # x0_known                  :(vector) list of known initial conditions, for the sake of convenience the model 
  # mB_val                    :(positive integer) Defined in section 2, maximal length of the  minimal paths connecting  the rough components to the all the smooth ones
  # linear_equation           :(integer, default= 0) indicate if the SDE is linear (when linear_equation =1) or not. Not mandatory, but allows to avoid unnecessary computation for linear SDEs. 
  # OUTPUT:
  # Return a list  res_est  composed of:
  # -res_est$estimation: the estimated value of theta and volatility
  # -res_est$State_opt: The optimal trajectory corresponding to the estimation 
  # -res_est$control_opt: The optimal control  corresponding to the estimation
  # -res_est$Crit_extern: External criteria value
  
  ############################################################
  mB_val_p1 = mB_val+1
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_Gamma(Times_obs[1],State_ini[,1],sigma_ini))
  
  inv_weight_obs = 1/weight_obs
  
  nb_obs = length(Times_obs)
  Yn_fin = Yn[,nb_obs];
  
  nb_param_struct = length(param_ini);
  
  nb_x0_known = length(x_0_known);
  
  if (nb_x0_known>0){
    x_0_known = matrix(x_0_known,nb_x0_known,1) 
  }
  
  func_eco_E <- function(param_t){
    
    if (length(sigma_ini) ==0){
    param_cur = param_t
    sigma_cur = c()
    }else{
      param_cur = param_t[1:nb_param_struct];
      sigma_cur = param_t[(nb_param_struct+1):length(param_t)]
    }

    State_i_m1 = State_ini;
    val_err_iter_cur = 10^100;
    stop_criteria = 0;
    nb_iter = 0
    nb_iter_max =50
    if (linear_equation==1){
        nb_iter_max = 1
    }

    
    while (stop_criteria ==0){
      ensemble_Ricatti_comp =   compute_Riccati_State_sde_oca(Times_obs,Yn,State_i_m1,mat_A,vect_r,mat_Gamma,mat_C,param_cur,sigma_cur,inv_weight_obs ,x_0_known)
      
      val_err_iter = sum((ensemble_Ricatti_comp$State_opt -  State_i_m1)^2)
      State_i_m1 = ensemble_Ricatti_comp$State_opt
      nb_iter = nb_iter+1
      
      
      if (val_err_iter <10^-4){
        stop_criteria = 1;
      }
      
      if (nb_iter  >= nb_iter_max){
        stop_criteria = 1;
      }
      
    }
    
    State_finale = ensemble_Ricatti_comp$State_opt
    Val_sum_square_control = ensemble_Ricatti_comp$Val_sum_square_control
    
    fval_seq = matrix(rep(0,nb_obs-mB_val_p1),1,nb_obs-mB_val_p1)
    for (nkk in mB_val_p1: (nb_obs-1)){
      Ykk_pmB= Yn[,nkk+1]
      
      Val_product_A = diag(dim_syst);
      Val_State_j = State_finale[,nkk-mB_val_p1+1]
      Val_Eps =matrix(rep(0,dim_obs^2),dim_obs,dim_obs);
      mat_sum_F = matrix(0,dim_obs,1)
      
      for (r in 1 : mB_val_p1){
        delta_k = Times_obs[nkk-r+2] -Times_obs[nkk-r+1]
        
        State_pred_cur = State_finale[,nkk-r+1]
        matAi = mat_A(Times_obs[nkk-r+1],State_pred_cur,param_cur)
        
        mat_A_par_discr = matAi*delta_k+diag(dim_syst)
        
        mat_F_cur =  delta_k*mat_C%*%Val_product_A%*%vect_r(Times_obs[nkk-r+1],param_cur)  
        mat_G_cur = sqrt(delta_k)*mat_C%*%Val_product_A%*%mat_Gamma(Times_obs[nkk-r+1],State_pred_cur,sigma_cur)                     
        Val_Eps = Val_Eps+mat_G_cur%*%t(mat_G_cur)
        
        mat_sum_F =  mat_sum_F+ mat_F_cur
                                     
       Val_product_A = Val_product_A%*%mat_A_par_discr;
      }
   
       X_bar =  Ykk_pmB -mat_C%*%Val_product_A%*%Val_State_j - mat_sum_F 
                  
      
       fval_seq[nkk] =  log(det(Val_Eps)) + t(X_bar)%*%solve(Val_Eps)%*%X_bar;
    }
    fval = sum(fval_seq)
    
   
    return(list(ensemble_Ricatti_comp,fval= fval))
  }
  
  func_eco <- function(param){
    ensemble_Ricatti_comp = func_eco_E(param)
    return(ensemble_Ricatti_comp$fval)
  }
  
  if (length(sigma_ini) ==0){
    param_ini_t = param_ini;
  }else{
    param_ini_t = c(param_ini,sigma_ini)
  }
 

  res_estim = optimr(param_ini_t, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,method="Nelder-Mead",control = list(trace=1))
  ensemble_Ricatti_final = func_eco_E(res_estim$par)
  State_est_fin = ensemble_Ricatti_final[[1]]$State_opt
  control_seq_fin = ensemble_Ricatti_final[[1]]$control_opt
  Crit_extern = ensemble_Ricatti_final[[1]]$Crit_extern
  
  residual_est = sum((Yn - mat_C%*%State_est_fin[1:dim_syst,])^2)
  Sum_square_control = matrix(rep(0,dim_control),dim_control,1)
  for (nnkk in 1: (nb_obs-1)){
    delta_k = Times_obs[nnkk+1] -Times_obs[nnkk]               
    Sum_square_control = control_seq_fin[,nnkk]^2 + Sum_square_control  
  }
     
  
  law_large_number_khi2 = sum((Sum_square_control - (nb_obs-1))^2)
  
  return(list(estimation = res_estim$par,law_large_number_khi2=law_large_number_khi2,residual_est  =residual_est,State_opt=State_est_fin,control_opt= control_seq_fin,Crit_extern=Crit_extern))
}
