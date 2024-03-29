compute_Riccati_State_sde_oca <- function(Times_obs,Yn,State_i_m1,mat_A,vect_r,mat_Gamma,mat_C,param_cur,sigma_cur,inv_weight_obs ,x_0_known){
  ############################################################
  # DESCRIPTION: compute the Riccati  solution at a given iteration l of the algorithm described in section 4.2
  # 
  # FUNCTION SIGNATURE:  compute_Riccati_State_sde_oca <- function(Times_obs,Yn,State_i_m1,mat_A,vect_r,mat_Gamma,mat_C,param_cur,sigma_cur,inv_weight_obs ,x_0_known)
  # 
  # INPUT :
  # Times_obs                 :(vector) observation time points
  # Yn                        :(dim_obs x nb_obs sized matrix) observations matrix
  # State_i_m1                :(dim_syst x nb_obs sized matrix) computed optimal trajectory at iteration (l-1)
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
  # param_cur                 :(nb_param sized vector)  current value for theta 
  # sigma_cur                 :(nb_volatility sized vector) current value for volatility 
  # inv_weight_obs            :(positive real) inverse of the penalization parameter for Brownian increment
  # x0_known                  :(vector) list of known initial conditions, for the sake of convenience the model should be written such that the corresponding states variables are written last
  # OUTPUT:
  # Return a list set_Ricatti_State_opt_value composed of the following elements computed at iteration l:
  #       -set_Ricatti_State_opt_value$control_opt: The optimal control
  #       -set_Ricatti_State_opt_value$State_opt: The optimal trajectory 
  #       -set_Ricatti_State_opt_value$Ricatti: The Riccati  and adjoint equation solution 
  #       -set_Ricatti_State_opt_value$Val_sum_square_control: Residual sum of square of the control
  #       -set_Ricatti_State_opt_value$Crit_extern: External criteria value
  
  ############################################################
  
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_Gamma(Times_obs[1],State_i_m1[,1],sigma_cur))
  
  nb_obs = length(Times_obs)
  nb_x0_known = length(x_0_known)
  
  
  Y_fin= Yn[,nb_obs]
  mat_E_k = t(mat_C)%*%mat_C
  vect_h_k =  -t(mat_C)%*%Y_fin
  
  History_E  = list()
  History_E = list(mat_E_k)
  
  History_h  = list()
  History_h = list(vect_h_k)
  
  for (nk in 1:(nb_obs-1)){
    delta_k = Times_obs[nb_obs-nk+1] -Times_obs[nb_obs-nk]
    Yk= Yn[,nb_obs-nk] 
    
    
    State_pred_cur = State_i_m1[,nb_obs-nk]
    matAi =mat_A(Times_obs[nb_obs-nk], State_pred_cur,param_cur)
    vect_ri = vect_r(Times_obs[nb_obs-nk], param_cur)
    
    mat_A_par_discr = matAi*delta_k+diag(dim_syst)
    
    mat_Gammak= mat_Gamma(Times_obs[nb_obs-nk],State_pred_cur,sigma_cur)
    
    mat_Pk = diag(inv_weight_obs,dim_control,dim_control)
   
    mat_G = solve(mat_Pk+ delta_k*t(mat_Gammak)%*%mat_E_k%*%mat_Gammak)
    val_inter_k = t(mat_Gammak)%*%mat_E_k%*%mat_A_par_discr
    
    mat_E_k_1 = t(mat_C)%*%mat_C + t(mat_A_par_discr)%*%mat_E_k%*%mat_A_par_discr - delta_k*t(val_inter_k)%*% mat_G%*%val_inter_k
    
    vect_h_k_1 = -t(mat_C)%*%Yk+ delta_k*t(mat_A_par_discr)%*%mat_E_k%*%vect_ri+t(mat_A_par_discr)%*%vect_h_k
    vect_h_k_1 = vect_h_k_1- delta_k*t(val_inter_k)%*%mat_G%*%t(mat_Gammak)%*%(vect_h_k+delta_k*mat_E_k%*%vect_ri)
    mat_E_k = mat_E_k_1
    vect_h_k = vect_h_k_1
    History_E = append(list(mat_E_k),History_E)
    History_h = append(list( vect_h_k ),History_h)
  }
  
  mat_R0 =mat_E_k
  vect_h0 = vect_h_k
  
  if (nb_x0_known > 0){
    if (nb_x0_known < dim_syst){
      dim_unk = dim_syst - nb_x0_known;
      mat_R0_u = mat_R0[1:dim_unk,1:dim_unk]
      mat_R0_un = mat_R0[1:dim_unk,(dim_unk+1):dim_syst]                
      vect_h0_u = vect_h0[1:dim_unk]
      est_x0_u = -solve(mat_R0_u)%*%(mat_R0_un%*%x_0_known+vect_h0_u);
      est_x0 = c(est_x0_u,x_0_known)   
    }else{
    est_x0 = x_0_known
    }
  }else{
    est_x0 = -solve(mat_R0)%*%vect_h0
  }
  
 
  State_i = matrix(rep(0,dim_syst*nb_obs),dim_syst,nb_obs)
  State_i[,1] =est_x0
 
  control_i = matrix(rep(0,dim_control*(nb_obs-1)),dim_control,nb_obs-1)

  Val_sum_square_control = 0
  Crit_extern =0
  for (nkk in 1:(nb_obs-1)){
    delta_k = Times_obs[nkk+1] -Times_obs[nkk]
    Yk= Yn[,nkk] 
    
    State_pred_cur = State_i_m1[,nkk]
    matAi = mat_A(Times_obs[nkk],State_pred_cur,param_cur)
    vect_ri = vect_r(Times_obs[nkk], param_cur)
    mat_A_par_discr = matAi*delta_k+diag(dim_syst)
    
    mat_Gammak =  mat_Gamma(Times_obs[nkk],State_pred_cur,sigma_cur)
   
    mat_Pk = diag(inv_weight_obs,dim_control,dim_control)
    
    mat_E_k_p1 =  History_E[[nkk+1]]
    vect_h_k_p1=  History_h[[nkk+1]]
    
   
    control_ukk = -sqrt(delta_k)*solve(mat_Pk+ delta_k*t(mat_Gammak)%*%mat_E_k_p1%*%mat_Gammak)%*%t(mat_Gammak)%*%(mat_E_k_p1%*%(mat_A_par_discr%*%State_i[,nkk]+delta_k*vect_ri )+ vect_h_k_p1)
   
    control_i[,nkk] = control_ukk
    State_i[,nkk+1] = mat_A_par_discr%*%State_i[,nkk]+ delta_k*vect_ri +sqrt(delta_k)*mat_Gammak%*%control_ukk
    
    control_i_normalized_nkk = sum(control_ukk^2)
    if (control_i_normalized_nkk>0){
    Crit_extern =  Crit_extern + log((control_i_normalized_nkk^(dim_control/2-1))*exp(-control_i_normalized_nkk/2))
    }
     
    Val_sum_square_control = Val_sum_square_control+ sum(control_ukk^2)
  }
  

  
  return(list(Riccati = list(History_E,History_h),State_opt =State_i,control_opt = control_i ,Val_sum_square_control=Val_sum_square_control, Crit_extern = Crit_extern ))
}
