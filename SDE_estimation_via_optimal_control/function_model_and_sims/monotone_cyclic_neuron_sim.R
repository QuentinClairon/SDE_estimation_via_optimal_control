monotone_cyclic_neuron_sim <- function(deb,h,fin,theta,y0,c_k){

  Times_obs = seq(deb,fin, by=h) 

dim_syst =  3;

nu_k = theta[1];


Mat_A =  matrix(rep(0,dim_syst^2),dim_syst,dim_syst)
Mat_A[1,1]= -nu_k; 
Mat_A[1,2]= 1;
Mat_A[2,2]= -nu_k;
Mat_A[2,3]= 1;
Mat_A[3,3]=-nu_k;

Res_sim =  matrix(rep(0,dim_syst*length(Times_obs)),dim_syst,length(Times_obs)) 
Stoch_part = matrix(rep(0,(dim_syst-2)*length(Times_obs)),dim_syst-2,length(Times_obs)) 
mat_b = matrix(c(0,0,c_k),dim_syst,dim_syst-2) 

Res_sim[,1] =  matrix(y0,dim_syst,1) 

for (i in  1 : (length(Times_obs)-1)){
  delta_i = Times_obs[i+1]-Times_obs[i]
  Stoch_part[,i+1]  =Stoch_part[,i]+sqrt(delta_i)*rnorm(dim_syst-2)
  Res_sim[,i+1] = (Mat_A*delta_i+diag(dim_syst))%*%Res_sim[,i] +mat_b%*%(Stoch_part[,i+1] -Stoch_part[,i])  
}

return(list(Times_obs,Res_sim))

}