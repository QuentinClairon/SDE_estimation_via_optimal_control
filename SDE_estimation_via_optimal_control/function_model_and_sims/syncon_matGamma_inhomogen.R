syncon_matGamma_inhomogen <-function (t,State,volat){
  dim_syst = 3
  dim_u = 2
  
  res_bt  =matrix(rep(0,dim_syst*dim_u),dim_syst,dim_u)
  res_bt[2,] = c(volat[1]*sqrt(State[2]),0)
  res_bt[3,] = c(0,volat[2]*sqrt(State[3]))
  
  return(res_bt)
}
