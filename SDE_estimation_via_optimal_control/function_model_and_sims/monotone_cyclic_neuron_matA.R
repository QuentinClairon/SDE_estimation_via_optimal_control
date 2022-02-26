monotone_cyclic_neuron_matA <- function(t,theta){
dim_syst = 3
nu_k = theta[1]

Mat_A = matrix(rep(0,dim_syst^2),dim_syst,dim_syst)

Mat_A[1,1] = -nu_k
Mat_A[1,2] = 1
Mat_A[2,2] = -nu_k
Mat_A[2,3] = 1
Mat_A[3,3] = -nu_k

return(Mat_A)
}
