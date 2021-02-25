The code contained in this folder is here to reproduce the simulated section results given by the optimal control based methods presented in the article
 "OPTIMAL CONTROL FOR PARAMETER ESTIMATION IN PARTIALLY OBSERVED ELLIPTIC AND HYPOELLIPTIC STOCHASTIC DIFFERENTIAL EQUATIONS".

For the code to work you simply need to set this folder as the working directory. It is composed of several scripts and two folders. 

Script:
The scripts simulate data set according to the models and experimental designs presented in "Section: 5 Experiments"
and then proceed to estimation by using the optimal control based methods with the same hyperparameters as in the article:
1) Scripts "script_test_cyclic_neuron_estimation"  reproduces experimental designs and optimal control based estimation procedures of section 5.2.1
2) Scripts "script_test_FHN_estimation" reproduces experimental designs and optimal control based estimation procedures of section 5.2.2
3) Scripts "script_test_synaptic_conductance_estimation" reproduces experimental designs and optimal control based estimation procedures of section 5.2.3

Folders:
1) The folder "function_model_and_sims" contains all the needed model specific functions to generate the observations and estimate parameters from the previous scripts.

2) The folder "function_util_estimation" contains the generic functions for estimating parameters via our optimal control based approaches for any SDE models. 
They are commented in order to be reusable for other practitioner models 