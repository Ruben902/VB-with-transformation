This MATLAB set of routines implements Variational approximation
through transformations for the Spam, ionosphere, KRKP & mushroom dataset.

You can freely use this code for academic purposes as long as 
you refer to the publication below.

M. S. Smith, R. A. Loaiza-Maya and D. J. Nott (2019). High-dimensional Copula 
Variational Approximation through Transformation.  

If you have any related questions do not hesitate to email Ruben Loaiza-Maya at rloaizma@gmail.com.

Ruben Loaiza-Maya. 
12 November 2019.

### Main programs : S_VB_SKEW and VBtransf

To select from the four alternative data sets change the value of the vari-
able <data_pick> in the script file <S_VB_SKEW>.
For more details resort to the web-appendix.


### Acknowledgement ###

We would like to also thank the following authors for providing their 
Matlab codes online.

1) M. Schmidt for the MCMC codes to run Bayesian Logistic Regression
    (Most files inside MCMC folder)

2) V. M.-H. Ong, D. J. Nott, and M. S. Smith for providing their routines 
    for calibrating the variational approximation with factor covariance 
     structure.

### Note ###
The subfolder Derivatives contains the routines needed to compute the terms
in Tables 1 and 4 in the paper. 
The key function is <dtheta_dBDelta.m> which is coded to resemble Table 4
in the paper.