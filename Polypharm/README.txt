This MATLAB software implements Variational approximation
through transformations for the polypharm dataset.

You can freely use this code for academic purposes as long as 
you refer to the publication below.

Smith, Loaiza-Maya and Nott (2019) "High-dimensional Copula Variational Approximation through
Transformation".

If you have any related questions do not hesitate to email Ruben Loaiza-Maya at rloaizma@gmail.com.

Ruben Loaiza-Maya. 
12 November 2019.

### Main programs : polypharm_VB_SKEW and VBtransf

From this script you can select the number of factor in the covariance, the 
type of copula used in the approximation, and the type of transformation.
For more details resort to the web-appendix.

For the approximations based on the iGH transformation 40000 VB steps were 
used. The remaining approximations were obtained with 20000 VB steps

### Note ###
The subfolder Derivatives contains the routines needed to compute the terms
in Tables 1 and 4 in the paper. 
The key function is <dtheta_dBDelta.m> which is coded to resemble Table 4
in the paper.