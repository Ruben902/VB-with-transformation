This directory contains the set of MATLAB routines used for estimating the 
empirical examples in the paper:

Smith, Loaiza-Maya and Nott (2019) "High-dimensional Copula Variational Approximation through
Transformation".

Below you can find a brief description for the three folders in this directory.

If you have any related questions do not hesitate to email Ruben Loaiza-Maya at rloaizma@gmail.com.

Ruben Loaiza-Maya. 
12 November 2019.

/* ------------------------------ First Folder: Copula --------------------------------------------
This folder contains all the routines necessary to replicate the copula example in Section 3.2 of 
the paper.

Run script <S_VB_SKEW_copula_model.m> in the folder, to perform variational inference on the copula
using the Attempted Murder data.

Run script <Compare.m> in the folder to replicate Figures 1, 3 and 4 in the paper.

/* ------------------------------ Second Folder: Polypharm -----------------------------------------
This folder contains all the routines necessary to replicate the example in Section 4.2.1 of the paper.

Run script <polypharm_VB_SKEW.m> in the folder, to perform variational inference in the 
Polypharm data using a mixed logistic regression.

Run script <CompareResults.m> in the folder to replicate Figures 5 and 6 in the paper. It also produces
Table 2.

Go to <Polypharm\INLA> and then run script <YJ_PAR_densities_plus_MF_INLA.m> to add INLA to Figure 6 in 
the paper. 

/* ------------------------------ Third Folder: Spam_KRKP_iono_Mush ---------------------------------
This folder contains all the routines necessary to replicate the copula example in Section 4.2.2 of 
the paper.

Run script <S_VB_SKEW.m> in the folder, to perform variational inference on the logistic regression
using the spam, KRKP, ionosphere or mushroom data. 

Run script <Compare_ELBO_times.m> in the folder to replicate Figure 7 in the paper.
Run script <Compare_ELBO_times.m> in the folder to replicate Table 3 in the paper.


------------------------------------------------------------------------------------------------------
Note: All the folders mentioned above contain their own README files with more details about them.