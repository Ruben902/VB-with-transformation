%-- Using VBDA on a univariate time series; only implemented for YJ transformation --%
clear grad_lambda_logq                                                     %Be aware of persistent variables in this function
clear
family = 'gumbel_mix';    %Copula family for pair-copula components
p = 3;                    %D-Vine Markov order
k = 3;                    %Number of factors for covariance matrix
S = 500;                  %Number of evaluations to compute gradient
nVB = 100;               %Number of VB steps
y = pickdata('AttemptedMurder');    %Ordinal time series, change to your own data
mdist = selectedDist(y);     %Marginal distribution of y                  %
VA = 3;                   %Variational Bayes approach for estimation  
filename = ['wd_AttemptedMurder_p' mat2str(p) '_k' mat2str(k) '_' family '_VA' mat2str(VA) '.mat'];
VBDAobj = VBDAfit_uni(y,mdist,family,p,k,S,nVB,VA,filename);%Using VBDA
VBDAplots(VBDAobj)        %Producing plots from object VBDAobj
