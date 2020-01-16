close all
clear;
clear dtheta_dBDelta;   %Clears persistent variables in the function
rng(19900209);
workingDir = pwd;
addpath(genpath([workingDir '/Derivatives']));
%Defining the type of approximation
Skewlab = {'NoSkew' 'Skew'};    
p = 5;                                                                     % Number of factors on the covariance matrix of the approximation
SkewNorm = 0;                                                              % Set to zero if skew-normal will not be used. Set to one if skew-normal is to be used
Transf = 'YJ';                                                             % Type of transformation. For Yeo-Johnson transformation set to 'YJ'. For GH transformation set to 'GH', for no transformation set to '' or 'None'
niter = 40000;                                                              % Number of total iterations for the VB algorithm
%------ Defining log-posteriors and derivatives of the model---------------
load('polypharm.mat')
b_n = 9; %first index of the random effect
[Ntr, ~] = size(Xtr);
q = size(Xtr,2) + 1;
g = @(theta) log_logreg(theta(1:(end-1)), Xtr,Ytr,(Ytr+1)/2,false) + log_glmm(theta, zeros(q,1), 100*ones(q,1), b_n,false);            %Creating function to evaluate log-posterior
delta_logh = @(theta) [log_logreg(theta(1:(end-1)), Xtr,Ytr,(Ytr+1)/2,true);0] + log_glmm(theta, zeros(q,1), 100*ones(q,1), b_n,true); %Creating function to evaluate gradient of log-posterior
%--------------------------------------------------------------------------
[F_mu,F_d,F_B,eta,alphaphi,F_LB,StoreTime] = VBtransf(g,delta_logh,q,p,Transf,SkewNorm,niter,100,false);
if p==0
    save(['Poly_' Skewlab{SkewNorm+1} Transf 'VB_MF.mat'])
else
    save(['Poly_' Skewlab{SkewNorm+1} Transf 'VB.mat'])
end
%------ Visualising approximate posteriors --------------------------------
ind = [1 100 200 300]; %Indexes of approximate posteriors to plot                                      
plot_posterior(ind,eta,F_mu,F_B,F_d,alphaphi,Transf)
