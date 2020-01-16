%%% This function computes a variational approximation to a posterior 
%%% density. The approximation is defined through a transformation as in
%%% M. S. Smith, R. A. Loaiza-Maya, and D. J. Nott (2019). High-dimensional Copula 
%%% Variational Approximation through Transformation. 
%%% The arguments of the function are:
%logpost:   A function that evalutes the log-posterior of the model
%dlogpost:  A function provides the gradients of the log posterior with
%           respect to the parameters of the model
%q:         The total number of parameters of the model
%p:         The number of factors for the factor structure of the
%           approximation
%Transf:    The type of transformation to be used for the approximation
%           Set to 'None', 'GH' or 'YJ'
%SkewNorm:  Indicates if skew-normal distribution is to be used. 1 for true
%           0 otherwise
%niter:     Total number of iterations for the VB algorithm
%time_inter:Determines how often iteration time is measured
%do_plot:   Indicates if LB plot should be produced
%%% The outputs of the function are the VB calibrated
%F_mu:      Location parameter of the VB approximation
%F_d:       Diagonal parameter values of D  in represenation BB'+D^2
%F_B:       Parameter values in factor matrix B of represenation BB'+D^2
%F_eta:     Parameter values of the transformations used in the
%           approximation
%F_alphaphi:Asymmetry parameter of the skew-normal distribution
function [F_mu,F_d,F_B,F_eta,F_alphaphi,F_LB,StoreTime] = VBtransf(logpost,dlogpost,q,p,Transf,SkewNorm,niter,time_inter,do_plot)
%time_inter = 100 ; do_plot = true;
%%% Initializing variational parameters starting point
B = zeros(q,p) + 0.001;
B(logical(eye(size(B)))) = 0;
B(find(triu(B,1))) = 0;
mu = sparse(q,1);
d = ones(q,1)*0.1;
alphaphi = sparse(q,1);
switch Transf
    case 'YJ'
        eta = ones(q,1);
    case 'GH'
        eta = zeros(q,2);
        eta(:,1)=0.001;
        eta(:,2)=0.1;
    case 'None'
        eta = ones(q,1);
    case ''
        eta = ones(q,1);
    case 'iGH'
        eta = zeros(q,2);
        eta(:,1)=0.001;
        eta(:,2)=0.01;
end
%%% ADADELTA Learning rate parameters
Edelta2_mu = zeros(length(mu),1);
Eg2_mu = zeros(length(mu),1);
Edelta2_B = zeros(length(B(:)),1);
Eg2_B = zeros(length(B(:)),1);
Edelta2_d = zeros(length(d),1);
Eg2_d = zeros(length(d),1);
Edelta2_alphaphi = zeros(length(alphaphi),1);
Eg2_alphaphi = zeros(length(alphaphi),1);
Edelta2_eta = zeros(length(eta(:)),1);
Eg2_eta = zeros(length(eta(:)),1);
Edelta2_tau = zeros(length(eta(:)),1);
Eg2_tau = zeros(length(eta(:)),1);
ADA.rho = 0.95;
ADA.eps_step = 10^-6;
ADA.Edelta2_mu = Edelta2_mu;
ADA.Eg2_mu = Eg2_mu;
ADA.Edelta2_B = Edelta2_B;
ADA.Eg2_B = Eg2_B;
ADA.Edelta2_d = Edelta2_d;
ADA.Eg2_d = Eg2_d;
ADA.Edelta2_alphaphi = Edelta2_alphaphi;
ADA.Eg2_alphaphi = Eg2_alphaphi;
ADA.Edelta2_eta = Edelta2_eta;
ADA.Eg2_eta = Eg2_eta;
ADA.Edelta2_tau = Edelta2_tau;
ADA.Eg2_tau = Eg2_tau;
%%%----------------------------------------

Dinv2B = bsxfun(@times,B,1./d.^2);
Siginv = diag(1./d.^2) - Dinv2B/(speye(p) + B'*Dinv2B)*Dinv2B';

StoreLB = zeros(niter,1);
StoreTime = zeros(niter,1);
mu_sum = zeros(size(mu));
d_sum = zeros(size(d));
B_sum = zeros(size(B));
alphaphi_sum = zeros(size(alphaphi));
eta_sum = zeros(size(eta));
n_store = 1000;
if niter<n_store
    n_store = 1;
end
tic
for iter = 1:niter
    [LowerB,B,mu,d,eta,alphaphi,ADA,Siginv] = VB_step(B, mu, d, eta, alphaphi, logpost,dlogpost,Siginv,ADA,p,SkewNorm,Transf);
    StoreLB(iter) = LowerB;
    StoreTime(iter) = toc;
    if ~(mod(iter,time_inter))
        disp(['Iter = ' mat2str(iter) '; LB = ' mat2str(round(StoreLB(iter)*100)/100) '; Time = ' mat2str(round(toc*100)/100)])
        if do_plot
            plot(StoreLB(1:iter),'r')
            drawnow
        end
    end
    if iter >(niter-n_store)
        mu_sum = mu_sum + mu;
        d_sum = d_sum + d;
        B_sum = B_sum + B;
        alphaphi_sum = alphaphi_sum + alphaphi;
        eta_sum = eta_sum + eta;
    end
end

F_mu = mu_sum/n_store;
F_d = d_sum/n_store;
F_B = B_sum/n_store;
F_eta = eta_sum/n_store;
F_alphaphi = alphaphi_sum/n_store;

F_LB = StoreLB;
