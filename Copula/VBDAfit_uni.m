function VBDAobj = VBDAfit_uni(y,mdist,family,p,k,S,nVB,VA,filename)
%-------------------------------------------------------------------------
%%% Input arguments:
%--- y:      Sample data vector.
%--- mdist:  Marginal distribution of y (must be ordinal). Distribution
%            object for which functions <cdf> and <icdf> can be evaluated.
%--- family: Bivariate copula family to use for pair-copulas.
%--- p:      Markov order of D-Vine copula
%--- k:      Number of factors in the covariance matrix of the normal
%            approximation.
%--- S:      Number of evaluations of the D-Vine copula used to compute the
%            gradient needed for VB estimation.
%--- nVB:    Number of VB steps used for inference.
%--- VA:     Variational Bayes approximation to the distribution of u
%-------------------------------------------------------------------------
%%% Output arguments:
%--- VBDAobj.gamma_mean: Variational posterior mean of copula parameters.
%--- VBDAobj.gamma_sd:   Variational posterior standard deviation of copula
%                        parameters.
%--- VBDAobj.mu:         Mean vector of approximating density q, on the
%                        transformed copula parameters.
%--- VBDAobj.B:          Factors of the covariance matrix of approximating
%                        density q, on the transformed copula parameters
%                        (Sigma = B*B'+D.^2).
%--- VBDAobj.D:          Idiosyncratic terms of the covariance matrix of
%                        approximating density q, on the transformed copula
%                        parameters (Sigma = B*B'+D.^2).
%--- VBDAobj.eta:        Calibrated variational parameters of the transfor-
%                        for the copula parameters.
%--- VBDAobj.muz         If VA = 2, or VA = 3, is the mean of the variatio-
%                        nal approximation to z = Phi^(-1)(u).
%--- VBDAobj.logsigmaz   If VA = 2, is log(sigma) of the variatio-nal appr-
%                        oximation to z = Phi^(-1)(u).
%--- VBDAobj.C           If VA = 3, is the matrix for which Sigma^(-1)=CC',
%                        where Sigma = var(z).
%--- VBDAobj.etau:       Calibrated variational parameters of the transfor-
%                        for auxilary copula data used for augmentation.
%--- VBDAobj.lambda:     Matrix with alternative reparametrization of mu
%                        and Sigma
%--- VBDAobj.LB:         Lower Bound values at each VB step.
%--- VBDAobj.dist:       CDF object for the ordinal data in y.
%--- VBDAobj.y:          Sample data vector.
%--- VBDAobj.family:     Bivariate copula family to use for pair-copulas.
%--- VBDAobj.p:          Markov order of D-Vine copula
%--- VBDAobj.k:          Number of factors in the covariace matrix of the
%                        normal approximation.
%--- VBDAobj.S:          Number of evaluations of the D-Vine copula used to
%                        compute the gradient needed for VB estimation.
%--- VBDAobj.nVB:        Number of VB steps used for inference.
%--- VBDAobj.VA:        Number of VB steps used for inference.
%-------------------------------------------------------------------------
disp('----> VB estimation underway')
[ai,bi,n,LB,ct,lambda,t,mu,~,muz,logsigmaz,C,tau,tauu,T,B,D,Transf] = VBinival(y,mdist,family,p,k,S,nVB,VA);
ADA = [];tic;
while(t<nVB)
    t = t+1;
    [theta,Gamma] = q_lambda(mu,S,family,B,D,tau2eta(tau,Transf),Transf);
    [u,z,vareps] = q_u(ai,bi,muz,logsigmaz,C,S,VA,tau2eta(tauu,Transf),Transf);
    [H,LogLhat,log_theta,log_q_lambda,grad_lambda] = gradient_compute(theta,Gamma,u,z,vareps,ai,bi,family,lambda,k,VA,ct,Transf);
    ct = ct_constant(LogLhat,log_theta,log_q_lambda,grad_lambda);
    [Change_delta,ADA] = ADADELTA(ADA,H);
    lambda = lambda + Change_delta;
    ind_d = (n+k*n-k*(k-1)*0.5+1):(n+k*n-k*(k-1)*0.5+n);
    lambda_d = lambda(ind_d);
    lambda(ind_d)=lambda_d.*(abs(lambda_d)>=1e-4)+1e-4.*(abs(lambda_d)<1e-4);   %Lower Limit for d
    [mu,~,B,D,muz,logsigmaz,C,tau,tauu] = lambda2musigma(lambda,n,T,k,VA,Transf);
    LB(t) = mean(LogLhat+log_theta-log_q_lambda);
    if mod(t,10)==0
        disp(['      VB step = ' mat2str(t) ', Number of factors = ' mat2str(k) ', Time = ' mat2str(round(toc*1000)/1000) ', VA = ' mat2str(VA) '.'])
        plot(LB(1:t))
        drawnow
        tic
    end
end
[~,gamma] = q_lambda(mu,10000,family,B,D,tau2eta(tau,Transf),Transf);
%%% --Output
VBDAobj.gamma_mean = mean(gamma);
VBDAobj.gamma_sd = std(gamma);
VBDAobj.mu = mu;
VBDAobj.B = B;
VBDAobj.D = D;
VBDAobj.muz = muz;
VBDAobj.logsigmaz = logsigmaz;
VBDAobj.tau = tau;
VBDAobj.eta = tau2eta(tau,Transf);
VBDAobj.tauu = tauu;
VBDAobj.etau = tau2eta(tauu,Transf);
VBDAobj.C = C;
VBDAobj.lambda = lambda;
VBDAobj.LB = LB;
VBDAobj.mdist = mdist;
VBDAobj.y = y;
VBDAobj.family = family;
VBDAobj.p = p;
VBDAobj.k = k;
VBDAobj.S = S;
VBDAobj.nVB = nVB;
VBDAobj.VA = VA;
VBDAobj.Transf = Transf;

disp('----> VB estimation finished')
if nargin == 9 && ~isempty(filename)
    save(filename,'VBDAobj')
end
tau = kendalltau(gamma,family);
tau_mean = mean(tau);
tau_sd = std(tau);
numpar = n/p;
temp = zeros(p*2,numpar);
for i = 1:p
 temp((i-1)*2+1:i*2,:) = [round(VBDAobj.gamma_mean((i-1)*numpar+1:numpar*i)*1000)/1000 ; round(VBDAobj.gamma_sd((i-1)*numpar+1:numpar*i)*100)/100];
end
disp('--------------------------------------------------------------------------------------------------------------------')
disp(['----------------------------- Estimates of ' mat2str(numpar) ' parameters for each of the ' mat2str(p) ' pair-copulas -----------------------------'])
disp('--------------------------------------------------------------------------------------------------------------------')
temp1 = num2str([repmat(1234567891,p*2,1) temp]);
for i = 1:p
    temp1((i-1)*2+1,1:10) = ['c_' mat2str(i) '   Mean'];
    temp1(i*2,1:10) =      '      SD  ';
end
disp(temp1)
disp('--------------------------------------------------------------------------------------------------------------------')
tempp = [round(1000*tau_mean)/1000;round(100*tau_sd)/100];
tempp = tempp(:);
temp2 = num2str([repmat(1234567891,p*2,1) tempp]);
for i = 1:p
    temp2((i-1)*2+1,1:10) = ['c_' mat2str(i) '   Mean'];
    temp2(i*2,1:10) =      '      SD  ';
end
disp(['----------------------- Estimates of the Kendalls tau for each of the ' mat2str(p) ' pair-copulas -------------------------------'])
disp('--------------------------------------------------------------------------------------------------------------------')
disp(temp2)
disp('--------------------------------------------------------------------------------------------------------------------')
disp('--------------------------------------------------------------------------------------------------------------------')
clear grad_lambda_logq                                                     %removes persistent variables
