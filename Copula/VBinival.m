function [ai,bi,n,LB,ct,lambda,t,mu,Sigma,muz,logsigmaz,C,tau,tauu,T,B,D,Transf] = VBinival(y,mdist,family,p,k,S,nVB,VA)
[ai,bi] = ai_bi_compute(y,mdist);
LB = zeros(nVB,1);
T = size(ai,1);
if strcmp(family,'gumbel_mix')
    Inipref = [-3 0 -3 0 3];
elseif strcmp(family,'t_mix')

end
if VA == 1
    muz = [];
    logsigmaz = [];
    C = [];
elseif VA == 2
    muz = zeros(T,1);
    logsigmaz = zeros(T,1);
    C = [];
elseif VA == 3
    muz = zeros(T,1);
    logsigmaz = [];
    C = sparse(1:T,1:T,1,T,T);
end
n = length(Inipref)*p;
mu = repmat(Inipref(:),p,1);
Sigma = eye(n)*0.1; 
Transf = 'YJ';
q = size(mu,1);
switch Transf
    case 'YJ'
        eta = ones(q,1);
        etau = ones(T,1);
        if VA==1
            etau=[];
        end
    case 'GH'
        eta = zeros(q,2);
        eta(:,1)=0.001;
        eta(:,2)=0.1;
        etau = zeros(T,2);
        etau(:,1)=0.001;
        etau(:,2)=0.1;
    case 'none'
        eta =[];
        etau = [];
end
[lambda,B,D,tau,tauu] = musigma2lambda(mu,Sigma,muz,logsigmaz,C,k,eta,etau,Transf);
[theta,Gamma] = q_lambda(mu,S,family,B,D,eta,Transf);
[u,z,vareps] = q_u(ai,bi,muz,logsigmaz,C,S,VA,etau,Transf);
[~,LogLhat,log_theta,log_q_lambda,grad_lambda] = gradient_compute(theta,Gamma,u,z,vareps,ai,bi,family,lambda,k,VA,[],Transf);
ct = ct_constant(LogLhat,log_theta,log_q_lambda,grad_lambda);
t = 0;