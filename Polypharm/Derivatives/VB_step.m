function [LowerB,B,mu,d,eta,alphaphi,ADA,Siginv] = VB_step(B, mu, d, eta, alphaphi, logpost,dlogpost,Siginv,ADA,p,SkewNorm,Transf)
if strcmp(Transf,'None') || strcmp(Transf,'none') || strcmp(Transf,'') || isempty(Transf)
    Transformation = 0;
    Transf = 'YJ';
else
    Transformation = 1;
end

rho = ADA.rho;
eps_step = ADA.eps_step;
oldEdelta2_mu = ADA.Edelta2_mu;
oldEg2_mu = ADA.Eg2_mu;

oldEdelta2_B = ADA.Edelta2_B;
oldEg2_B = ADA.Eg2_B;

oldEdelta2_d = ADA.Edelta2_d;
oldEg2_d = ADA.Eg2_d;

oldEdelta2_alphaphi = ADA.Edelta2_alphaphi;
oldEg2_alphaphi = ADA.Eg2_alphaphi;

oldEdelta2_tau = ADA.Edelta2_tau;
oldEg2_tau = ADA.Eg2_tau;

q = length(mu); %%% asssume each b_i is 1

zeps = randn(p+q,1)';
ugamma = randn(2,1)';
u =  ugamma(1);
gamma =  ugamma(2);
z = zeps(1:p)';
eps = zeps((p+1):end)';

Sigma = B*B'+diag(d.^2);
S_phi = sparse(1:q,1:q,diag(Sigma),q,q);
S_phi_invhalf = sparse(1:q,1:q,diag(S_phi).^(-1/2),q,q);
Omega = S_phi_invhalf*(Sigma)*S_phi_invhalf;
dphi = alphaphi2deltaphi(alphaphi,Omega);
dphitilde = (S_phi.^0.5)*dphi;
tau = eta2tau(eta,Transf);

phi = mu + dphitilde*abs(u)+(speye(q)-dphitilde*dphitilde'*Siginv)*(B*z + d.*eps)+gamma*sqrt(1-dphitilde'*Siginv*dphitilde)*dphitilde;
switch Transf
    case 'YJ'
        theta = tYJi(phi,eta);
    case 'GH'
        theta = igh(phi,eta(:,1),eta(:,2));
    case 'iGH'
        theta = gh(phi,eta(:,1),eta(:,2));
end


[L_mu,L_B,L_d,L_alphaphi,L_tau,g] = gradient_compute(theta,mu,B,z,d,eps,dphi,tau,u,gamma,logpost,dlogpost,Siginv,Transf,phi);


L_B(~tril(ones(size(L_B))))  = 0;

%% mu update

ADA.Eg2_mu = rho*oldEg2_mu + (1-rho)*L_mu.^2;
Change_delta_mu = sqrt(oldEdelta2_mu + eps_step)./sqrt(ADA.Eg2_mu + eps_step).*L_mu;


mu = mu + Change_delta_mu;
ADA.Edelta2_mu = rho*oldEdelta2_mu + (1- rho)*Change_delta_mu.^2;


%% B update

vecL_B = L_B(:);

ADA.Eg2_B = rho*oldEg2_B + (1-rho)*vecL_B.^2;
Change_delta_B = sqrt(oldEdelta2_B + eps_step)./sqrt(ADA.Eg2_B + eps_step).*vecL_B;


B = B + vec2mat(Change_delta_B,q,p);
ADA.Edelta2_B = rho*oldEdelta2_B + (1- rho)*Change_delta_B.^2;


%% d update

ADA.Eg2_d = rho*oldEg2_d + (1-rho)*L_d.^2;
Change_delta_d = sqrt(oldEdelta2_d + eps_step)./sqrt(ADA.Eg2_d + eps_step).*L_d;

d = d + Change_delta_d;
ADA.Edelta2_d = rho*oldEdelta2_d + (1- rho)*Change_delta_d.^2;
%% alpha_phi update

ADA.Eg2_alphaphi = rho*oldEg2_alphaphi + (1-rho)*L_alphaphi.^2;
Change_delta_alphaphi = sqrt(oldEdelta2_alphaphi + eps_step)./sqrt(ADA.Eg2_alphaphi + eps_step).*L_alphaphi;

alphaphi = alphaphi + Change_delta_alphaphi*SkewNorm;
ADA.Edelta2_alphaphi = rho*oldEdelta2_alphaphi + (1- rho)*Change_delta_alphaphi.^2;
%% tau update
ADA.Eg2_tau = rho*oldEg2_tau + (1-rho)*L_tau.^2;
Change_delta_tau = sqrt(oldEdelta2_tau + eps_step)./sqrt(ADA.Eg2_tau + eps_step).*L_tau;
taustep = reshape(Change_delta_tau,size(tau));
tau = tau + taustep*Transformation;
ADA.Edelta2_tau = rho*oldEdelta2_tau + (1- rho)*Change_delta_tau.^2;
eta = tau2eta(tau,Transf);
%% Lowerbound
loghtheta = g;
phi = mu + dphitilde*abs(u)+(eye(q)-dphitilde*dphitilde'*Siginv)*(B*z + d.*eps)+gamma*sqrt(1-dphitilde'*Siginv*dphitilde)*dphitilde;
switch Transf
    case 'YJ'
        theta = tYJi(phi,eta);
        dt_dtheta = dtYJ_dtheta(theta,eta);
    case 'GH'
        theta = igh(phi,eta(:,1),eta(:,2));
        dt_dtheta = dgh(theta,eta(:,1),eta(:,2));
    case 'iGH'
        theta = gh(phi,eta(:,1),eta(:,2));
        dt_dtheta = digh(theta,eta(:,1),eta(:,2),phi);
end
phimmu = phi-mu;
Dinv2B = bsxfun(@times,B,1./d.^2);
Siginv = diag(1./d.^2) - Dinv2B/(eye(p) + B'*Dinv2B)*Dinv2B';
Sigma = B*B'+diag(d.^2);
S_phi = sparse(1:q,1:q,diag(Sigma),q,q);
S_phi_invhalf = sparse(1:q,1:q,diag(S_phi).^(-1/2),q,q);
%Omega = S_phi_invhalf*(B*B'+diag(d.^2))*S_phi_invhalf;
PhiNorm = normcdf(alphaphi'*S_phi_invhalf*phimmu);
logphiNorm = logmvnpdf(phi',mu',Sigma);
logJacobian = sum(log(dt_dtheta));
LowerB = loghtheta - logJacobian - log(2) - logphiNorm - log(PhiNorm);



