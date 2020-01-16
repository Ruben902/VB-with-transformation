function [L_mu,L_B,L_d,L_alphaphi,L_tau,g] = gradient_compute(theta,mu,B,z,d,eps,dphi,tau,u,gamma,logpost,dlogpost,Siginv,Transf,phi)
q = size(B,1);
p = size(B,2);
eta = tau2eta(tau,Transf);

g = logpost(theta);
delta_logh = dlogpost(theta);

delta_logq = grad_theta_logq(theta,eta,mu,B,d,dphi,Transf,phi);
L_mu = dtheta_dmu(phi,eta,Transf,theta)*(delta_logh - delta_logq);
[dtheta_dB,dtheta_dd] = dtheta_dBDelta(eta,B,d,dphi,u,z,eps,gamma,Siginv,Transf,phi,theta);
L_B = reshape(dtheta_dB'*(delta_logh - delta_logq),q,p);
L_d = dtheta_dd'*(delta_logh - delta_logq);
L_alphaphi = dtheta_dalphaphi(eta,B,d,dphi,u,z,eps,gamma,Siginv,Transf,phi,theta)'*(delta_logh - delta_logq);
L_tau =  (dtheta_dtau(phi,tau,theta,Transf)')*(delta_logh - delta_logq);    %% Derivatives with respect to tau, the fisher transformation of eta







