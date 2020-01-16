function  dthetadalphaphi = dtheta_dalphaphi(eta,B,Delta,delta_phi,u,z,epsilon,gamma,Sigma_phiinv,Transf_type,phi,theta)
if isempty(Transf_type)
   Transf_type ='YJ'; 
end
n = size(epsilon,1);
S_phi = sparse(1:n,1:n,diag(B*B'+diag(Delta.^2)),n,n);
S_phi_invhalf = sparse(1:n,1:n,diag(S_phi).^(-1/2),n,n);
Omega = S_phi_invhalf*(B*B'+diag(Delta.^2))*S_phi_invhalf;
delta_phi_tilde = (S_phi.^0.5)*delta_phi;
Term1 = Omega\delta_phi;
alpha_phi = Term1/sqrt(1-delta_phi'*Term1);

xi = (B*z+Delta.*epsilon);
In = sparse(1:n,1:n,ones(n,1),n,n);

switch Transf_type
    case 'YJ'
        dthetadphi=dtheta_dphi(phi,eta);      %%Equivalent to dt(phi)/dphi
    case 'GH'
        dthetadphi = digh(phi,eta(:,1),eta(:,2),theta);
    case 'iGH'
        dthetadphi = dgh(phi,eta(:,1),eta(:,2));
end

dthetadphi = sparse(1:n,1:n,dthetadphi,n,n);
Omegaalpha = Omega*alpha_phi;


M1 = (1-delta_phi_tilde'*Sigma_phiinv*delta_phi_tilde)^0.5;
M2 = -(1/M1)*delta_phi_tilde'*Sigma_phiinv;
M3 = kron(gamma*In(:),M2);
M4 = kron(delta_phi_tilde',In);
M5 = M4*M3+gamma*M1*In;
M6 = xi'*Sigma_phiinv*delta_phi_tilde*In+kron(xi'*Sigma_phiinv,delta_phi_tilde);
M7 = abs(u)*In-M6+M5;
M8 = dthetadphi*M7*(S_phi.^0.5);
M9 = (1+alpha_phi'*Omegaalpha);
M10 = M9^(-3/2)*(M9*Omega-Omegaalpha*Omegaalpha');
dthetadalphaphi = M8*M10;

