%This function provides the derivatives in Table 4 of
%Smith, Loaiza-Maya and Nott (2019) "High-dimensional Copula Variational Approximation through
%Transformation".
function [dthetadB,dthetadDelta] = dtheta_dBDelta(eta,B,Delta,delta_phi,u,z,epsilon,gamma,Sigma_phiinv,Transf_type,phi,theta)
%Common constants
n = size(B,1);
K = size(B,2);
D = sparse(1:n,1:n,Delta,n,n);
Sigma_phi = B*B'+D.^2;
S_phi = sparse(1:n,1:n,diag(Sigma_phi),n,n);
S_phi_invhalf = sparse(1:n,1:n,diag(S_phi).^(-1/2),n,n);
delta_phi_tilde = (S_phi.^0.5)*delta_phi;
persistent  Kcomnp
if  isempty(Kcomnp)
    Kcomnp = communtationM(ones(n,K));
end
xi = (B*z(:)+Delta.*epsilon(:));

In = sparse(1:n,1:n,1,n,n);
P = sparse(1:(n+1):n^2,1:n,1);

delta_Sigmainv = delta_phi_tilde'*Sigma_phiinv;
xi_Sigmainv = xi'*Sigma_phiinv;
Sigmainv_B = Sigma_phiinv*B;

const2 = 1-(delta_Sigmainv)*delta_phi_tilde;

M1 = delta_phi_tilde*delta_phi_tilde';
M2 =  xi_Sigmainv*delta_phi_tilde*In+kron(xi_Sigmainv,delta_phi_tilde);
M3 = 0.5*gamma*(const2.^(-1/2))*delta_phi_tilde;
M4 = gamma*sqrt(const2);
M5 = delta_phi_tilde'*kron(Sigmainv_B,delta_Sigmainv);
M6 = delta_phi_tilde'*kron(Sigma_phiinv,delta_Sigmainv*B)*Kcomnp;
M7 = sparse(repmat(1:n,1,K),1:(n*K),B(:),n,n*K);
M8 = sparse(1:n,1:n,delta_phi,n,n)*S_phi_invhalf*M7;
M9 = kron(-2*M3,delta_Sigmainv*M8);
M10 = kron(M3,M5)+kron(M3,M6);
M11 = M2*M8;
M12 = -kron(xi_Sigmainv*B,Sigma_phiinv)-kron(xi_Sigmainv,Sigmainv_B)*Kcomnp;  
TB0 = abs(u)*M8;
TB1 = sparse(kron(z(:)',In));
M13 = M12+Sigma_phiinv*TB1;
TB2 = M11+M1*M13;
TB3 = M9+M10+M4*M8;

[M11D,M5D] = const8deri3simple(Delta,B,z,epsilon,delta_phi_tilde);
M6D = sparse(repmat(1:n,1,n),1:(n*n),D(:),n,n*n);
M7D = (sparse(1:n,1:n,delta_phi,n,n)*S_phi_invhalf*M6D)*P;
M8D = kron(-2*M3,delta_Sigmainv*M7D);
M9D = -kron(M3,delta_phi_tilde'*M5D);
M10D = M2*M7D;
TD0 = abs(u)*M7D;
TD1 = sparse(kron(epsilon(:)',In))*P;
M12D = M11D+Sigma_phiinv*TD1;
TD2 = M10D+M1*M12D;
TD3 = M8D+M9D+M4*M7D;

switch Transf_type
    case 'YJ'
        dthetdphi = dtheta_dphi(phi,eta);
    case 'GH'
        dthetdphi = digh(phi,eta(:,1),eta(:,2),theta);
    case 'iGH'
        dthetdphi = dgh(phi,eta(:,1),eta(:,2));
end

dthetadphi = sparse(1:n,1:n,dthetdphi,n,n);      %%Equivalent to dt(phi)/dphi
dthetadB = dthetadphi*(TB0+TB1-TB2+TB3);
dthetadDelta = dthetadphi*(TD0+TD1-TD2+TD3);
