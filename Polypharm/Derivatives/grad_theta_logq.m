%Refer to notes on derivatives
function T = grad_theta_logq(theta,eta,mu_phi,B,Delta,delta_phi,Transf_type,phi)
if strcmp(Transf_type,'GH')
    if size(eta,2)~=2
        error('GH transformation requires 2 parameters. Make sure eta has two columns')
    end
end

if strcmp(Transf_type,'iGH')
    if size(eta,2)~=2
        error('iGH transformation requires 2 parameters. Make sure eta has two columns')
    end
end

n = length(theta);
K = size(B,2);
D = sparse(1:n,1:n,Delta,n,n);
Tq1vec = sparse(n,1);
Sigma_phi = B*B'+D.^2;
Dm2 = sparse(1:n,1:n,1./(Delta.^2),n,n);
Sigma_phiinv = Dm2-Dm2*B/(sparse(1:K,1:K,1,K,K)+B'*Dm2*B)*B'*Dm2;  
S_phi = diag(diag(Sigma_phi));
S_phi_invhalf = diag(diag(S_phi).^(-1/2));
Omega_phi = S_phi_invhalf*Sigma_phi*S_phi_invhalf;%
switch Transf_type
    case 'YJ'
        dt_dtheta = dtYJ_dtheta(theta,eta);
        %Tq1
        c2 = theta<0;
        c3 = 0<=theta;
        Tq1vec(c2) = (eta(c2)-1)./(-theta(c2)+1);
        Tq1vec(c3) = (eta(c3)-1)./(theta(c3)+1);
    case 'GH'
        dt_dtheta = dgh(theta,eta(:,1),eta(:,2));
        ddt_dtheta = ddgh(theta,eta(:,1),eta(:,2));
        Tq1vec = ddt_dtheta./dt_dtheta;
    case 'iGH'
        dt_dtheta = digh(theta,eta(:,1),eta(:,2),phi);
        ddt_dtheta = ddigh(theta,eta(:,1),eta(:,2),phi);
        Tq1vec = ddt_dtheta./dt_dtheta;
end

Temp1 = Omega_phi\delta_phi;
alpha_phi = Temp1/sqrt(1-delta_phi'*Temp1);
arg1 = alpha_phi'*S_phi_invhalf*(phi-mu_phi);
phiPhi = normpdf(arg1)./normcdf(arg1);
dtdtheta = sparse(1:n,1:n,dt_dtheta,n,n);
Tq2vec = -dtdtheta'*Sigma_phiinv*(phi-mu_phi);
Tq3vec = dtdtheta*S_phi_invhalf*alpha_phi*phiPhi;

T = Tq1vec+Tq2vec+Tq3vec;
