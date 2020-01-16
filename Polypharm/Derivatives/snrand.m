function theta = snrand(S,eta,mu_phi,B,Delta,alpha_phi,Transf_type)
if nargin ==6
   Transf_type ='YJ'; 
end
if strcmp(Transf_type,'')
   Transf_type ='YJ'; 
end
if strcmp(Transf_type,'GH')
    if size(eta,2)~=2
       error('GH transformation requires 2 parameters. Make sure eta has two columns') 
    end
end
n = size(B,1);
K = size(B,2);
z = randn(K,S);
u = randn(1,S);
gamma = randn(1,S);
epsilon = randn(n,S);
In = sparse(1:n,1:n,1,n,n);
D = sparse(1:n,1:n,Delta,n,n);
Sigma_phi = B*B'+D.^2;
S_phi = sparse(1:n,1:n,diag(Sigma_phi),n,n);
S_phiinvhalf = sparse(1:n,1:n,diag(Sigma_phi).^(-0.5),n,n);
Dm2 = sparse(1:n,1:n,1./(Delta.^2),n,n);
Sigma_phiinv = Dm2-Dm2*B/(sparse(1:K,1:K,1,K,K)+B'*Dm2*B)*B'*Dm2;                                     %Woodbury formula
Omega = S_phiinvhalf*Sigma_phi*S_phiinvhalf;
delta_phi = alphaphi2deltaphi(alpha_phi,Omega);
delta_phi_tilde = (S_phi.^0.5)*delta_phi;
Term1 = repmat(mu_phi,1,S)+repmat(delta_phi_tilde,1,S).*abs(repmat(u,n,1));
Term2 = (In - delta_phi_tilde*delta_phi_tilde'*Sigma_phiinv)*(B*z+repmat(Delta,1,S).*epsilon);
Term3 = repmat(gamma,n,1).*repmat(delta_phi_tilde*(1-delta_phi_tilde'*Sigma_phiinv*delta_phi_tilde)^(0.5),1,S);

switch Transf_type
    case 'YJ'
         theta = tYJi(Term1+Term2+Term3,repmat(eta,1,S));
    case 'GH'
        Theta = igh(Term1(:)+Term2(:)+Term3(:),repmat(eta(:,1),S,1),repmat(eta(:,2),S,1));
        theta = reshape(Theta,size(eta,1),S);
    case 'iGH'
        Theta = gh(Term1(:)+Term2(:)+Term3(:),repmat(eta(:,1),S,1),repmat(eta(:,2),S,1));
        theta = reshape(Theta,size(eta,1),S);
end
