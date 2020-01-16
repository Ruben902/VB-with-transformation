function [theta,T_theta] = q_lambda(mu,S,family,B,D,eta,Transf)
npar = length(mu);
if strcmp(family,'gumbel_mix')
    UB = 0.99*ones(S,npar);
    LB = zeros(S,npar);
end
if strcmp(family,'t')
    lag = npar/2;
    UB = repmat(repmat([0.99 40],1,lag),S,1);
    LB = repmat(repmat([-0.99 0.03],1,lag),S,1);
end
k = size(B,2);
z = randn(k,S);
eps = randn(npar,S);
switch Transf
    case 'YJ'
        phi = repmat(mu,1,S)+B*z+D*eps;
        Eta = repmat(eta,1,S);
        theta = tYJi(phi,Eta);
        theta = theta';
    case 'GH'
        phi = repmat(mu,1,S)+B*z+D*eps;
        Eta_g = repmat(eta(:,1),1,S);
        Eta_h = repmat(eta(:,2),1,S);
        theta = igh(phi,Eta_g,Eta_h);
        theta = theta';
    case 'none'
        theta = repmat(mu,1,S)+B*z+D*eps;
        theta = theta';
end
T_theta = (UB-LB)./(exp(-theta)+1)+LB;

