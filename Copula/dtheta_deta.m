function dthetadeta=dtheta_deta(phi,eta)
if ~isequal(size(phi), size(eta)) && numel(eta)~=1
    error('The number of inputs should be equal to the number of tranformation parameters')
end
if numel(eta)==1
    eta = eta.*ones(size(phi));
end
eta(abs(eta)<0.0000001) = 0.0000001;
thetaB = thetaBound();
theta_u = thetaB;
theta_l = -thetaB;
phi_u = tYJ(zeros(size(eta))+theta_u,eta);
phi_l = tYJ(zeros(size(eta))+theta_l,eta);
c1 =   phi<=phi_l;
c2 = phi_l< phi & phi<0;
c3 =     0<=phi & phi<phi_u;
c4 = phi_u<=phi;
dthetadeta = zeros(size(phi));
%% c1
dthetadeta(c1) = -((2-eta(c1)).*((1-theta_l).^(2-eta(c1)))*(log(-theta_l+1))-(1-theta_l).^(2-eta(c1))+1)./((2-eta(c1)).^2);
%% c2
dthetadeta(c2) =-((1-phi(c2).*(2-eta(c2))).^(1./(2-eta(c2))).*(phi(c2)./((2-eta(c2)).*(1-phi(c2).*(2-eta(c2))))+((log(1-phi(c2).*(2-eta(c2))))./((2-eta(c2)).^2))));
%% c3
dthetadeta(c3) = (1+phi(c3).*eta(c3)).^(1./eta(c3)).*((phi(c3)./(eta(c3).*(1+phi(c3).*eta(c3))))-((log(1+phi(c3).*eta(c3)))./(eta(c3).^2)));
%% c4
dthetadeta(c4) = -((eta(c4).*(1+theta_u).^eta(c4)*log(theta_u+1)-((1+theta_u).^(eta(c4))-1))./(eta(c4).^2));


