%Inverse of the Yeo Jonhston transformation
function theta = tYJi(phi,eta)
if ~isequal(size(phi), size(eta)) && numel(eta)~=1
    error('The number of inputs should be equal to the number of tranformation parameters')
end
if numel(eta)==1
    eta = eta.*ones(size(phi));
end
thetaB = thetaBound();
theta_u = thetaB;
theta_l = -thetaB;
phi_u = tYJ(zeros(size(eta))+theta_u,eta);
phi_l = tYJ(zeros(size(eta))+theta_l,eta);
c1 =   phi<=phi_l;
c2 = phi_l< phi & phi<0;
c3 =     0<=phi & phi<phi_u;
c4 = phi_u<=phi;
theta = zeros(size(phi));
%% c1
thetac1 = phi(c1)-phi_l(c1)+theta_l;
%% c2
thetac2 =1-(1-phi(c2).*(2-eta(c2))).^(1./(2-eta(c2)));
%% c3
thetac3 = (phi(c3).*eta(c3)+1).^(1./eta(c3))-1;
%% c4
thetac4 = phi(c4)-phi_u(c4)+theta_u;
%%
eta0c3 = (eta==0 & c3);
eta2c2 = (eta==2 & c2);
theta(c1)=thetac1;
theta(c2)=thetac2;
theta(c3)=thetac3;
theta(c4)=thetac4;
theta(eta0c3) = exp(phi(eta0c3))-1;
theta(eta2c2) = 1-exp(-phi(eta2c2));