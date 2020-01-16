function phi = tYJ(theta,eta)
if ~isequal(size(theta), size(eta)) && numel(eta)~=1
    error('The number of inputs should be equal to the number of tranformation parameters')
end
if numel(eta)==1
    eta = eta.*ones(size(theta));
end
thetaB = thetaBound();
theta_u = thetaB;
theta_l = -thetaB;
c1 =   theta<=theta_l;
c2 = theta_l<theta & theta<0;
c3 =       0<=theta & theta<theta_u;
c4 =  theta_u<=theta;
phi = zeros(size(theta));
%% c2
phic1 =-((-theta_l+1).^(2-eta(c1))-1)./(2-eta(c1))+(theta(c1)-theta_l);
%% c2
phic2 =-((-theta(c2)+1).^(2-eta(c2))-1)./(2-eta(c2));
%% c3
phic3 = ((theta(c3)+1).^eta(c3)-1)./eta(c3);
%% c4
phic4 = ((theta_u+1).^eta(c4)-1)./eta(c4)+(theta(c4)-theta_u);
eta0c3 = (eta==0 & c3);
eta2c2 = (eta==2 & c2);
eta0c4 = (eta==0 & c4);
eta2c1 = (eta==2 & c1);
phi(c1)=phic1;
phi(c2)=phic2;
phi(c3)=phic3;
phi(c4)=phic4;
phi(eta0c3) = log(theta(eta0c3)+1);
phi(eta2c2) = -log(-theta(eta2c2)+1);
phi(eta0c4) = log(theta(eta0c4)+1)+(theta(eta0c4)-theta_u);
phi(eta2c1) = -log(-theta(eta2c1)+1)+(theta(eta2c1)-theta_l);