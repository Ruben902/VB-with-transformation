%Inverse of the Yeo Jonhston transformation
function theta = tYJi(phi,eta)
if ~isequal(size(phi), size(eta)) && numel(eta)~=1
    error('The number of inputs should be equal to the number of tranformation parameters')
end
if numel(eta)==1
    eta = eta.*ones(size(phi));
end

c2 = phi<0;
c3 = 0<=phi;

theta = zeros(size(phi));
%% c2
thetac2 =1-(1-phi(c2).*(2-eta(c2))).^(1./(2-eta(c2)));
%% c3
thetac3 = (phi(c3).*eta(c3)+1).^(1./eta(c3))-1;
%%
eta0c3 = (eta==0 & c3);
eta2c2 = (eta==2 & c2);
theta(c2)=thetac2;
theta(c3)=thetac3;
theta(eta0c3) = exp(phi(eta0c3))-1;
theta(eta2c2) = 1-exp(-phi(eta2c2));