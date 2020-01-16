function phi = tYJ(theta,eta)
if ~isequal(size(theta), size(eta)) && numel(eta)~=1
    error('The number of inputs should be equal to the number of tranformation parameters')
end
if numel(eta)==1
    eta = eta.*ones(size(theta));
end

c2 = theta<0;
c3 = 0<=theta;
phi = zeros(size(theta));
%% c2
phic2 =-((-theta(c2)+1).^(2-eta(c2))-1)./(2-eta(c2));
%% c3
phic3 = ((theta(c3)+1).^eta(c3)-1)./eta(c3);

eta0c3 = (eta==0 & c3);
eta2c2 = (eta==2 & c2);

phi(c2)=phic2;
phi(c3)=phic3;

phi(eta0c3) = log(theta(eta0c3)+1);
phi(eta2c2) = -log(-theta(eta2c2)+1);
