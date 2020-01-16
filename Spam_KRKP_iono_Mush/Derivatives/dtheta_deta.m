function dthetadeta=dtheta_deta(phi,eta)
if ~isequal(size(phi), size(eta)) && numel(eta)~=1
    error('The number of inputs should be equal to the number of tranformation parameters')
end
if numel(eta)==1
    eta = eta.*ones(size(phi));
end
eta(abs(eta)<0.0000001) = 0.0000001;

c2 = phi<0;
c3 = 0<=phi;
dthetadeta = zeros(size(phi));
%% c2
dthetadeta(c2) =-((1-phi(c2).*(2-eta(c2))).^(1./(2-eta(c2))).*(phi(c2)./((2-eta(c2)).*(1-phi(c2).*(2-eta(c2))))+((log(1-phi(c2).*(2-eta(c2))))./((2-eta(c2)).^2))));
%% c3
dthetadeta(c3) = (1+phi(c3).*eta(c3)).^(1./eta(c3)).*((phi(c3)./(eta(c3).*(1+phi(c3).*eta(c3))))-((log(1+phi(c3).*eta(c3)))./(eta(c3).^2)));
