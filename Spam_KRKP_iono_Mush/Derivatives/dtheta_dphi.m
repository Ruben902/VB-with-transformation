%%Equivalent to dt(phi)/dphi
function dthetadphi=dtheta_dphi(phi,eta)
if ~isequal(size(phi), size(eta)) && numel(eta)~=1
    error('The number of inputs should be equal to the number of tranformation parameters')
end
if numel(eta)==1
    eta = eta.*ones(size(phi));
end
eta(abs(eta)<0.0000001) = 0.0000001;

c2 = phi<0;
c3 = 0<=phi ;
dthetadphi = zeros(size(phi));
dthetadphi(c2) = (1-phi(c2).*(2-eta(c2))).^((eta(c2)-1)./(2-eta(c2)));
dthetadphi(c3) = (1+phi(c3).*eta(c3)).^((1-eta(c3))./(eta(c3)));
