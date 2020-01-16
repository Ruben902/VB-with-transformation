% phi = tYJ(theta,eta)
% theta = tYJ(phi,eta)
% This function computes dt(theta)/dtheta 
%Jacobian of the transformation
function dphi_dtheta = dtYJ_dtheta(theta,eta)
if ~isequal(size(theta), size(eta)) && numel(eta)~=1
    error('The number of inputs should be equal to the number of tranformation parameters')
end
if numel(eta)==1
    eta = eta.*ones(size(theta));
end
eta(abs(eta)<0.0000001) = 0.0000001;

c2 = theta<0;
c3 = 0<=theta ;
dphi_dtheta = zeros(size(theta));
%% c2
dphi_dtheta(c2) = (-theta(c2)+1).^(1-eta(c2));
%% c3
dphi_dtheta(c3) = (theta(c3)+1).^(eta(c3)-1);
