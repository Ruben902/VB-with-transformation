% phi = tYJ(theta,eta)
% theta = tYJ(phi,eta)
% This function computes dt(theta)/dtheta 
%Jacobian of the transformation
function ddphi_dthetadeta = ddtYJ_dthetadeta(theta,eta)
if ~isequal(size(theta), size(eta)) && numel(eta)~=1
    error('The number of inputs should be equal to the number of tranformation parameters')
end
if numel(eta)==1
    eta = eta.*ones(size(theta));
end
eta(abs(eta)<0.0000001) = 0.0000001;

c2 = theta<0;
c3 = 0<=theta;
ddphi_dthetadeta = zeros(size(theta));
%% c2
ddphi_dthetadeta(c2) = -(-theta(c2)+1).^(1-eta(c2)).*log(abs(-theta(c2)+1));
%% c3
ddphi_dthetadeta(c3) = (theta(c3)+1).^(eta(c3)-1).*log(theta(c3)+1);
