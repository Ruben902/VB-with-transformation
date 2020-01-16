% tau is the Fisher transformation of eta
function tau = eta2tau(eta,Transf_type)
if strcmp(Transf_type,'GH')
    if size(eta,2)~=2
       error('GH transformation requires 2 parameters. Make sure eta has two columns') 
    end
end

switch Transf_type
    case 'YJ'
        eta_max = 1.5;
        eta_min = 0.5;
        %tau = log(eta./(2-eta));
        tau=-log(((eta_max-eta_min)./(eta-eta_min))-1);
    case 'GH'
        tau(:,1) = eta(:,1);
        tau(:,2) = log(eta(:,2)-1.0000e-6);
end

