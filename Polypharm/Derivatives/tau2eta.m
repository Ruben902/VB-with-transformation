% tau is the Fisher transformation of eta
function eta = tau2eta(tau,Transf_type)
if strcmp(Transf_type,'GH')
    if size(tau,2)~=2
        error('GH transformation requires 2 parameters. Make sure eta has two columns')
    end
end

if strcmp(Transf_type,'iGH')
    if size(tau,2)~=2
        error('iGH transformation requires 2 parameters. Make sure eta has two columns')
    end
end

switch Transf_type
    case 'YJ'
        eta_max = 2;
        eta_min = 0;
        %eta = 2./(exp(-tau)+1);
        eta = (eta_max-eta_min)./(exp(-tau)+1)+eta_min;
    case 'GH'
        eta(:,1) = tau(:,1);
        eta(:,2) = exp(tau(:,2))+1.0000e-6;
    case 'iGH'
        eta(:,1) = tau(:,1);
        eta_max = 1;
        eta_min = 1.0000e-6;
        eta(:,2) = (eta_max-eta_min)./(exp(-tau(:,2))+1)+eta_min;
end


