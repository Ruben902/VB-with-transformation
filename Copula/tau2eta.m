% tau is the Fisher transformation of eta
function eta = tau2eta(tau,Transf_type)
if strcmp(Transf_type,'GH')
    if size(tau,2)~=2
       error('GH transformation requires 2 parameters. Make sure eta has two columns') 
    end
end

switch Transf_type
    case 'YJ'
        eta_max = 1.5;
        eta_min = 0.5;
        %eta = 2./(exp(-tau)+1);
        eta = (eta_max-eta_min)./(exp(-tau)+1)+eta_min;
    case 'GH'
        eta(:,1) = tau(:,1);
        eta(:,2) = exp(tau(:,2))+1.0000e-6;
    case 'none'
        eta=[];
end


