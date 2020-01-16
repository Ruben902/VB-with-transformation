function detadtau = deta_dtau(tau,Transf_type)
switch Transf_type
    case 'YJ'
        eta_max = 2;
        eta_min = 0;
        detadtau = (eta_max-eta_min)*exp(-tau)./((exp(-tau)+1).^2);
    case 'GH'
        detadtau(:,1) = ones(size(tau(:,1)));
        detadtau(:,2) = exp(tau(:,2));
    case 'iGH'
        eta_max = 1;
        eta_min = 1.0000e-6;
        detadtau(:,1) = ones(size(tau(:,1)));
        detadtau(:,2) = (eta_max-eta_min)*exp(-tau(:,2))./((exp(-tau(:,2))+1).^2);
end