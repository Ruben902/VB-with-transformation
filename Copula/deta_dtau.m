function detadtau = deta_dtau(tau,Transf_type)
switch Transf_type
    case 'YJ'
        eta_max = 1.5;
        eta_min = 0.5;
        detadtau = (eta_max-eta_min)*exp(-tau)./((exp(-tau)+1).^2);
    case 'GH'
        detadtau(:,1) = ones(size(tau(:,1)));
        detadtau(:,2) = exp(tau(:,2));
end