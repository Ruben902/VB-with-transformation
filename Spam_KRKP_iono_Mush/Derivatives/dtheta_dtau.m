function dthetadtau = dtheta_dtau(phi,tau,theta,Transf_type)
if nargin ==3
   Transf_type ='YJ'; 
end

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
        q = size(tau,1);
        eta = tau2eta(tau,'YJ');
        dthetadeta = dtheta_deta(phi,eta);
        detadtau = deta_dtau(tau,'YJ');
        dthetadtau = sparse(1:q,1:q,dthetadeta.*detadtau,q,q);
    case 'GH'
        q = size(tau,1);
        eta = tau2eta(tau,'GH');
        dthetadphi = digh(phi,eta(:,1),eta(:,2),theta);
        dghg = dgh_g(theta,eta(:,1),eta(:,2));
        dghh = dgh_h(theta,eta(:,1),eta(:,2));
        dphideta = [dghg dghh];
        detadtau = deta_dtau(tau,'GH');
        dthetadtautemp = -repmat(dthetadphi,1,2).*dphideta.*detadtau;
        dthetadtau = sparse(1:q,1:q,dthetadtautemp(:,1),q,q*2)+sparse(1:q,(q+1):q*2,dthetadtautemp(:,2),q,q*2);
        %The negative value comes from the Triple Product Rule, or cyclic
        %chain rule
    case 'iGH'
        q = size(tau,1);
        eta = tau2eta(tau,'iGH');
        dthetadg = dgh_g(phi,eta(:,1),eta(:,2));
        dthetadh = dgh_h(phi,eta(:,1),eta(:,2));
        dthetadeta = [dthetadg dthetadh];
        detadtau = deta_dtau(tau,'iGH');
        dthetadtautemp = dthetadeta.*detadtau;
        dthetadtau = sparse(1:q,1:q,dthetadtautemp(:,1),q,q*2)+sparse(1:q,(q+1):q*2,dthetadtautemp(:,2),q,q*2);
end

