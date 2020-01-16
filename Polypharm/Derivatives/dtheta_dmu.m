function dthetadmu=dtheta_dmu(phi,eta,Transf_type,theta)
if nargin ==2
    Transf_type ='YJ';
end

if strcmp(Transf_type,'GH')
    if size(eta,2)~=2
        error('GH transformation requires 2 parameters. Make sure eta has two columns')
    end
end

if strcmp(Transf_type,'iGH')
    if size(eta,2)~=2
        error('iGH transformation requires 2 parameters. Make sure eta has two columns')
    end
end

switch Transf_type
    case 'YJ'
        dthetadphi = dtheta_dphi(phi,eta);
    case 'GH'
        dthetadphi = digh(phi,eta(:,1),eta(:,2),theta);
    case 'iGH'
        dthetadphi = dgh(phi,eta(:,1),eta(:,2));
end

dthetadmu=sparse(1:size(eta,1),1:size(eta,1),dthetadphi,size(eta,1),size(eta,1));