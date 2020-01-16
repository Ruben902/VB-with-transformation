function grad_lambda = grad_lambda_logq(theta,z,vareps,lambda,k,T,VA,Transf)
persistent indPositiveC id1 id2 IT2 IT Kcom;    % Be aware!, persistent variables stay saved in the function. To reset use the MATLAB command: clear grad_lambda_logq
S = size(theta,1);
d = size(theta,2);
[mu,Sigma,B,D,muz,logsigmaz,C,tau,tauu] = lambda2musigma(lambda,d,T,k,VA,Transf);
switch Transf
    case 'YJ'
        eta = tau2eta(tau,Transf);
        phi = tYJ(theta,repmat(eta',S,1));
        logJacob = log(dtYJ_dtheta(theta,repmat(eta',S,1)));
        dtdtheta = exp(logJacob);
        dthetadeta = dtheta_deta(phi,repmat(eta',S,1));
        ddt_dthetadeta = ddtYJ_dthetadeta(theta,repmat(eta',S,1));
        detadtau = deta_dtau(tau,Transf);
        ddt_dthetadtau = ddt_dthetadeta.*repmat(detadtau',S,1);
        if VA>1
            etau = tau2eta(tauu,Transf);
            logJacob_u = log(dtYJ_dtheta(vareps,repmat(etau,1,S)));
            dtdvareps = exp(logJacob_u);
            ddt_dvarepsdeta = ddtYJ_dthetadeta(vareps,repmat(etau,1,S));
            detaudtauu = deta_dtau(tauu,Transf);
            ddt_dvarepsdtau = ddt_dvarepsdeta.*repmat(detaudtauu,1,S);
        end
    case 'GH'
        eta = tau2eta(tau,Transf);
        phi = gh(theta,repmat(eta(:,1)',S,1),repmat(eta(:,2)',S,1));
        logJacob = log(dgh(theta,repmat(eta(:,1)',S,1),repmat(eta(:,2)',S,1)));
        dtdtheta = exp(logJacob);
        
        etau = tau2eta(tauu,Transf);
        logJacob_u = log(dgh(vareps,repmat(etau(:,1),1,S),repmat(etau(:,2),1,S)));
        dtdvareps = exp(logJacob_u);
    case 'none'
        phi = theta;
        logJacob=0;
        dtdtheta = exp(logJacob);
        detadtau = zeros(d,1);
        dthetadeta = zeros(S,d);
        ddt_dthetadtau= zeros(S,d);
end
Dm2 = sparse(1:d,1:d,1./(diag(D).^2),d,d);
Sigmainv = Dm2-(Dm2*B)/(sparse(1:k,1:k,1,k,k)+B'*(Dm2*B))*(B'*Dm2);
%Sigmainv = inv(B*B'+D.^2);
SigmainvB = Sigmainv*B;
SigmainvD = Sigmainv*D;
Term1 = (phi-repmat(mu',S,1))';
dtdeta = -dtdtheta.*dthetadeta;
dtdtau = dtdeta.*repmat(detadtau',S,1);

dmu = Sigmainv*Term1;
dtau = -(dtdtau').*dmu+(ddt_dthetadtau.*(1./(dtdtheta)))';
dB = zeros(size(B2vechB(B,k),1),S);
dD = zeros(size(Sigma,1),S);
for s = 1:S
    dBtemp = -SigmainvB+Sigmainv*(phi(s,:)'-mu)*(phi(s,:)'-mu)'*SigmainvB;
    dB(:,s) = B2vechB(dBtemp,k);
    dD(:,s) = diag(-SigmainvD+Sigmainv*(phi(s,:)'-mu)*(phi(s,:)'-mu)'*SigmainvD);
end
if VA == 1
    dmuz = [];
    dlogsigmaz = [];
    dC = [];
    dtauz = [];
elseif VA == 2
    sigma2 = exp(2*logsigmaz);
    sigma2rep = repmat(sigma2,1,S);
    muZ = repmat(muz,1,S);
    dlogsigmaz = -1+(1./sigma2rep).*(z-muZ).^2;
    dmuz = (1./sigma2rep).*(z-muZ);
    dC = [];
    dtauz = (ddt_dvarepsdtau.*(1./(dtdvareps)))';
elseif VA == 3
    if isempty(indPositiveC) || isempty(id1) || isempty(id2) || isempty(IT2) || isempty(IT) || isempty(Kcom)
        PositiveC = diag(ones(T,1))+diag(ones(T-1,1),-1);
        indPositiveC = (PositiveC(:)==1);
        [ind1,ind2] = meshgrid(1:T,1:T);
        id1 = ind1(:)';
        id2 = ind2(:)';
        IT2 = sparse(1:T^2,1:T^2,1);
        IT = sparse(1:T,1:T,1);
        Kcom = communtationM(IT);
    end
    Cinvt = inv(C)';
    muZ = repmat(muz,1,S);
    zMmuzt = (z-muZ)';
    zMmuzKronzMmuz = zMmuzt(:,id1).*zMmuzt(:,id2);
    dSigmaInvdC = (IT2+Kcom)*kron(C,IT);
    dC = repmat(Cinvt(:)',S,1)-0.5*zMmuzKronzMmuz*dSigmaInvdC;
    dC = dC(:,indPositiveC);
    dC = dC';
    dmuz = (zMmuzt*C)*C';
    dmuz = dmuz';
    dlogsigmaz = [];
    dtauz = (ddt_dvarepsdtau.*(1./(dtdvareps)))';
end

if strcmp(Transf,'none')
    dtau =[];
    dtauz = [];
end
grad_lambda = [dmu; dB; dD; dmuz; dlogsigmaz;dC;dtau;dtauz'];

