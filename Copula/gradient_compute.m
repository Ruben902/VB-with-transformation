function [H,LogLhat,log_theta,log_q_lambda,grad_lambda] = gradient_compute(theta,Gamma,u,z,vareps,ai,bi,family,lambda,k,VA,ct,Transf)
d = size(theta,2);       %Number of copula parameters
T = size(ai,1);
[mu,Sigma,~,~,muz,logsigmaz,C,tau,tauu] = lambda2musigma(lambda,d,T,k,VA,Transf);
LogLhat = mix_cop_DVine_llParalel_u(Gamma,ai,bi,family,u);
S = length(LogLhat);
Ai = repmat(ai,1,S);
Bi = repmat(bi,1,S);
switch Transf
    case 'YJ'
        eta = tau2eta(tau,Transf);
        etau = tau2eta(tauu,Transf);
        phi = tYJ(theta,repmat(eta',S,1));
        logJacob = log(dtYJ_dtheta(theta,repmat(eta',S,1)));
        if VA>1
            logJacob_u = log(dtYJ_dtheta(vareps,repmat(etau,1,S)));
        else
            logJacob_u=[];
        end
    case 'GH'
        eta = tau2eta(tau,Transf);
        etau = tau2eta(tauu,Transf);
        phi = gh(theta,repmat(eta(:,1)',S,1),repmat(eta(:,2)',S,1));
        logJacob = log(dgh(theta,repmat(eta(:,1)',S,1),repmat(eta(:,2)',S,1)));
        if VA>1
            logJacob_u = log(dgh(vareps,repmat(etau(:,1),1,S),repmat(etau(:,2),1,S)));
        else
            logJacob_u=[];
        end
    case 'none'
        phi = theta;
        logJacob=0;
        logJacob_u=0;
end
log_theta = logpri(theta,family,tau2eta(tau,Transf),Transf);
log_q_lambda =  log_mvnpdf(phi',repmat(mu,1,S),Sigma)'+sum(logJacob,2)';

if  VA == 1
    log_q_u = repmat(sum(-log(bi-ai)),1,S);
elseif VA ==2
    sigmaz2 = repmat(exp(2*logsigmaz),1,S);
    muZ = repmat(muz,1,S);
    log_q_u = sum(-0.5*log(2*pi*sigmaz2)-0.5*(1./(sigmaz2)).*(z-muZ).^2-log(Bi-Ai)+0.5*log(2*pi)+0.5*vareps.^2)+sum(logJacob_u);
elseif VA == 3
    muZ = repmat(muz,1,S);
    zMmuztC = (z-muZ)'*C;
    log_q_u = -(T/(2))*log(2*pi)+log(abs(det(C)))-0.5*sum(zMmuztC.*zMmuztC,2)'+sum(-log(Bi-Ai)+0.5*log(2*pi)+0.5*vareps.^2)+sum(logJacob_u);
end

log_q_lambda = log_q_lambda + log_q_u;
grad_lambda = grad_lambda_logq(theta,z,vareps,lambda,k,T,VA,Transf);

nLambda = size(grad_lambda,1);
H = [];
if ~isempty(ct)
    H = zeros(nLambda,1);
    for i = 1:nLambda
        H(i) = mean((LogLhat+log_theta-log_q_lambda-ct(i)).*grad_lambda(i,:));
    end
end