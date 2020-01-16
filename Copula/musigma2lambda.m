function [lambda,B,D,tau,tauu] = musigma2lambda(mu,Sigma,muz,logsigmaz,C,k,eta,etau,Transf)
d = size(Sigma,1);
B = zeros(d,k)+0.001;
D = diag(sqrt(diag(Sigma)));
b = B2vechB(B,k);
T = size(muz,1);
indPositiveC = [];
if ~isempty(C)
    PositiveC = diag(ones(T,1))+diag(ones(T-1,1),-1);
    indPositiveC = (PositiveC(:)==1);
end
if strcmp(Transf,'none')
    tau =[];
    tauu =[];
else
    tau = eta2tau(eta,Transf);
    tauu = eta2tau(etau,Transf);
end
lambda = [mu;b;sqrt(diag(Sigma));muz;logsigmaz;C(indPositiveC);tau(:);tauu(:)];
