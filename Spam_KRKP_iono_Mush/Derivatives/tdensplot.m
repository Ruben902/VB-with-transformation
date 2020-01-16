function []=tdensplot(mu,sigma,eta,transf)
phiup = mu+5*sigma;
phimin = mu-5*sigma;
switch transf
    case 'YJ'
        jacobian =  @(x) dtYJ_dtheta(x,eta);
        tfunc = @(x) tYJ(x,eta);
        thetaup = tYJi(phiup,eta);
        thetamin = tYJi(phimin,eta);
    case 'GH'
        jacobian =  @(x) dgh(x,repmat(eta(1),size(x)),repmat(eta(2),size(x)));
        tfunc = @(x) gh(x,repmat(eta(1),size(x)),repmat(eta(2),size(x)));
        thetaup = igh(phiup,repmat(eta(1),size(2)),repmat(eta(2),size(2)));
        thetamin = igh(phimin,repmat(eta(1),size(2)),repmat(eta(2),size(2)));
end
X = thetamin:0.001:thetaup;
func = @(x) jacobian(x).*normpdf(tfunc(x),mu,sigma);
fx = func(X);
plot(X,fx,'LineWidth',1.5)