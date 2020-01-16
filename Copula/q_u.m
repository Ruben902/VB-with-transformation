function [u,z,vareps] = q_u(ai,bi,muz,logsigmaz,C,S,VA,etau,Transf)
T = size(ai,1);
Ai = repmat(ai,1,S);
Bi = repmat(bi,1,S);

switch Transf
    case 'YJ'
        transf_inv = @(x,eta,ignore) tYJi(x,eta);
        etau = [etau etau];
    case 'GH'
        transf_inv = @(x,g,h) igh(x,g,h);
    case 'none'
        transf_inv = @(x,ignore1,ignore2) x;
        etau = [NaN NaN];   %This term is ignored in the inversion
end

if VA == 1
    z = [];
    u = rand(T,S).*(Bi-Ai)+Ai;
    vareps=norminv((u-Ai)./(Bi-Ai));
elseif VA == 2
    z = normrnd(0,1,T,S).*repmat(exp(logsigmaz),1,S)+repmat(muz,1,S);
    vareps = transf_inv(z,repmat(etau(:,1),1,S),repmat(etau(:,2),1,S));
    u = normcdf(vareps).*(Bi-Ai)+Ai;
elseif VA == 3
    Cinvt = inv(C)';
    z = repmat(muz,1,S) + Cinvt*normrnd(0,1,T,S);
    vareps = transf_inv(z,repmat(etau(:,1),1,S),repmat(etau(:,2),1,S));
    u = normcdf(vareps).*(Bi-Ai)+Ai;
end


