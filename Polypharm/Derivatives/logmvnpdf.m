function [logp] = logmvnpdf(x,mu,Sigma,CholSigma)
if nargin ==3
    CholSigma = [];
end
    
    
% outputs log likelihood array for observations x  where x_n ~ N(mu,Sigma)
% x is NxD, mu is 1xD, Sigma is DxD

[N,D] = size(x);
const = -0.5 * D * log(2*pi);

xc = bsxfun(@minus,x,mu);

term1 = -0.5 * sum((xc / Sigma) .* xc, 2); % N x 1
term2 = const - 0.5 * logdet(Sigma,CholSigma);    % scalar
logp = term1' + term2;

end

function y = logdet(A,U)
if isempty(U)
    U = chol(A);
end
y = 2*sum(log(abs(diag(U))));
end