% Set <diff> to true to provide the gradient as the output
% Set <diff> to false to provide log h(theta) as the output
function val = log_normal(theta, m, K,diff)
%


[D D2] = size(K);

if D == D2  % the covariance matrix is full
    %
    L = chol(K)';
    d = theta - m;
    ld = L\d;
    if diff ==false
        val = - 0.5*D*log(2*pi) - sum(log(diag(L)))  - 0.5*(ld'*ld);
    else
        val = - (L'\ld);
    end
    %
else % the covariance matrix is diagonal
    %
    d = theta - m;
    if diff ==false
        val = - 0.5*D*log(2*pi) - 0.5*sum(log(K))  - 0.5*sum((d.^2)./K);
    else
        val = - d./K;
    end
    %
end


