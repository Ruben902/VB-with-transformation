% Set <diff> to true to provide the gradient as the output
% Set <diff> to false to provide log h(theta) as the output
function val = log_glmm(theta, m, K, b_n,diff)
% m is the mean
% K is the covariance
% We consider only diagonal covariance here

b_size = length(m) - b_n;

W = exp(theta(end));
if diff == false
    val = - b_size*log(W) - 0.5/W^2*sum(theta(b_n:(end-1)).^2) - 1/(2*K(1))*sum(theta(1:(b_n-1)).^2) - 1/(2*K(end))*(theta(end).^2);
else
    dg_b = - theta(b_n:(end-1))/W^2;
    dg_beta = -1/(K(2))*theta(1:(b_n-1));
    dg_xi = - b_size + sum(theta(b_n:(end-1)).^2)/W^2 -  1/(K(end))*theta(end);
    val = [dg_beta ; dg_b  ; dg_xi];
end


