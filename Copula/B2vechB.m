function b = B2vechB(B,K)
ind = zeros(K,K)+tril(ones(K,K));
ind = ind(:);
nozero = B(1:K,1:K);
nozero = nozero(:);
b1 = nozero(ind~=0);
b2 = B(K+1:end,1:K);
b2 = b2(:);
b = [b1;b2];
end