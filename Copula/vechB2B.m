function B = vechB2B(b,n,K)
B = zeros(n,K);
Dd = DuplicationM(K);
B(1:K,1:K) = tril(reshape(Dd*b(1:(0.5*K*(K+1))),K,K));
B(K+1:end,:) = reshape(b((0.5*K*(K+1))+1:end),n-K,K);
%B = vech(B)
end