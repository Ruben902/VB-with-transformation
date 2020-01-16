function K = communtationM(A)
[m,n] = size(A);
M = reshape(1:m*n,m,n);
Mt = M';
vecMt = Mt(:);
K = sparse(1:n*m,vecMt,1,n*m,n*m);