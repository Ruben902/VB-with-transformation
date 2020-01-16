%Check David's notes <notes_Ruben_Delta>
% Also check appendix of document <sncopula.tex>
function [const8deri3DP,kronderi3DP] = const8deri3simple(Delta,B,z,epsilon,delta_phi_tilde)
n = size(B,1);
K = size(B,2);
D = diag(Delta);
j = 1:(n+1):n^2;
%k = floor(j/r)+1;
k = floor((j-1)/n)+1;                                                      % Correct Davids notes
l = j-(k-1)*n;
Dm2 = sparse(1:n,1:n,1./(Delta.^2),n,n);
Sigma_phiinv = Dm2-Dm2*B/(sparse(1:K,1:K,1,K,K)+B'*Dm2*B)*B'*Dm2;                                     %Woodbury formula  
Term1 = (B*z+Delta.*epsilon)'*Sigma_phiinv;
Term2 = Term1*D;
Term3 = delta_phi_tilde'*Sigma_phiinv;
Term4 = Term3*D;
kele_Term2 = Term2(k);
lcol_Sigmainv = Sigma_phiinv(:,l);
kcol_D = D(:,k);
lele_Term1 = Term1(l);
const8deri3DP = -(repmat(kele_Term2,n,1).*lcol_Sigmainv)-((Sigma_phiinv*kcol_D).*repmat(lele_Term1,n,1));
kronderi3DP = -(repmat(Term4(k),n,1).*lcol_Sigmainv)-((Sigma_phiinv*kcol_D).*repmat(Term3(l),n,1));










