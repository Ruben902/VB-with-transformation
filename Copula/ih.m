function u_i = ih(w_i,zeta,M,im1,l1,t,m,llim,ulim,copula_choice,q) 
% Evaluate u_i=h_{i,i-1}^{-1} o ... o h_{i,i-q}^{-1}(w_t|u_{i-q|i-1})
% If q=im1, then i-q=i-(i-1)=1.
i = im1+1;
u_i = w_i;

for j = i-q:im1
% Compute index relationship (i,j)-> (k,l1,l2)
s = ceil(j/m);
k = t-s;
l2 = j-m*(s-1);
u2 = reshape(M(j,im1,:),numel(M(j,im1,:)),1);
zeta_s = zeta(l1,l2,k+1,:);
u_i = invh2( u_i, u2, zeta_s ,copula_choice{l1,l2,k+1});
  if u_i>ulim  %check lines to keep things stable
   u_i=ulim;
%   print*,'ih: ulim used at i,j=',i,j
  elseif(u_i<llim)
%   print*,'ih: llim used at i,j=',i,j
   u_i=llim;
  end
end