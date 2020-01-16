function val=lnc_gen_gumbel(U,tau,w)
%   The parameter for this copula is the kendalls tau, and the weight for 0
%   and 180 rotated copula
if abs(tau)>1
   error('abs(tau) must be less than 1')
end

alpha = copulaparam('gumbel',abs(tau));
if tau>=0
    val = log(w*copulapdf('gumbel',U,alpha)+(1-w)*copulapdf('gumbel',[1-U(:,1) 1-U(:,2)],alpha));
elseif tau<0
    val = log(w*copulapdf('gumbel',[U(:,1) 1-U(:,2)],alpha)+(1-w)*copulapdf('gumbel',[1-U(:,1) U(:,2)],alpha));
end
end