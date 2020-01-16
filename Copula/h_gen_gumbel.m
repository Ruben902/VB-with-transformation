function h = h_gen_gumbel(u,v,tau,w)
if abs(tau)>1
   error('abs(tau) must be less than 1')
end

alpha = copulaparam('gumbel',abs(tau));

if tau>=0
    h = w*h_gumbel(u,v,alpha)+(1-w)*(1-h_gumbel(1-u,1-v,alpha)); 
elseif tau<0
    h = w*h_gumbel(u,1-v,alpha)+(1-w)*(1-h_gumbel(1-u,v,alpha)); 
end