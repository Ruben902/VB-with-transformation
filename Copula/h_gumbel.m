function h = h_gumbel(u,v,alpha)
if alpha<1
   error('alpha must be greater than 1')
end
u = U_transf(u);
v = U_transf(v);
t1 = (-log(u)).^alpha;
t2 = (-log(v)).^alpha;

h = -(exp(-(t1+t2).^(1/alpha)).*(t1+t2).^((1/alpha)-1).*t2)./(v.*log(v));

