function h = h_mix_cop(u1,u2,family,param,arg)
if nargin ==4
   arg =1; 
end
%arg, indicates if h function 1 or 2 is needed, only necesary for Archemedean copulas
u1 = U_transf(u1(:));
u2 = U_transf(u2(:));
if strcmp(family,'t_mix') 
    h = h_tcop_mix(u1,u2,param(1),param(2),param(3),param(4),param(5));
elseif strcmp(family,'sjc_mix') 
    h = h_sjc_mix( u1,u2,param(1),param(2),param(3),param(4),param(5));
elseif strcmp(family,'gumbel_mix') 
    h = h_gen_gumbel_mix( u1,u2,param(1),param(2),param(3),param(4),param(5),arg);
elseif strcmp(family,'t')
    h = h_tcop(u1,u2,param(1),param(2)); 
end

           

