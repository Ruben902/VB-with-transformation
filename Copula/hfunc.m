function h = hfunc(u1,u2,family,param,arg)
if nargin == 4
   arg = 1; 
end
 if strcmp(family,'t_mix') || strcmp(family,'sjc_mix') || strcmp(family,'gumbel_mix')  || strcmp(family,'skewt_mix')
     h = h_mix_cop(u1,u2,family,param,arg);
 elseif strcmp(family,'Clayton') || strcmp(family,'Gaussian') || strcmp(family,'t')  || strcmp(family,'Gumbel') || strcmp(family,'Galambos')
     h = h_func(u1,u2,param,family);
 elseif strcmp(family,'Gen_gumbel') 
     h = h_gen_gumbel(u1,u2,param(1),param(2));
 elseif strcmp(family,'Gen_galambos') 
     h = h_gen_galambos(u1,u2,param(1),param(2));
 end
