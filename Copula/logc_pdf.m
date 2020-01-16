function l_c = logc_pdf(u1,u2,family,param)
 if strcmp(family,'t_mix') || strcmp(family,'sjc_mix') || strcmp(family,'gumbel_mix')  || strcmp(family,'skewt_mix')
     [~,l_c] = mix_copulapdf([u1(:) u2(:)],param,family);
 elseif strcmp(family,'Clayton') || strcmp(family,'Gaussian')  || strcmp(family,'Gumbel')
     l_c = log(copulapdf(family,[u1(:) u2(:)],param(1)));
 elseif  strcmp(family,'t') 
     l_c = log(copulapdf(family,[u1(:) u2(:)],param(1),param(2)));
 elseif  strcmp(family,'Gen_gumbel') 
     l_c = lnc_gen_gumbel([u1(:) u2(:)],param(1),param(2));
 elseif  strcmp(family,'Galambos') 
     l_c = lnc_galambos(u1(:),u2(:),param(1));
 elseif  strcmp(family,'Gen_galambos')
     l_c = lnc_galambos_gen(u1(:),u2(:),param(1),param(2));
 end

 
