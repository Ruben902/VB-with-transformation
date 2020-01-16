function mix_cop_plot(Param,u_lim,scale,family,Log)
cc= 25;
v = linspace(u_lim(1),u_lim(2),cc);
[U1,U2] = meshgrid(v,v);
if nargin ==4
    Log = false;
end
w = Param(5);
 if strcmp(family,'t_mix')
   F = w*copulapdf('t',[U1(:) U2(:)],Param(1),Param(2))+(1-w)*copulapdf('t',[1-U1(:) U2(:)],Param(3),Param(4));
 elseif strcmp(family,'t')
   F = copulapdf('t',[U1(:) U2(:)],Param(1),Param(2));
 elseif strcmp(family,'sjc_mix')
   F = exp(lnc_sjc_mix2(U1(:),U2(:),[Param w]));
 elseif strcmp(family,'gumbel_mix')
   F = exp(lnc_gen_gum_mix([U1(:) U2(:)],[Param(1:4) w]));
    elseif strcmp(family,'clayton_mix')
        F = exp(lnc_gen_clay_mix([U1(:) U2(:)],[Param(1:4) w]));
 elseif strcmp(family,'skewt_mix')
 
 end
 Zlabel = 'c(u_{t-1},u_t)';
 Z = F;
 if Log == true
    Z = log(F); 
    Zlabel = 'log(c(u_{t-1},u_t))';
 end
 ZZ = scale*reshape(Z,cc,cc);
    surf(U1,U2,ZZ);
xlabel('u_{t-1}')
ylabel('u_t')
zlabel(Zlabel)