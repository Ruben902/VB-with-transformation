function [pdf,log_pdf] = mix_copulapdf(U,Param,family)
U =  U_transf(U);
if strcmp(family,'t_mix')
    w = Param(5);
    %log_pdf = log(w*copulapdf('t',U,Param(1),Param(2))+(1-w)*copulapdf('t',[(1-U(:,1)) U(:,2)],Param(3),Param(4)));
    logc_1 = log(copulapdf('t',U,Param(1),Param(2)));
    logc_2 = log(copulapdf('t',[(1-U(:,1)) U(:,2)],Param(3),Param(4)));
    log_pdf = log(w)+logc_1+log(1+expma(log(1-w)+logc_2-log(w)-logc_1));
 elseif strcmp(family,'sjc_mix')
    log_pdf = lnc_sjc_mix2(U(:,1),U(:,2),Param);
 elseif strcmp(family,'gumbel_mix')
    log_pdf = lnc_gen_gum_mix(U,Param);
elseif strcmp(family,'skewt_mix')
  %invh = @(w_t,u_t_1) 
elseif strcmp(family,'t')
   log_pdf = log(copulapdf('t',U,Param(1),Param(2)));
end
pdf = expma(log_pdf);

