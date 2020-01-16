function log_pdf=lnc_gen_gum_mix(U,Param)
% pdf copula of a mix of 2 bidirectional gumbels
    logc_1 = lnc_gen_gumbel(U,Param(1),Param(2));
    logc_2 = lnc_gen_gumbel([1-U(:,1) U(:,2) ],Param(3),Param(4));
    w = Param(5);
    log_pdf = log(w)+logc_1+log(1+expma(log(1-w)+logc_2-log(w)-logc_1));
end