function [H1,H2] = hfunctions_t_multiple(u,xparam,family,numparpercop)
lag = length(xparam)/numparpercop;
H1 =  rand(length(u),lag);
H2 =  rand(length(u),lag);
h1 = u(1:(end-1));
h2 = u(2:end);

for i = 1:lag
    H1((i+1):end,lag-i+1) = h1;   %     [ut-4|ut-1, ut-3|ut-1, ut-2|ut-1, ut-1]      %lag = 4
    H2((i+1):end,lag-i+1) = h2;    %    [ut|ut-3  , ut|ut-2  , ut|ut-1  , ut]
    par_temp = xparam((numparpercop*(i-1)+1):(numparpercop*i));
    a = h_mix_cop(h1(1:end-1),h2(1:end-1),family,par_temp,2); %u1/u2
    b = h_mix_cop(h2(2:end),h1(2:end),family,par_temp,1); %u2/u1
    h1 = a;
    h2 = b;
end
