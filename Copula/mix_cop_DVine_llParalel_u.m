%xparam: Theta, the parameters whose posterior distiburion we want tocompute. This is of dimension (S)X(#parameters), where S is the number of points for gradient computation
%ai = F(xi-1), and bi = F(xi).
%family: copula famility to be used for dependence
%N: Number of particles used to copute likelihood approximation
function ll = mix_cop_DVine_llParalel_u(xparam,ai,bi,family,u)
%T = size(ai,1);
if strcmp(family,'gumbel_mix')
numparpercop = 5;
end
if strcmp(family,'t')
numparpercop = 2;
end
lag = (size(xparam,2))/numparpercop;
if mod(lag,1)~=0
    error('Number of Lags should be discrete')
end
%S = size(xparam,1);
%Ai = repmat(ai,1,S);
%Bi = repmat(bi,1,S);
%u = (Bi-Ai).*rand(size(Ai))+Ai;
xparam1 = xparam(:,1:numparpercop*lag);
S1 = size(u,2);
hini1 = u(1:(end-1),:);
hini2 = u(2:end,:);
ll = zeros(1,S1);
parfor s = 1:S1
    h1 = hini1(:,s);
    h2 = hini2(:,s);
    Xparam = xparam1(s,:);
    for i = 1:lag
        par_temp = Xparam((numparpercop*(i-1)+1):(numparpercop*i));
        [~,temp] = mix_copulapdf([h1 h2],par_temp,family);
        ll(s) = ll(s) + sum(temp);
        if i<lag
        a = h_mix_cop(h1(1:end-1),h2(1:end-1),family,par_temp,2); %ut-1/ut
        b = h_mix_cop(h2(2:end),h1(2:end),family,par_temp,1); %ut/u{t-1}
        h1 = a;
        h2 = b;
        end
    end
end