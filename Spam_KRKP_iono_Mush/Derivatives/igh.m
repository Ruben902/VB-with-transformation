function [ z ] = igh( tau,g,h )
%igh : compute the inverse g & h transformation
% inputs: 
% g : real valued vector
% h > 0: positive vector
% tau : vector valued argument
inisize = size(tau);
tau = tau(:);
g = g(:);
h = h(:);

if((any(size(tau)~=size(g)))||any(size(tau)~=size(h))) 
    error('tau,g,h must be the same size in function igh');
end
aprox = 0.00001;
iter_max = 100;
N_obs = size(tau,1);
z = zeros(N_obs,1); 
noConverge = ones(N_obs,1)==1;
tau_est = zeros(N_obs,1);
i=0;
while(sum(noConverge)>0) && (i<=iter_max)
    tau_est(noConverge) = gh(z(noConverge),g(noConverge),h(noConverge));
    dt_dtheta = dgh(z(noConverge),g(noConverge),h(noConverge));
    z(noConverge) = z(noConverge) - ((tau_est(noConverge)-tau(noConverge))./dt_dtheta);
    noConverge = logical(abs(tau_est-tau)>=aprox);
    i =i+1;
end
if sum(noConverge)~=0
   warning(['A number of ' mat2str(sum(noConverge)) 'did not converge on first inversion']) 
end
u0=0;
for i=1:N_obs
    if noConverge(i) || isnan(z(i))
        f1=@(z) (gh(z,g(i),h(i))-tau(i));
        [x,~,~,~] = fzero(f1,u0);
        z(i)=x;
        u0=x;
    else
        u0 = z(i);
    end
end


z = reshape(z,inisize);

    


