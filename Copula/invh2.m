function u = invh2(w,v,Param,family)
iter_max = 1000;
aprox = 0.0000001;
w = w(:);
v = v(:);
n = size(w,1);
a = zeros(n,1);
b = ones(n,1);
Fa = zeros(n,1);
Fb = ones(n,1);
Fx = ones(n,1);
i = 1;
Conditions = ones(n,1)==1;   % Created to make inversions in parallel
while(sum(Conditions)>0 && i<=iter_max)
    u_est = a+(b-a).*(w-Fa)./(Fb-Fa);
    Fx(Conditions) = hfunc(u_est(Conditions),v(Conditions),family,Param);
    Conditions1 = logical(((Fx-w)>0).*(abs(Fx-w)>=aprox));
    Conditions2 = logical(((Fx-w)<0).*(abs(Fx-w)>=aprox));
    b(Conditions1) = u_est(Conditions1);
    Fb(Conditions1) = Fx(Conditions1);
    a(Conditions2) = u_est(Conditions2);
    Fa(Conditions2) = Fx(Conditions2);
    Conditions = logical(Conditions1 + Conditions2);
    i =i+1;
end
u = u_est;

end


















