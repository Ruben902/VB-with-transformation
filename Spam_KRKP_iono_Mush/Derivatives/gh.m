function [ tau ] = gh( z,g,h )
%gh : compute the g & h transformation of z
% inputs: 
% g : real valued array
% h > 0: positive array
% z : array valued argument

if((any(size(z)~=size(g)))||any(size(z)~=size(h))) 
    error('z,g,h must be the same size in function gh');
end

A=exp((h.*z.^2)/2);
k1=find(g==0);

if(isempty(k1))
    t2=(exp(z.*g)-1)./g;
    tau=t2.*A;
else
    tau=zeros(size(z));
    k2=find(g~=0);
    tau(k1)=z(k1).*A(k1);
    t2=(exp(z(k2).*g(k2))-1)./g(k2);
    tau(k2)=t2.*A(k2);
end

end

