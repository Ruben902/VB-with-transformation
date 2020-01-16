function [ out ] = dgh( z,g,h )
%dgh: differential of the g-&-h function wrt to its argument
% inputs: 
% g : real valued array
% h > 0: positive array
% z : array valued argument

if((any(size(z)~=size(g)))||any(size(z)~=size(h))) 
    error('z,g,h must be the same size in function dgh');
end

A=exp((h.*z.^2)/2);
k1=find(g==0);

if(isempty(k1))
    out = A.*exp(g.*z) + (h.*z.*A.*(exp(g.*z) - 1))./g;
else
    out=zeros(size(z));
    k2=find(g~=0);
    z1=z(k1); h1=h(k1); A1=A(k1);
    g2=g(k2); z2=z(k2); h2=h(k2); A2=A(k2);

    out(k1) = A1 + h1.*(z1.^2).*A1;
    out(k2) = A2.*exp(g2.*z2) + (h2.*z2.*A2.*(exp(g2.*z2) - 1))./g2;

end

