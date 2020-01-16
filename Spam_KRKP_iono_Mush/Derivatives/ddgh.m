function [ out ] = ddgh( z,g,h )
%ddgh: second differential of the g-&-h function wrt to its argument
% inputs: 
% g : real valued array
% h > 0: positive array
% z : array valued argument

if((any(size(z)~=size(g)))||any(size(z)~=size(h))) 
    error('z,g,h must be the same size in function ddgh');
end

A=exp((h.*z.^2)/2);
B=(exp(g.*z) - 1)./g;
k1=find(g==0);

if(isempty(k1))
    C=exp(g.*z);
    out = g.*A.*C + h.*A.*B + (2*h).*(z.*A.*C) + (h.^2).*(z.^2).*A.*B;
else
    out=zeros(size(z));
    k2=find(g~=0);
    z1=z(k1); h1=h(k1); A1=A(k1);
    g2=g(k2); z2=z(k2); h2=h(k2); A2=A(k2); B2=B(k2);
    
    out(k1) = (h1.^2).*(z1.^3).*A1 + (3*h1).*z1.*A1;
    C2=exp(g2.*z2);
    out(k2) = g2.*A2.*C2 + h2.*A2.*B2 + (2*h2).*(z2.*A2.*C2) + (h2.^2).*(z2.^2).*A2.*B2;

end