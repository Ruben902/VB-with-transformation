function [ out ] = dgh_g( z,g,h )
%dgh_g : compute the first differential of g & h transformation wrt g
% inputs: 
% g : real valued array
% h > 0: positive array
% z : array valued argument

if((any(size(z)~=size(g)))||any(size(z)~=size(h))) 
    error('z,g,h must be the same size in function dgh_g');
end

A=exp((h.*z.^2)/2);
k1=find(g==0);
    
if(isempty(k1))
    out = (z.*A.*exp(g.*z))./g - (A.*(exp(g.*z) - 1))./(g.^2);
else
    out=zeros(size(z));
    k2=find(g~=0);
    g2=g(k2); z2=z(k2); A2=A(k2);

    out(k2) = (z2.*A2.*exp(g2.*z2))./g2 - (A2.*(exp(g2.*z2) - 1))./(g2.^2);
end