function [ out ] = dgh_h( z,g,h )
%dgh_g : compute the first differential of g & h transformation wrt g
% inputs: 
% g : real valued array
% h > 0: positive array
% z : array valued argument

if((any(size(z)~=size(g)))||any(size(z)~=size(h))) 
    error('z,g,h must be the same size in function dgh_h');
end

A=exp((h.*z.^2)/2);
k1=find(g==0);

if(isempty(k1))
   out = (z.^2.*A.*(exp(g.*z) - 1))./(2*g); 
else
    out=zeros(size(z));
    k2=find(g~=0);
    z1=z(k1); A1=A(k1);
    g2=g(k2); z2=z(k2); A2=A(k2);

    out(k1)=((z1.^3).*A1)/2;
    out(k2) = (z2.^2.*A2.*(exp(g2.*z2) - 1))./(2*g2);
end

