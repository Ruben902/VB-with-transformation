function [ out ] = digh( tau,g,h,z)
%digh: compute the first derivative of the inverse g & h transformation 
% inputs: 
% g : real valued array
% h > 0: positive array
% tau : array valued argument

if((any(size(tau)~=size(g)))||any(size(tau)~=size(h))) 
    error('tau,g,h must be the same size in function digh');
end
if nargin ==3
    z=igh(tau,g,h);
end
out=ones(size(z))./dgh(z,g,h);

end

