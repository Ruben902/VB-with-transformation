function [ out ] = ddigh( z,g,h,zinv)
%ddigh: second differential of the g-&-h function inverse wrt to its argument
% inputs:
% g : real valued array
% h > 0: positive array
% z : array valued argument
if((any(size(z)~=size(g)))||any(size(z)~=size(h)))
    error('z,g,h must be the same size in function ddgh');
end
if nargin==3
    zinv = igh(z,g,h);
end
term1 = ddgh( zinv,g,h );
term2 = -(dgh( zinv,g,h ).^(-3));
out = term2.*term1;