function [ hval ] = h_gen_gumbel_mix( u1,u2,alpha1,w1,alpha2,w2,w,arg)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
if arg ==1  %v is second argument as in paper so this is u_t|u_{t-1}
%u1 = v, u2 = u 
h1 = h_gen_gumbel(u1,u2,alpha1,w1);
h2 = h_gen_gumbel(u1,1-u2,alpha2,w2);
hval=w*h1+(1-w)*h2;
else
%u1 = u, u2 = v 
h1 = h_gen_gumbel(u1,u2,alpha1,w1);
h2 = (1-h_gen_gumbel(1-u1,u2,alpha2,w2));
hval=w*h1+(1-w)*h2;
end
end
