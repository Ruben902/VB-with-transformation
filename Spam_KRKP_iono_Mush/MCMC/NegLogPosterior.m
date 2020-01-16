function [f] = NegLogPosterior(w,X,y,vInv)
[n,p] = size(X);
Xw = X*w;
yXw = y.*Xw;
sig = 1./(1+exp(-yXw));
f = sum(mylogsumexp([zeros(n,1) -yXw])) + (1/2)*w'*vInv*w;
end