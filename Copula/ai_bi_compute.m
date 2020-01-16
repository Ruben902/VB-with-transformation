%Tis code is only valid for counts greater or equal to zero.
function [ai,bi] = ai_bi_compute(y,mdist,type)
if nargin ==2
    type = 'discrete';
end
switch type
    case 'discrete'                         % Empirical distribution function
        bi = cdf(mdist,y(:));
        ai = cdf(mdist,y(:)-1);
    case 'continous'
        bi = cdf(mdist,y(:));           % ssv kernel
        ai = bi;
    case 'continuous'
        bi = cdf(mdist,y(:));           % ssv kernel
        ai = bi;
end