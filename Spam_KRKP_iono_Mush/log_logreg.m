% Set <diff> to true to provide the gradient as the output
% Set <diff> to false to provide log h(theta) as the output
function val = log_logreg(theta, X, Y, T,diff)
F = X*theta;

YF = -Y.*F;

m = max(0,YF);
if diff == false
    val = - sum( m + log( exp(-m) + exp( YF - m )) );
else
    S = sigmoid(F);
    val = X'*(T - S);
end



