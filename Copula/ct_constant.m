function ct = ct_constant(LogLhat,log_theta,log_q_lambda,grad_lambda)
nLambda = size(grad_lambda,1);
ct = zeros(nLambda,1);
for i = 1:nLambda
    v = cov((LogLhat+log_theta-(log_q_lambda)).*grad_lambda(i,:),grad_lambda(i,:));
    ct(i) = v(1,2)/v(2,2);
end


