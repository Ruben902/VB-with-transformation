function log_pdf = log_mvnpdf(mu,mu_star,Sigma) 
m = size(mu,1);
log_pdf = -0.5*m*log(2*pi)-0.5*log(det(Sigma))-0.5*sum(((mu - mu_star)'/Sigma).*(mu - mu_star)',2);  % Smart way to do it, dont forget, trust is correct, I have  checked many times



