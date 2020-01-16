function [] = plot_posterior(ind,eta,F_mu,F_B,F_d,alphaphi,Transf)
%Only first four posterior indicated by ind are plotted
ind = ind(1:min(4,length(ind)));
theta = snrand(10000,eta,F_mu,F_B,F_d,alphaphi,Transf);

for i = 1:length(ind)
    subplot(2,ceil(length(ind)/2),i)
    dist = fit_ssvkernel(theta(ind(i),:));
    lim = icdf(dist,[0.0001 1-0.0001]);
    llim = lim(1);
    ulim = lim(2);
    x = linspace(llim,ulim,100);
    plot(x,pdf(dist,x))
    ylabel(['f(\theta_{' mat2str(ind(i)) '}|y)'])
    xlabel(['\theta_{' mat2str(ind(i)) '}'])
end

