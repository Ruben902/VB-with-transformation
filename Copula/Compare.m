%%% Comparing Lower Bound Values (Figure 2 in paper)
clear
example = 'AttemptedMurder';
p=3;
VA = 2;
k=3;
family='gumbel_mix';
VBskew = load(['wd_' example '_p' mat2str(p) '_k' mat2str(k) '_' family '_VA' mat2str(VA) '.mat']);
VBnoskew = load([ 'wd_' example '_p' mat2str(p) '_k' mat2str(k) '_' family '_VA' mat2str(VA) '_noskew.mat']);
VBnoskewnok = load(['wd_' example '_p' mat2str(p) '_k' mat2str(0) '_' family '_VA' mat2str(VA) '.mat']);
ini=400;
subplot(1,1,1)
plot(ini:5000,VBskew.VBDAobj.LB(ini:5000))
hold on
plot(ini:5000,VBnoskew.VBDAobj.LB(ini:5000))
plot(ini:5000,VBnoskewnok.VBDAobj.LB(ini:5000))
hold off
title('Lower bound for attempted murder')
xlabel('VB step')
ylabel('Lower Bound')
legend('Copula Approximation','Gaussian Approximation (k=3)','Gaussian Approximation (k=0)','location','SouthEast')
legend boxoff


clear
p=3;
VA = 2;
k=3;
family='gumbel_mix';
y = pickdata('AttemptedMurder');    %Ordinal time series, change to your own data
mdist = selectedDist(y);     %Marginal distribution of y 
[ai,bi] = ai_bi_compute(y,mdist);

%%%Loading MCMC posterior
MCMC = load(['wd_AttemptedMurder_p' mat2str(p) '_gumbel_mix_fitMCMC2.mat'],'G','Lb','Ub','U','ai','bi');
theta = log(MCMC.G./(MCMC.Ub(1)-MCMC.G));
u_mcmc = MCMC.U';
J = size(u_mcmc,2);Ai = repmat(MCMC.ai,1,J);Bi = repmat(MCMC.bi,1,J);
vareps_mcmc = norminv((u_mcmc-Ai).*(1./(Bi-Ai)));

%%% Loading VBDA posterior with asymmetry induced via Yeo-Johnson transformation.
VBskew = load(['wd_AttemptedMurder_p' mat2str(p) '_k' mat2str(k) '_' family '_VA' mat2str(VA) '.mat']);
eta_va_skew = tau2eta(VBskew.VBDAobj .tau,'YJ');
mu_va_skew = VBskew.VBDAobj.mu;
B_skew = VBskew.VBDAobj.B;
D_skew = VBskew.VBDAobj.D;
sigma_va_skew =  sqrt(diag(B_skew*B_skew'+D_skew.^2));
muz_skew = VBskew.VBDAobj.muz;
logsigmaz_skew = VBskew.VBDAobj.logsigmaz;
C_skew = VBskew.VBDAobj.C;
etau_skew = tau2eta(VBskew.VBDAobj.tauu,'YJ');
[thetaVB,~] = q_lambda(mu_va_skew,50000,VBskew.VBDAobj.family,B_skew,D_skew,eta_va_skew,'YJ');
[u,~,vareps] = q_u(ai,bi,muz_skew,logsigmaz_skew,C_skew,50000,VBskew.VBDAobj.VA,etau_skew,'YJ');

%%% Loading VBDA posterior with no asymmetry.
VBnoskew = load(['wd_AttemptedMurder_p' mat2str(p) '_k' mat2str(k) '_' family '_VA' mat2str(VA) '_noskew.mat']);
mu_va_noskew = VBnoskew.VBDAobj.mu;
B_noskew = VBnoskew.VBDAobj.B;
D_noskew = VBnoskew.VBDAobj.D;
sigma_va_noskew =  sqrt(diag(B_noskew*B_noskew'+D_noskew.^2));
muz_noskew = VBnoskew.VBDAobj.muz;
logsigmaz_noskew = VBnoskew.VBDAobj.logsigmaz;
[thetaVBnos,T_theta] = q_lambda(mu_va_noskew,50000,VBnoskew.VBDAobj.family,B_noskew,D_noskew,eta_va_skew*0+1,'YJ');
[u_nos,z,vareps_no] = q_u(ai,bi,VBnoskew.VBDAobj.muz,VBnoskew.VBDAobj.logsigmaz,VBnoskew.VBDAobj.C,50000,VBnoskew.VBDAobj.VA,ones(size(VBnoskew.VBDAobj.muz)),'YJ');

%%% Fitting kernel density to MCMC posteriors for plotting
mcmcdist=cell(size(thetaVB,2),1);
for i =1:size(thetaVB,2)
    mcmcdist{i} = fit_ssvkernel(theta(:,i),200);
end

%%% Plot comparing posterior densities of three most assymetric parameters (Figure 4 in paper).
figure
ParSkewness=skewness(theta);
[~,ix]=sort(abs(ParSkewness));
Ind = ix(end-3:end);
letter = {'a)' 'b)' 'c)' 'd)' 'e)' 'f)' 'g)'};
for j =1:4
    i = Ind(j);
    subplot(2,2,j)
    x = icdf(mcmcdist{i},0.00001):0.01:icdf(mcmcdist{i},0.99999);
    mcmcpdf = pdf(mcmcdist{i},x);
    vbpdf2 = normpdf(x,mu_va_noskew(i),sigma_va_noskew(i));
    plot(x,mcmcpdf,'b','LineWidth',1.5)
    hold on
    tdensplot(mu_va_skew(i),sigma_va_skew(i),eta_va_skew(i),'YJ','r',1.5)   % VB with skewness
    plot(x,vbpdf2,'g','LineWidth',1.5)                                      % VB no skewness
    hold off
    title([letter{j} ' \theta_{' mat2str(i) '}, Skew = ' mat2str(round(100*ParSkewness(i))/100)])
    ylabel(['p(\theta_{' mat2str(i) '}|y)'])
    %xlim([-3 10])
end
LG = legend('MCMC','Copula VA','Gaussian VA','Location','NorthEast');
LG.FontSize =6;
legend boxoff

%%% Plot comparing posterior densities across all parameters (Figure 3 in paper).
figure
vareps_mcmc_mean = mean(vareps_mcmc,2);
vareps_vs_skew = mean(vareps,2);
vareps_vs_no_skew = mean(vareps_no,2);
max_mean = max([vareps_mcmc_mean(:);vareps_vs_skew(:);vareps_vs_no_skew(:)]);
min_mean = min([vareps_mcmc_mean(:);vareps_vs_skew(:);vareps_vs_no_skew(:)]);
vareps_mcmc_std = std(vareps_mcmc');
vareps_vs_skew_std = std(vareps');
vareps_vs_no_skew_std = std(vareps_no');
max_std = max([vareps_mcmc_std(:);vareps_vs_skew_std(:);vareps_vs_no_skew_std(:)]);
min_std = min([vareps_mcmc_std(:);vareps_vs_skew_std(:);vareps_vs_no_skew_std(:)]);
vareps_mcmc_skew = skewness(vareps_mcmc');
vareps_vs_skew_skew = skewness(vareps');
vareps_vs_no_skew_skew = skewness(vareps_no');
max_skew = max([vareps_mcmc_skew(:);vareps_vs_skew_skew(:);vareps_vs_no_skew_skew(:)]);
min_skew = min([vareps_mcmc_skew(:);vareps_vs_skew_skew(:);vareps_vs_no_skew_skew(:)]);

theta_mcmc_mean = mean(theta,1);
theta_vs_skew = mean(thetaVB,1);
theta_vs_no_skew = mean(thetaVBnos,1);
thetamax_mean = max([theta_mcmc_mean(:);theta_vs_skew(:);theta_vs_no_skew(:)]);
thetamin_mean = min([theta_mcmc_mean(:);theta_vs_skew(:);theta_vs_no_skew(:)]);
theta_mcmc_std = std(theta);
theta_vs_skew_std = std(thetaVB);
theta_vs_no_skew_std = std(thetaVBnos);
theta_max_std = max([theta_mcmc_std(:);theta_vs_skew_std(:);theta_vs_no_skew_std(:)]);
theta_min_std = min([theta_mcmc_std(:);theta_vs_skew_std(:);theta_vs_no_skew_std(:)]);
theta_mcmc_skew = skewness(theta);
theta_vs_skew_skew = skewness(thetaVB);
theta_vs_no_skew_skew = skewness(thetaVBnos);
theta_max_skew = max([theta_mcmc_skew(:);theta_vs_skew_skew(:);theta_vs_no_skew_skew(:)]);
theta_min_skew = min([theta_mcmc_skew(:);theta_vs_skew_skew(:);theta_vs_no_skew_skew(:)]);

subplot(3,2,2)
plot(vareps_mcmc_mean,vareps_vs_skew,'o','MarkerSize',4.5,'LineWidth',1.3)
hold on
plot(vareps_mcmc_mean,vareps_vs_no_skew,'+','MarkerSize',5.5,'LineWidth',1.3)
plot([min_mean max_mean],[min_mean max_mean],'k-')
hold off
title('(b) Posterior mean')
xlabel('MCMC')
ylabel('VB')
axis square

subplot(3,2,4)
plot(vareps_mcmc_std,vareps_vs_skew_std,'o','MarkerSize',4.5,'LineWidth',1.3)
hold on
plot(vareps_mcmc_std,vareps_vs_no_skew_std,'+','MarkerSize',5.5,'LineWidth',1.3)
plot([min_std max_std],[min_std max_std],'k-')
hold off
title('(d) Posterior SD')
xlabel('MCMC')
ylabel('VB')
axis square

subplot(3,2,6)
plot(vareps_mcmc_skew,vareps_vs_skew_skew,'o','MarkerSize',4.5,'LineWidth',1.3)
hold on
plot([min_skew max_skew],[min_skew max_skew],'k-')
hold off
title('(f) Posterior skew')
xlabel('MCMC')
ylabel('VB')
axis square

subplot(3,2,1)
plot(theta_mcmc_mean,theta_vs_skew,'o','MarkerSize',4.5,'LineWidth',1.3)
hold on
plot(theta_mcmc_mean,theta_vs_no_skew,'+','MarkerSize',5.5,'LineWidth',1.3)
plot([thetamin_mean thetamax_mean],[thetamin_mean thetamax_mean],'k-')
hold off
title('(a) Posterior mean')
xlabel('MCMC')
ylabel('VB')
LG=legend('Copula', 'Gaussian','MCMC','Location','SouthEast');
LG.FontSize=6;
axis square

subplot(3,2,3)
plot(theta_mcmc_std,theta_vs_skew_std,'o','MarkerSize',4.5,'LineWidth',1.3)
hold on
plot(theta_mcmc_std,theta_vs_no_skew_std,'+','MarkerSize',5.5,'LineWidth',1.3)
plot([theta_min_std theta_max_std],[theta_min_std theta_max_std],'k-')
hold off
title('(c) Posterior SD')
xlabel('MCMC')
ylabel('VB')
axis square

subplot(3,2,5)
plot(theta_mcmc_skew,theta_vs_skew_skew,'o','MarkerSize',4.5,'LineWidth',1.3)
hold on
plot([theta_min_skew theta_max_skew],[theta_min_skew theta_max_skew],'k-')
hold off
title('(e) Posterior skew')
xlabel('MCMC')
ylabel('VB')
axis square

