clear
rng(19900209);
workingDir = pwd;
addpath(genpath([workingDir '/Derivatives']));
S =10000;
MF = false;
Mfchar = {'' '_MF'};
A1 = load('Poly_NoSkewVB_MF.mat');
A3 = load('Poly_NoSkewVB.mat');
A4 = load('Poly_SkewVB.mat');
A5 = load('Poly_NoSkewYJVB.mat');
A6 = load('Poly_SkewYJVB.mat');
A7 = load('Poly_NoSkewiGHVB.mat');
A8 = load('Poly_SkewiGHVB.mat');


MCMC = table2array(readtable('mcmc_polypharm_v2.txt'));
MCMC_all = table2array(readtable('mcmc_polypharm.txt'))';
MCMC_all= MCMC_all([1:7 9:end 8],:);

thetaA1 = snrand(S,A1.eta,A1.F_mu,A1.F_B,A1.F_d,A1.alphaphi,A1.Transf);
thetaA3 = snrand(S,A3.eta,A3.F_mu,A3.F_B,A3.F_d,A3.alphaphi,A3.Transf);
thetaA4 = snrand(S,A4.eta,A4.F_mu,A4.F_B,A4.F_d,A4.alphaphi,A4.Transf);
thetaA5 = snrand(S,A5.eta,A5.F_mu,A5.F_B,A5.F_d,A5.alphaphi,A5.Transf);
thetaA6 = snrand(S,A6.eta,A6.F_mu,A6.F_B,A6.F_d,A6.alphaphi,A6.Transf);
thetaA7 = snrand(S,A7.eta,A7.F_mu,A7.F_B,A7.F_d,A7.alphaphi,A7.Transf);
thetaA8 = snrand(S,A8.eta,A8.F_mu,A8.F_B,A8.F_d,A8.alphaphi,A8.Transf);


%%%The marginal posteriors of the model parameters for the mixed logistic regression model to the polypharmacy dataset (Figure 6 in paper)
Letters = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)'};
Labels = {'\beta_0' '\beta_{gender}' '\beta_{race}' '\beta_{age}' '\beta_{M1}' '\beta_{M2}' '\beta_{M3}' '\beta_{IM}' '\zeta'};
s = 0;
for i  = [1:8 509]
    s=s+1;
    distvb = fit_ssvkernel(thetaA3(i,:));
    distvb3 = fit_ssvkernel(thetaA5(i,:));
    distvbMF = fit_ssvkernel(thetaA1(i,:));
    if i <9
        subplot(3,3,i)
        distmcmc = fit_ssvkernel(MCMC(:,i));
    else
        subplot(3,3,9)
        distmcmc = fit_ssvkernel(MCMC(:,9));
        s = 9;
    end
    Lims = [min(icdf(distmcmc,0.0005),icdf(distvb,0.0005))    max(icdf(distmcmc,1-0.0005),icdf(distvb,1-0.0005))] ;
    xx = Lims(1):0.005:Lims(2);
    pdfmcmc = pdf(distmcmc,xx);
    pdfvb = pdf(distvb,xx);
    pdfvb3 = pdf(distvb3,xx);
    pdfvbMF = pdf(distvbMF,xx);
    plot(xx,pdfmcmc,'LineWidth',1.3)
    hold on
    plot(xx,pdfvb,'LineWidth',1.3)
    plot(xx,pdfvb3,'LineWidth',1.3)
    plot(xx,pdfvbMF,'--','LineWidth',1.3)
    hold off
    title(Letters{s})
    xlabel(Labels{s})
    drawnow
end
legend('MCMC','Normal','Copula','Normal-MF')
legend boxoff
set(gcf, 'Position',  [351,61,1459,912])
saveas(gcf,'INLA/YJ_PAR_densities_plus_MF.fig')

figure


Transf = 'YJ';
Skew = skewness(MCMC_all');
theta1Skew = skewness(thetaA3');
theta2Skew = skewness(thetaA4');
if strcmp(Transf,'YJ')
    theta3Skew = skewness(thetaA5');
    theta4Skew = skewness(thetaA6');
else
    theta3Skew = skewness(thetaA7');
    theta4Skew = skewness(thetaA8');
end

McmcMean = mean(MCMC_all');
theta1Mean = mean(thetaA3');
theta2Mean = mean(thetaA4');
if strcmp(Transf,'YJ')
    theta3Mean = mean(thetaA5');
    theta4Mean = mean(thetaA6');
else
    theta3Mean = mean(thetaA7');
    theta4Mean = mean(thetaA8');
end

McmcSD = std(MCMC_all');
theta1SD = std(thetaA3');
theta2SD = std(thetaA4');
if strcmp(Transf,'YJ')
    theta3SD = std(thetaA5');
    theta4SD = std(thetaA6');
else
    theta3SD = std(thetaA7');
    theta4SD = std(thetaA8');
end


MeanLims=[min([McmcMean(:)' theta1Mean(:)' theta2Mean(:)' theta3Mean(:)' theta4Mean(:)']) max([McmcMean(:)' theta1Mean(:)' theta2Mean(:)' theta3Mean(:)' theta4Mean(:)'])];
SDLims=[min([McmcSD(:)' theta1SD(:)' theta2SD(:)' theta3SD(:)' theta4SD(:)']) max([McmcSD(:)' theta1SD(:)' theta2SD(:)' theta3SD(:)' theta4SD(:)'])];
SKEWLims=[min([Skew(:)' theta1Skew(:)' theta2Skew(:)' theta3Skew(:)' theta4Skew(:)']) max([Skew(:)' theta1Skew(:)' theta2Skew(:)' theta3Skew(:)' theta4Skew(:)'])];


theta1Skew = skewness(thetaA3');
theta2Skew = skewness(thetaA4');
if strcmp(Transf,'YJ')
    theta3Skew = skewness(thetaA5');
    theta4Skew = skewness(thetaA6');
else
    theta3Skew = skewness(thetaA7');
    theta4Skew = skewness(thetaA8');
end

McmcMean = mean(MCMC_all');
theta1Mean = mean(thetaA3');
theta2Mean = mean(thetaA4');
if strcmp(Transf,'YJ')
    theta3Mean = mean(thetaA5');
    theta4Mean = mean(thetaA6');
else
    theta3Mean = mean(thetaA7');
    theta4Mean = mean(thetaA8');
end

McmcSD = std(MCMC_all');
theta1SD = std(thetaA3');
theta2SD = std(thetaA4');
if strcmp(Transf,'YJ')
    theta3SD = std(thetaA5');
    theta4SD = std(thetaA6');
else
    theta3SD = std(thetaA7');
    theta4SD = std(thetaA8');
end
subplot(3,4,1)
plot(McmcMean,theta1Mean,'o','MarkerSize',4.5,'LineWidth',1.3)
title('(a) Normal')
hold on
plot([MeanLims(1) MeanLims(2)],[MeanLims(1) MeanLims(2)],'k-')
hold off
xlabel('Posterior Mean')
ylabel('VA mean')
subplot(3,4,2)
plot(McmcMean,theta2Mean,'o','MarkerSize',4.5,'LineWidth',1.3)
title('(b) Skew Normal')
hold on
plot([MeanLims(1) MeanLims(2)],[MeanLims(1) MeanLims(2)],'k-')
hold off
xlabel('Posterior Mean')
ylabel('VA mean')
subplot(3,4,3)
plot(McmcMean,theta3Mean,'o','MarkerSize',4.5,'LineWidth',1.3)
title(['(c) Gaussian Copula (' Transf ')'])
hold on
plot([MeanLims(1) MeanLims(2)],[MeanLims(1) MeanLims(2)],'k-')
hold off
xlabel('Posterior Mean')
ylabel('VA mean')
subplot(3,4,4)
plot(McmcMean,theta4Mean,'o','MarkerSize',4.5,'LineWidth',1.3)
title(['(d) Skew Normal Copula (' Transf ')'])
hold on
plot([MeanLims(1) MeanLims(2)],[MeanLims(1) MeanLims(2)],'k-')
hold off
xlabel('Posterior Std.Dev.')
ylabel('VA Std.Dev.')
subplot(3,4,5)
plot(McmcSD,theta1SD,'o','MarkerSize',4.5,'LineWidth',1.3)
title('(e) Normal')
hold on
plot([SDLims(1) SDLims(2)],[SDLims(1) SDLims(2)],'k-')
hold off
xlabel('Posterior Std.Dev.')
ylabel('VA Std.Dev.')
subplot(3,4,6)
plot(McmcSD,theta2SD,'o','MarkerSize',4.5,'LineWidth',1.3)
title('(f) Skew Normal')
hold on
plot([SDLims(1) SDLims(2)],[SDLims(1) SDLims(2)],'k-')
hold off
xlabel('Posterior Std.Dev.')
ylabel('VA Std.Dev.')
subplot(3,4,7)
plot(McmcSD,theta3SD,'o','MarkerSize',4.5,'LineWidth',1.3)
title(['(g) Gaussian Copula (' Transf ')'])
hold on
plot([SDLims(1) SDLims(2)],[SDLims(1) SDLims(2)],'k-')
hold off
xlabel('Posterior Std.Dev.')
ylabel('VA Std.Dev.')
subplot(3,4,8)
plot(McmcSD,theta4SD,'o','MarkerSize',4.5,'LineWidth',1.3)
title(['(h) Skew Normal Copula (' Transf ')'])
hold on
plot([SDLims(1) SDLims(2)],[SDLims(1) SDLims(2)],'k-')
hold off
xlabel('Posterior Std.Dev.')
ylabel('VA Std.Dev.')
subplot(3,4,9)
plot(Skew,theta1Skew,'o','MarkerSize',4.5,'LineWidth',1.3)
title('(i) Normal')
hold on
plot([SKEWLims(1) SKEWLims(2)],[SKEWLims(1) SKEWLims(2)],'k-')
hold off
xlabel('Posterior Skew')
ylabel('VA Skew')
subplot(3,4,10)
plot(Skew,theta2Skew,'o','MarkerSize',4.5,'LineWidth',1.3)
title('(j) Skew Normal')
hold on
plot([SKEWLims(1) SKEWLims(2)],[SKEWLims(1) SKEWLims(2)],'k-')
hold off
xlabel('Posterior Skew')
ylabel('VA Skew')
subplot(3,4,11)
plot(Skew,theta3Skew,'o','MarkerSize',4.5,'LineWidth',1.3)
title(['(k) Gaussian Copula (' Transf ')'])

hold on
plot([SKEWLims(1) SKEWLims(2)],[SKEWLims(1) SKEWLims(2)],'k-')
hold off
xlabel('Posterior Skew')
ylabel('VA Skew')
subplot(3,4,12)
plot(Skew,theta4Skew,'o','MarkerSize',4.5,'LineWidth',1.3)
title(['(l) Skew Normal Copula (' Transf ')'])
hold on
plot([SKEWLims(1) SKEWLims(2)],[SKEWLims(1) SKEWLims(2)],'k-')
hold off
xlabel('Posterior Skew')
ylabel('VA Skew')

set(gcf,'Position',[680,116,1132,862])


%%%Reproducing Table 2 in paper
clear
addpath(genpath('Polypharm/Derivatives/'))
A1 = load('Poly_NoSkewVB_MF.mat');
A2 = load('Poly_NoSkewYJVB_MF.mat');
A3 = load('Poly_NoSkewVB.mat');
A4 = load('Poly_SkewVB.mat');
A5 = load('Poly_NoSkewYJVB.mat');
A6 = load('Poly_SkewYJVB.mat');
A7 = load('Poly_NoSkewiGHVB.mat');
A8 = load('Poly_SkewiGHVB.mat');


LB1 = mean(A1.F_LB((end-1000+1):end));
LB2 = mean(A2.F_LB((end-1000+1):end));
LB3 = mean(A3.F_LB((end-1000+1):end));
LB4 = mean(A4.F_LB((end-1000+1):end));
LB5 = mean(A5.F_LB((end-1000+1):end));
LB6 = mean(A6.F_LB((end-1000+1):end));
LB7 = mean(A7.F_LB((end-1000+1):end));
LB8 = mean(A8.F_LB((end-1000+1):end));

nLambdaA1 = numel(A1.F_mu(:))+numel(A1.F_B(:))-(A3.p)*(A1.p-1)*0.5*(A1.p~=0)+numel(A1.F_d(:))+numel(A1.alphaphi(:))*A1.SkewNorm+numel(A1.eta(:))*~isempty(A1.Transf);
nLambdaA2 = numel(A2.F_mu(:))+numel(A2.F_B(:))-(A5.p)*(A2.p-1)*0.5*(A2.p~=0)+numel(A2.F_d(:))+numel(A2.alphaphi(:))*A2.SkewNorm+numel(A2.eta(:))*~isempty(A2.Transf);
nLambdaA3 = numel(A3.F_mu(:))+numel(A3.F_B(:))-(A3.p)*(A3.p-1)*0.5*(A3.p~=0)+numel(A3.F_d(:))+numel(A3.alphaphi(:))*A3.SkewNorm+numel(A3.eta(:))*~isempty(A3.Transf);
nLambdaA4 = numel(A4.F_mu(:))+numel(A4.F_B(:))-(A4.p)*(A4.p-1)*0.5*(A4.p~=0)+numel(A4.F_d(:))+numel(A4.alphaphi(:))*A4.SkewNorm+numel(A4.eta(:))*~isempty(A4.Transf);
nLambdaA5 = numel(A5.F_mu(:))+numel(A5.F_B(:))-(A5.p)*(A5.p-1)*0.5*(A5.p~=0)+numel(A5.F_d(:))+numel(A5.alphaphi(:))*A5.SkewNorm+numel(A5.eta(:))*~isempty(A5.Transf);
nLambdaA6 = numel(A6.F_mu(:))+numel(A6.F_B(:))-(A6.p)*(A6.p-1)*0.5*(A6.p~=0)+numel(A6.F_d(:))+numel(A6.alphaphi(:))*A6.SkewNorm+numel(A6.eta(:))*~isempty(A6.Transf);
nLambdaA7 = numel(A7.F_mu(:))+numel(A7.F_B(:))-(A7.p)*(A7.p-1)*0.5*(A7.p~=0)+numel(A7.F_d(:))+numel(A7.alphaphi(:))*A7.SkewNorm+numel(A7.eta(:))*~isempty(A7.Transf);
nLambdaA8 = numel(A8.F_mu(:))+numel(A8.F_B(:))-(A8.p)*(A8.p-1)*0.5*(A8.p~=0)+numel(A8.F_d(:))+numel(A8.alphaphi(:))*A8.SkewNorm+numel(A8.eta(:))*~isempty(A8.Transf);

time1 = A1.StoreTime(end)*(1000/A1.niter)/60;
time2 = A2.StoreTime(end)*(1000/A2.niter)/60;
time3 = A3.StoreTime(end)*(1000/A3.niter)/60;
time4 = A4.StoreTime(end)*(1000/A4.niter)/60;
time5 = A5.StoreTime(end)*(1000/A5.niter)/60;
time6 = A6.StoreTime(end)*(1000/A6.niter)/60;
time7 = A7.StoreTime(end)*(1000/A7.niter)/60;
time8 = A8.StoreTime(end)*(1000/A8.niter)/60;


LBall = round([LB1 LB2 LB3 LB4 LB5 LB6 LB7 LB8]'*100)/100;
nLambda = round([nLambdaA1 nLambdaA2 nLambdaA3 nLambdaA4 nLambdaA5 nLambdaA6 nLambdaA7 nLambdaA8]'*100)/100;
Time = round(100*[time1 time2 time3 time4 time5 time6 time7 time8]')/100;
Table2 = [nLambda LBall Time]; %Columns  two and three of table two in the paper

