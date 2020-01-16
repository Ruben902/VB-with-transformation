close all
clear
load('polypharm_INLA.mat')
newstr = cell(9,1);
newstr{1} = PosteriorsFixed.Int;
newstr{2} = PosteriorsFixed.gender;
newstr{3} = PosteriorsFixed.race;
newstr{4} = PosteriorsFixed.age;
newstr{5} = PosteriorsFixed.M1;
newstr{6} = PosteriorsFixed.M2;
newstr{7} = PosteriorsFixed.M3;
newstr{8} = PosteriorsFixed.IM;
newstr{9} = PosteriorsXi;
uiopen('YJ_PAR_densities_plus_MF.fig',1);


obj = gcf;

for i = 1:9
    ax = subplot(3,3,i);
    Xlim = xlim;
    Ylim = ylim;
    hold on
    plot(newstr{i}(:,1),newstr{i}(:,2),':b','LineWidth',2)
    hold off
    xlim(Xlim)
    ylim(Ylim)
end

legend('MCMC','Normal','Copula','Normal-MF','INLA')