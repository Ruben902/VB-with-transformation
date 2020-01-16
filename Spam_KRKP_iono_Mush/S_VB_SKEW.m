clear
clear dtheta_dBDelta;   %Clears persistent variables in the function
rng(1990)
workingDir = pwd;
addpath(genpath([workingDir '/Derivatives']));
%Select data
%Defining the type of approximation
Labs = {'Spam' 'Krkp' 'Iono' 'Mush'};
Skewlab = {'NoSkew' 'Skew'};                                               
SkewNorm = 0;                                                              % Set to zero if skew-normal will not be used. Set to one if skew-normal is to be used
Transf = 'YJ';  
niter = 40000;
p = 3; %% Number of factors
%Selecting the data set
data_pick = 2;  
%%%%%%%%%%%%Delete from here
obj = load([Labs{data_pick} '_' Skewlab{SkewNorm+1} Transf 'k' mat2str(p) 'VB.mat']);
plot(obj.StoreLB)
xlim([1 1000])
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ytr,Xtr,q] =  select_data(data_pick);
%------ Defining log-posteriors and derivatives of the model---------------
logpost = @(theta) log_logreg(theta, Xtr,Ytr,(Ytr+1)/2,false) +  log_normal(theta, zeros(q,1),10*ones(q,1),false);
dlogpost = @(theta) log_logreg(theta, Xtr,Ytr,(Ytr+1)/2,true) +  log_normal(theta, zeros(q,1),10*ones(q,1),true);
%------ Computing approximation -------------------------------------------
[mu,d,B,eta,alphaphi,StoreLB,StoreTime] = VBtransf(logpost,dlogpost,q,p,Transf,SkewNorm,niter,100,false);
%------ Saving results -------------------------------------------
save([Labs{data_pick} '_' Skewlab{SkewNorm+1} Transf 'k' mat2str(p) 'VB.mat'])
%------ Visualising approximate posteriors --------------------------------
ind = [1 10 20 30]; %Indexes of approximate posteriors to plot                                      
plot_posterior(ind,eta,mu,B,d,alphaphi,Transf)
