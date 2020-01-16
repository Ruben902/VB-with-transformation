clear all;

for data_pick = 1:4

    
switch data_pick 
    case 1
        disp('Running Spam data...')
        load('Data/spamdata.mat')
    case 2
        disp('KRKP data')
        load('Data/krkpdata.mat')
    case 3
        disp('ionosphere data')
        load('Data/ionospheredata.mat')
    case 4
        disp('Mushroom data')
        load('Data/mushroomdata.mat')
    otherwise
        disp('Please pick one of the three datasets')
end


spam_train = fulldata;

Xtr = spam_train(:,2:end);
Xtr = [ones(size(Xtr,1),1) spam_train(:,2:end)];
Ytr = 2*spam_train(:,1) - 1;


q = size(Xtr,2);
tic
MCMCsample = logist2SampleMH(Xtr,Ytr,inv(diag(10*ones(q,1))),500001);
toc
index = floor(linspace(400000,500000,5000));
MCMCmean = mean(MCMCsample(:,index),2);
MCMCsd = sqrt(var(MCMCsample(:,index)'));
MCMCskew = skewness(MCMCsample(:,index)');

MCMCsample2 = MCMCsample(:,10:10:end);

switch data_pick 
    case 1
        disp('Running Spam data...')
        save('Output/MCMC_spam.mat','MCMCmean','MCMCsd')
    case 2
        disp('KRKP data')
        save('Output/MCMC_KRKP.mat','MCMCmean','MCMCsd')
    case 3
        disp('ionosphere data')
        save('Output/MCMC_iono.mat','MCMCmean','MCMCsd')
    case 4
        disp('Mushroom data')
        save('MCMC_mush.mat','MCMCmean','MCMCsd','MCMCskew','MCMCsample2')
    otherwise
        disp('Please pick one of the three datasets')
end
end
