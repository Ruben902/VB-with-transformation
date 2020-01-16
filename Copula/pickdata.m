function [y,type] = pickdata(example)
type = [];
switch example
    case 'Murder'
        Data = load('NSWmonthlyOffencesProcessed.mat');
        y = Data.AggregateCrimes(:,1);
    case 'AttemptedMurder'
        Data = load('NSWmonthlyOffencesProcessed.mat');
        y = Data.AggregateCrimes(:,2);
    case 'Manslaughter'
        Data = load('NSWmonthlyOffencesProcessed.mat');
        y = Data.AggregateCrimes(:,4);
    case 'Logit'
        Data = load('wd_logit_simul.mat');
        y = Data.Y;
    case 'Homicide'
        Data = load('NSWmonthlyOffencesProcessed.mat','AggregateCrimes');
        y = Data.AggregateCrimes(:,[1 2 4]);
        type = {'discrete','discrete','discrete'};
    case 'BankVix'
        Data = load('Bankruptcy.mat');
        YdataB = Data.BankruptMonth(:);
        YdataV=csvread('VIX.csv');
        y = [YdataB((end-length(YdataV)+1):end) YdataV(:)];
        type = {'discrete','continuous'};
    case 'Bankruptcy'
        Data = load('Bankruptcy.mat');
        y = Data.BankruptMonth(:);
end
