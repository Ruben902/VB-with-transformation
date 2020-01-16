function [Ytr,Xtr,q] =  select_data(data_pick)
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
%% All data
Xtr = fulldata(:,2:end);
Xtr = [ones(size(Xtr,1),1) fulldata(:,2:end)];
Ytr = 2*fulldata(:,1) - 1;
q = size(Xtr,2);