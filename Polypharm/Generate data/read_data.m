clear all;
fulldata = xlsread('polypharm.xls');


gender = fulldata(:,11);
age = fulldata(:,14);
age = log(age/10);
%age = -1+2*(age-min(age))/(max(age)-min(age)) 
%age = (age-min(age))/(max(age)-min(age))
race = fulldata(:,12);
race(race == 2) = 1;
MHV1 = (fulldata(:,3) == 1);
MHV2 = (fulldata(:,3) == 2);
MHV3 = (fulldata(:,3) == 3);
INPTMHV = fulldata(:,4);
INPTMHV(INPTMHV == 2) = 1;

intercept = kron(eye(500),ones(7,1));

Xtr = [ones(3500,1) gender race age MHV1 MHV2 MHV3 INPTMHV intercept];



Ytr = fulldata(:,2)*2 - 1;

save('polypharm.mat','Ytr','Xtr');
