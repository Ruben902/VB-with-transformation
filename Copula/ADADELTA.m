function [Change_delta,ADA] = ADADELTA(ADA,L) 
if isempty(ADA)
    Edelta2 = zeros(length(L),1);
    Eg2 = zeros(length(L),1);
    ADA.rho = 0.95;
    ADA.eps_step = 10^-6;
    ADA.Edelta2 = Edelta2;
    ADA.Eg2 = Eg2;
end
rho = ADA.rho;
eps_step = ADA.eps_step;
oldEdelta2 = ADA.Edelta2;
oldEg2 = ADA.Eg2;
%% mu update
ADA.Eg2 = rho*oldEg2 + (1-rho)*L.^2;
Change_delta = sqrt(oldEdelta2 + eps_step)./sqrt(ADA.Eg2 + eps_step).*L;
ADA.Edelta2 = rho*oldEdelta2 + (1- rho)*Change_delta.^2;