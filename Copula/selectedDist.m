function mdist = selectedDist(y,type)
if nargin ==1
    type = 'discrete';
end
switch type
    case 'discrete'                         
        [Fx,x] = ecdf([min(y(:))-1; y(:);max(y(:))+1]);             % Empirical distribution function
        mdist = makedist('PiecewiseLinear','x',[x(1)-1 x(2:end)'],'Fx',Fx(:)');
        %[Fx,x] = ecdf(y(:));             % Empirical distribution function
        %llim=3.d-03;ulim=1.d0-llim;       % Limits on the u's matching those of the MCMC
        %Fx(1)=llim;
        %Fx(end)=ulim;
        %mdist = makedist('PiecewiseLinear','x',[x(1)-2 x(1)-1 x(2:end)' x(end)+1],'Fx',[0 Fx(:)' 1]);
    case 'continous'
        mdist = fit_ssvkernel(y);           % ssv kernel
    case 'continuous'
        mdist = fit_ssvkernel(y);           % ssv kernel
end