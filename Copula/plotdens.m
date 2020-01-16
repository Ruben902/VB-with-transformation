function p = plotdens(obj,varargin)
eps = 0.00001;
yl = icdf(obj,eps);
yu = icdf(obj,1-eps);
ygrid = yl:((yu-yl)/(199)):yu;
p = plot(ygrid,pdf(obj,ygrid),varargin{:});


