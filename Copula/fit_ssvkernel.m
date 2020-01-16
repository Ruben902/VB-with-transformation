function [obj1,E_points] = fit_ssvkernel(y,Npoints)
% Fit KDE with adaptive bandwith
if nargin == 1
   Npoints = 500; 
end
a=std(y);
x_max = max(y);
x_min = min(y);
gap = (x_max-x_min)/(Npoints-1);
x_middle = x_min:gap:x_max;                 % Guarantees 500 points between min and max
x_low = (x_min-gap):-gap:(x_min-a);         % Guarantees equal spacing from x_min-a to x_min
x_up = (x_max+gap):gap:(x_max+a);           % Guarantees equal spacing from x_max to x_max + a
if gap<a
    x = [fliplr(x_low) x_middle x_up];
else
    x = [(x_min-gap) x_middle (x_max+gap)];
end
E_points = length(x);    %  Effective number of points
[fx,x,~] = ssvkernel(y,x);
Fx = cumtrapz(x,fx);
if abs((1-Fx(end)))<0.0001
Fx = Fx/Fx(end);  
else
   error('The tolerance level is violated') 
end
obj1=makedist('PiecewiseLinear','x',x,'Fx',Fx);