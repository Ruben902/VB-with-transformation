function ddt_dtheta = ddtYJ_dtheta(theta,eta)
dt_dtheta = dtYJ_dtheta(theta,eta);
c2 = theta<0;
c3 = 0<=theta;
Tq1vec(c2) = (eta(c2)-1)./(-theta(c2)+1);
Tq1vec(c3) = (eta(c3)-1)./(theta(c3)+1);
ddt_dtheta = Tq1vec.*dt_dtheta;