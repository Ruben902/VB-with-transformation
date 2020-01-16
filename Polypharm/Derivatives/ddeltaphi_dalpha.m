function ddeltaphidalpha = ddeltaphi_dalpha(Omega,alpha_phi)
Omegaalpha = Omega*alpha_phi;
const1 = (1+alpha_phi'*Omegaalpha);
ddeltaphidalpha = (1./(const1.^(3/2))).*(const1*Omega-Omegaalpha*Omegaalpha');
