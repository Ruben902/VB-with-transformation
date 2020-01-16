function delta_phi = alphaphi2deltaphi(alpha_phi,Omega)
Omegaalpha = Omega*alpha_phi;
delta_phi = Omegaalpha./sqrt(1+alpha_phi'*Omegaalpha);