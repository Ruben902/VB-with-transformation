function log_theta = logpri(theta,family,eta,Transf)
if strcmp(family,'gumbel_mix')
    log_theta = sum(log(0.99)-theta-2*log(exp(-theta)+1),2)';        %log prior on theta, implied by uniform on gamma in line three
end
if strcmp(family,'t')

end
