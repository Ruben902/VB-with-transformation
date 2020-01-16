function tau = kendalltau(Par,family)
numdim = length(size(Par));
if numdim ==2
    switch family
        case 'gumbel_mix'
            p = size(Par,2)/5;
            tau = zeros(size(Par,1),p);
            for cop = 1:p
                tau(:,cop) = Par(:,(cop-1)*5+1).*Par(:,cop*5)-Par(:,(cop-1)*5+3).*(1-Par(:,cop*5));
            end
        case 't_mix'
            error('Kendalls tau computations needs to be coded')
        case 't'
            error('Kendalls tau computations needs to be coded')
    end
elseif numdim >2
    switch family{1,1,1}
        case 'gumbel_mix'
            [~,r,pp1,numpar,J] = size(Par);
            tau = zeros(r,r,pp1,J);
            for i = 1:r
                for j = 1:r
                    for p = 1:pp1
                        tau(i,j,p,:) = Par(i,j,p,1,:).*Par(i,j,p,5,:)-Par(i,j,p,3,:).*(1-Par(i,j,p,5,:));
                    end
                end
            end
        case 't_mix'
            error('Kendalls tau computations needs to be coded')
        case 't'
            error('Kendalls tau computations needs to be coded')
    end
end