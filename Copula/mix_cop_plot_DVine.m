function mix_cop_plot_DVine(Param,m,p,u_lim,scale,family,Log)
ParamM = Gamma2Param(Param,m,p);
numplots = length(Param)/5;
cc= 25;
v = linspace(u_lim(1),u_lim(2),cc);
[U1,U2] = meshgrid(v,v);
if nargin ==4
    Log = false;
end
letters = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' '(k)' '(l)' '(m)' '(n)' '(o)' '(p)' '(q)' '()' '()' '()' '()' '()' '()' '()' '()' '()' };
if strcmp(family,'t_mix')
    %F = w*copulapdf('t',[U1(:) U2(:)],Par(1,i),Par(2,i))+(1-w)*copulapdf('t',[1-U1(:) U2(:)],Par(3,i),Par(4,i));
    FF = @(U1,U2,w,Par) w*copulapdf('t',[U1(:) U2(:)],Par(1),Par(2))+(1-w)*copulapdf('t',[1-U1(:) U2(:)],Par(3),Par(4));
elseif strcmp(family,'sjc_mix')
    %F = exp(lnc_sjc_mix2(U1(:),U2(:),[Par(:,i) w]));
    FF = @(U1,U2,w,Par) exp(lnc_sjc_mix2(U1(:),U2(:),[Par(:) w]));
elseif strcmp(family,'gumbel_mix')
    %F = exp(lnc_gen_gum_mix([U1(:) U2(:)],[Par(1:4,i)' w]));
    FF = @(U1,U2,w,Par) exp(lnc_gen_gum_mix([U1(:) U2(:)],[Par(1:4)' w]));
elseif strcmp(family,'clayton_mix')
    %F = exp(lnc_gen_clay_mix([U1(:) U2(:)],[Par(1:4,i)' w]));
    FF = @(U1,U2,w,Par) exp(lnc_gen_clay_mix([U1(:) U2(:)],[Par(1:4)' w]));
elseif strcmp(family,'skewt_mix')
    
end
Zlabel = 'c(u_{t-1},u_t)';
Tite = ' c';
if Log == true
    Zlabel = 'log(c(u_{t-1},u_t))';
    Tite = ' log c';
end

for j =1:m
    for i =(j+1):m
        %subplot(ceil(numplots/ceil(sqrt(numplots))),ceil(sqrt(numplots)),s)
        pos = subplotij(m,m*(p+1)-1,j,i-1);
        Par = squeeze(ParamM(i,j,1,:));
        w = Par(5);
        F = FF(U1,U2,w,Par);
        Z = F;
        if Log == true
            Z = log(F);
        end
        ZZ = scale*reshape(Z,cc,cc);
        surf(U1,U2,ZZ);
        xlabel('u_{t-1}')
        ylabel('u_t')
        %zlabel(Zlabel)
        title([letters{pos-j*(j-1)*0.5} Tite '_{' mat2str(i) ',' mat2str(j) '}^{(' mat2str(0) ')}'])
        %axis square
    end
end

for t = 1:p
    for i =1:m
        for j =1:m
            %subplot(ceil(numplots/ceil(sqrt(numplots))),ceil(sqrt(numplots)),s)
            pos = subplotij(m,m*(p+1)-1,j,i+t*m-1);
            Par = squeeze(ParamM(i,j,t+1,:));
            w = Par(5);
            F = FF(U1,U2,w,Par);
            Z = F;
            if Log == true
                Z = log(F);
            end
            ZZ = scale*reshape(Z,cc,cc);
            surf(U1,U2,ZZ);
            xlabel('u_{t-1}')
            ylabel('u_t')
            %zlabel(Zlabel)
            title([letters{pos-j*(j-1)*0.5} Tite '_{' mat2str(i) ',' mat2str(j) '}^{(' mat2str(t) ')}'])
            %axis square
        end
    end
end