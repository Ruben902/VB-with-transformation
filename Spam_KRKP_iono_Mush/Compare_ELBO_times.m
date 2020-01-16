%%Plotting Lower Bound vs Time for alternative variational approximations (Figure 7, Table 3 in paper)
clear
Data_sets = {'Spam' 'Krkp' 'Iono' 'Mush'};
Labels = {'Spam' 'Krkp' 'Ionosphere' 'Mushroom'};
letters = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)'};
UB = [0.7 0.5 1 0.7]*0+0.3;
Table3 = zeros(4,6);
for var = 1:4
    vari = Data_sets{var};
    vari2 = Labels{var};
    S =50000;
    A3 = load([vari '_NoSkewk3VB.mat']);
    A4 = load([vari '_Skewk3VB.mat']);
    A5 = load([vari '_NoSkewYJk3VB.mat']);
    A6 = load([vari '_SkewYJk3VB.mat']);
    A7 = load([vari '_NoSkewiGHk3VB.mat']);
    A8 = load([vari '_SkewiGHk3VB.mat']);
    subplot(4,2,(var-1)*2+1)
    nini = 200;
    nend = 15000;
    uni = 60; %seconds
    plot(A3.StoreTime(nini:nend)/uni,A3.StoreLB(nini:nend))
    hold on
    plot(A7.StoreTime(nini:nend)/uni,A6.StoreLB(nini:nend))
    plot(A5.StoreTime(nini:nend)/uni,A5.StoreLB(nini:nend))
    hold off
    ylabel('Lower Bound')
    xlabel('Clock time (Minutes)')
    title([letters{(var-1)*2+1} ' ' vari2])
    xlim([0 UB(var)])
    legend('Gaussian','Gaussian Copula (iGH)','Gaussian Copula (YJ)','Location','SouthEast','FontSize',7)
    legend boxoff
        
    subplot(4,2,var*2)
    plot(A4.StoreTime(nini:nend)/uni,A3.StoreLB(nini:nend))
    hold on
    plot(A8.StoreTime(nini:nend)/uni,A6.StoreLB(nini:nend))
    plot(A6.StoreTime(nini:nend)/uni,A5.StoreLB(nini:nend))
    hold off
    ylabel('Lower Bound')
    xlabel('Clock time (Minutes)')
    title([letters{var*2} ' ' vari2])
    xlim([0 UB(var)])
    Table3(var,:) = round(100*[mean(A3.StoreLB((end-1000+1):end)) mean(A4.StoreLB((end-1000+1):end)) mean(A5.StoreLB((end-1000+1):end)) mean(A6.StoreLB((end-1000+1):end)) mean(A7.StoreLB((end-1000+1):end)) mean(A8.StoreLB((end-1000+1):end))])/100;
    %Table3v2(var,:) = round(100*[max(A3.StoreLB((end-1000+1):end)) max(A4.StoreLB((end-1000+1):end)) max(A5.StoreLB((end-1000+1):end)) max(A6.StoreLB((end-1000+1):end)) max(A7.StoreLB((end-1000+1):end)) max(A8.StoreLB((end-1000+1):end)) ])/100;
    %Table3v3(var,:) = round(100*[median(A3.StoreLB((end-1000+1):end)) median(A5.StoreLB((end-1000+1):end)) median(A4.StoreLB((end-1000+1):end)) median(A6.StoreLB((end-1000+1):end)) median(A7.StoreLB((end-1000+1):end)) median(A8.StoreLB((end-1000+1):end))])/100;

    legend('Skew Normal','Skew Normal Copula (iGH)','Skew Normal Copula (YJ)','Location','SouthEast','FontSize',7)
        legend boxoff
end


