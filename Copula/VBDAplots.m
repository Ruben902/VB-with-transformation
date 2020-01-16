function VBDAplots(VBDAobj)
letters = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)'};
figure
subplot(1,2,1)
plot(1:length(VBDAobj.y),VBDAobj.y)
ylabel('y_t')
title(letters{1})
xlim([0 length(VBDAobj.y)])
subplot(1,2,2)
histogram(VBDAobj.y,'Normalization','pdf')
ylabel('Relative Frequency')
xlabel('y')
title(letters{2})

figure
s = 1;
for k = 1:VBDAobj.p
    subplot(1,VBDAobj.p,k)
    mix_cop_plot(VBDAobj.gamma_mean((k-1)*5+1:(5*k)),[0.02 0.98],1,VBDAobj.family,true);
    title([letters{s} ' ' ': log c_' mat2str(k+1)])
    zlabel(['log(c_' mat2str(k+1) ')'])
    s = s+1;
end

figure
plot(1:VBDAobj.nVB,VBDAobj.LB)
ylabel('Lower Bound')
xlabel('VB steps')
title('Lower Bound')
