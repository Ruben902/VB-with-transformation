function u = uDaibi2u(isDiscrete,ai,bi,uD)
[T,r] = size(ai);
S = size(uD,2);
u = zeros(r*T,S);
discreteElements = zeros(T,r);
discreteElements(:,isDiscrete==1)=1;
discreteElementst = discreteElements';

bitc = bi(:,isDiscrete==0)';
bitc = bitc(:);

u(discreteElementst(:)==1,:)= uD;
u(discreteElementst(:)==0,:)= repmat(bitc,1,S);