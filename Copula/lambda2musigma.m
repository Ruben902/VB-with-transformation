function [mu,Sigma,B,D,muz,logsigmaz,C,tau,tauu] = lambda2musigma(lambda,d,T,k,VA,Transf)
if VA == 1
    mu = lambda(1:d);
    B =  vechB2B(lambda((d+1):((k+1)*d-(k-1)*k*0.5)),d,k);
    c = ((k+1)*d-(k-1)*k*0.5);
    D =  sparse(1:d,1:d,lambda((c+1):(c+d)),d,d);
    Sigma = B*B'+D.^2;
    muz = [];
    logsigmaz = [];
    C = [];
    c = c+d;
elseif VA == 2
    mu = lambda(1:d);
    B =  vechB2B(lambda((d+1):((k+1)*d-(k-1)*k*0.5)),d,k);
    c = ((k+1)*d-(k-1)*k*0.5);
    D =  sparse(1:d,1:d,lambda((c+1):(c+d)),d,d);
    c = c+d;
    muz = lambda((c+1):(c+T));
    c = c+T;
    logsigmaz = lambda((c+1):(c+T));
    c = c+T;
    Sigma = B*B'+D.^2;
    C = [];
elseif VA == 3
    c = 0;
    mu = lambda((c+1):(c+d));
    c = c+d;
    B =  vechB2B(lambda((c+1):(c+(k*d-0.5*(k-1)*k))),d,k);
    c = c + (k*d-0.5*(k-1)*k);
    D =  sparse(1:d,1:d,lambda((c+1):(c+d)),d,d);
    Sigma = B*B'+D.^2;
    c = c+d;
    muz = lambda((c+1):(c+T));
    c = c+T;
    PositiveC = diag(ones(T,1))+diag(ones(T-1,1),-1);
    indPositiveC = (PositiveC(:)==1);
    C = zeros(T,T);
    C(indPositiveC) = lambda((c+1):(c+2*T-1));
    c = c+2*T-1;
    C = sparse(C);
    logsigmaz = [];
end
switch Transf
    case 'YJ'
        tau = lambda((c+1):(c+d));
        c = c+d;
        tauu = [];
        if VA>1
        tauu= lambda((c+1):(c+T));
        end
    case 'GH'
        tau = lambda((c+1):(c+2*d));
        c = c+d;
        tauu= lambda((c+1):(c+2*T));
    case 'none'
        tau=[];
        tauu =[];
end


