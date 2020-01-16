function U = U_transf(U)
    llim=1.d-7;ulim=1.d0-llim;                                
    U(U<llim) =  llim;
    U(U>ulim) =  ulim;
    
