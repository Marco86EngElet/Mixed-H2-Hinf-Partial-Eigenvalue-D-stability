%% A) Compute the H_2 and H-infinity norms of the open-loop system   

    %A.1) H_2 norm of the controllable canonical system of Gody(s)
    
    H2_norm_Godys=norm(Canon_Godys,2);
    
    %A.2) H-infinity norm of the controllable canonical system of Godz(s)

    Hinf_norm_Godzs=norm(Canon_Godzs,'inf');
    
%% B) Compute the poles of the open-loop system
    
    Open_Loop_Poles=cplxpair(pole(Canon_Godys));
