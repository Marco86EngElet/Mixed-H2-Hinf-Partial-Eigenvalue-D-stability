%% A) Show transient specifications and parameters of the transient regions
disp('Parametros escolhidos para D-estabilidade')
    w_H,
    
    theta_s,
    
    r_d,
    
    q_d,
    
    beta_v,
    
    alpha_v,
    
    e_P,
disp('----------------------------------')    
    
%% B) Table to Compare Poles of the Open-Loop and Closed-Loop System in the Complex Plane

    
    Table_Poles=table(Open_Loop_Poles,...
               Poles_Gcd_classic,...
               Poles_Gcd_Partial);
    Table_Poles,

    
%% C) Make tables to compare norms


    %C.1) Column for H2 norms
    H2_norm=...
        [H2_norm_Godys;H2_norm_Gcdys_classic;H2_norm_Gcdys_partial];
    
    %C.2) Column for discrepancy between H_2 norms of the closed-loop and open-loop system
    H2_norm_percentage_discrepancy=H2_norm-H2_norm_Godys;
    
    H2_norm_percentage_discrepancy=...
        100*H2_norm_percentage_discrepancy/H2_norm_Godys;
    
    %C.3) Column for H_infinity norms
    Hinf_norm=...
        [Hinf_norm_Godzs;Hinf_norm_Gcdzs_classic;Hinf_norm_Gcdzs_partial];
    
    %C.4) Column for discrepancy between H_inf norms of the closed-loop and open-loop system
    Hinf_norm_percentage_discrepancy=Hinf_norm-Hinf_norm_Godzs;
    
    Hinf_norm_percentage_discrepancy=...
        100*Hinf_norm_percentage_discrepancy/Hinf_norm_Godzs;
    
    %C.5) Column to specify control designs

    Project=["Open-Loop";"Classic";"Theorem 7"];
    
    %C.6) Optimization time of each design

    Optimization_time=...
        [NaN('single');classic_optimization_time;partial_optimization_time];
  
	%C.7) Show Table
      Table_Norms=...
        table(Project,H2_norm,H2_norm_percentage_discrepancy,...
              Hinf_norm,Hinf_norm_percentage_discrepancy,...   
              Optimization_time),    
