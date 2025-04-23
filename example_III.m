clear all;
close all;
clc;

%   1) Criar matriz A_x formato de toeplitz
    A_x=toeplitz([0.1:0.01:0.5,0.01]);
    n_x=size(A_x,1);
    %---------------------------------
%   2)  Criar matriz B_u
    n_u=2;
    B_u=zeros(n_x,n_u);
    B_u(end-1:end,1:2)=eye(2);
    %---------------------------------
%   3) Criar matriz B_d
    n_d=5;
    B_d=0.1*eye(n_x,n_d);
    %---------------------------------
%   4) Criar matriz C_y
    n_y=1;
    C_y=zeros(n_y,n_x);
    C_y(1,end)=1;
    C_y(1,end-1)=-0.5;
    C_y(1,end-2)=0.25;
    %---------------------------------
%   5) Criar matriz D_y
    D_y=zeros(n_y,n_u);
    D_y(1,1)=0.01;
    %--------------------------------    
%   6) Criar matriz E_y
    E_y=zeros(n_y,n_d);
    %--------------------------------
%   7) Criar matriz C_z
    n_z=1;
    C_z=zeros(n_z,n_x);
    C_z(1,1)=1;
    C_z(1,end)=-1;
    %-----------------------------  
%   8) Criar matriz D_z
    D_z=zeros(n_z,n_u);
    D_z(1,1)=0.01;
    %----------------------------
%   9) Criar matriz E_z
    E_z=zeros(n_z,n_d);
    E_z(1,1)=0.01;
    E_z(1,2)=-0.01;
    E_z(1,3)=0.01;
    E_z(1,4)=-0.01;
    E_z(1,5)=0.01;
    %------------------------------

Canon_Godys=ss(A_x,B_d,C_y,0);
Canon_Godzs=ss(A_x,B_d,C_z,E_z);

%5) Input for cost function 

c_H2=1;
c_Hinf=2;    

%6) Input parameter for D-stability

theta_s=0.84045;
r_d=2.8074;
q_d=0;
alpha_v=0.89093;
beta_v=2.4;
    
%7)  Configure optimization settings for YALMIP
Yalmip_sdpsettings =...
     sdpsettings('verbose',1,'solver','sdpt3','debug',1);

%8)  Make Decision Variables for D-stability LMIs
X  = sdpvar(n_x,n_x,'symmetric');
W  = sdpvar(n_u,n_x,'full');
gamma = sdpvar(1,1,'symmetric'); 
Z = sdpvar(n_y,n_y,'symmetric');
rho = sdpvar(1,1,'symmetric');   

%13) Input Q_p
Q = sqrt(2)*eye(2);
n_p=2;
     
% )Make Closed Loop System

Open_Loop_Poles=cplxpair(eig(A_x));

% 14) Start optimization time count for partial allocation
Kpf_theorem_7=zeros(n_u,n_x);
Kpf_theorem_1=zeros(n_u,n_x);

% 15) Counting structure for partial allocation optimization

    % 15.1) Obtain the Eigenvalue and Eigenvector Matrices
    [Right_Eigenvectors, Eigenvalues, Left_Eigenvectors] = eig(A_x) ;

    % 15.7) Find the matrices for partial allocation
    [row_poles,col_poles]=find(real(Eigenvalues)>0);
    
    Lambda_p = Eigenvalues(row_poles, col_poles);
    L_p = Left_Eigenvectors(1:end, col_poles);
    
    % 15.8) Obtain partial controllability matrix
    controlability_matrix = ctrb(Lambda_p, (Q * L_p') * B_u / 2);

    % 15.9) Calculate the number of uncontrollable states in the reduced system
    number_of_non_controlable_states = ...
        2 - rank(controlability_matrix);

    % 15.10) Execute optimization for partial allocation
    if number_of_non_controlable_states == 0

        % 15.10.1) Compute similar transformations     
        Tilde_Lambda_p = Q*Lambda_p*Q'/2;
        Tilde_Bu = Q*L_p'*B_u/2;
        Tilde_Bd = Q*L_p'*B_d/2;
        proj_Lp=L_p/(L_p'*L_p);
        Tilde_Cy = C_y*proj_Lp*Q'/2;
        Tilde_Cz = C_z*proj_Lp*Q'/2;
        
        % 15.10.2) Make Decision Variables 
        
        Tilde_X  = sdpvar(n_p,n_p,'symmetric');
        Tilde_W  = sdpvar(n_u,n_p,'full');
        
        % 15.10.3) Make LMIs for theorem 6
        set_LMIs = Tilde_X>=eps*eye(n_p);

        % 15.10.4) Make LMI for Re(s)<=-alpha_v
            Tilde_Ax_Valpha = Tilde_Lambda_p+alpha_v*eye(n_p);
            set_LMIs = [ set_LMIs,...
                Tilde_Ax_Valpha*Tilde_X+Tilde_X*Tilde_Ax_Valpha'+...
                Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu'<=-eps*eye(n_p)  ];
        

        % 15.10.5) Make LMI for Re(s)>=-beta_v
        
            Tilde_Ax_Vbeta = -Tilde_Lambda_p-beta_v*eye(n_p);
            set_LMIs = [ set_LMIs,...
                Tilde_Ax_Vbeta*Tilde_X+Tilde_X*Tilde_Ax_Vbeta'-...
                Tilde_Bu*Tilde_W-Tilde_W'*Tilde_Bu'<=-eps*eye(n_p) ];
        

        % 15.10.6) Make LMI for abs(s)<=-r_d
        
            Tilde_Ax_Dqr = Tilde_Lambda_p+q_d*eye(n_p); 
            set_LMIs = [ set_LMIs,...
                [ -r_d*Tilde_X,...
                  Tilde_Ax_Dqr*Tilde_X+Tilde_Bu*Tilde_W;...
                  Tilde_X*Tilde_Ax_Dqr'+Tilde_W'*Tilde_Bu',...
                  -r_d*Tilde_X ]...
                  <=-eps*eye(2*n_p)    ];
        

        % 15.10.7) Make LMI for -theta_s<=angle(s)<=theta_s
        
            Tilde_Axsin = Tilde_Lambda_p*sin(theta_s);
            Tilde_Axcos = Tilde_Lambda_p*cos(theta_s);
            Tilde_Busin = Tilde_Bu*sin(theta_s);
            Tilde_Bucos = Tilde_Bu*cos(theta_s);
            set_LMIs =[ set_LMIs,...
                [ Tilde_Axsin*Tilde_X+Tilde_X*Tilde_Axsin'+...
                  Tilde_Busin*Tilde_W+Tilde_W'*Tilde_Busin',...
                  Tilde_Axcos*Tilde_X-Tilde_X*Tilde_Axcos'+...
                  Tilde_Bucos*Tilde_W-Tilde_W'*Tilde_Bucos';...
                  Tilde_X*Tilde_Axcos'-Tilde_Axcos*Tilde_X+...
                  Tilde_W'*Tilde_Bucos'-Tilde_Bucos*Tilde_W,...
                  Tilde_Axsin*Tilde_X+Tilde_X*Tilde_Axsin'+...
                  Tilde_Busin*Tilde_W+Tilde_W'*Tilde_Busin'...
                ]<=-eps*eye(2*n_p) ];
        

        % 15.10.10) Make LMI for H2-Control
           Z = sdpvar(n_y,n_y,'symmetric');
            rho = sdpvar(1,1,'symmetric');   
            
            set_LMIs = [ set_LMIs, ...
                Tilde_Lambda_p*Tilde_X+Tilde_X*Tilde_Lambda_p'+...
                Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu'+...
                Tilde_Bd*Tilde_Bd'<=-eps*eye(n_p) ];
            
            set_LMIs = [ set_LMIs,...
                        Z>=eps*eye(n_y),...
                        trace(Z)<=rho,...
                        rho>=eps ];
                    
            set_LMIs = [set_LMIs,...
                [ -Z, Tilde_Cy*Tilde_X+D_y*Tilde_W;...
                    Tilde_X'*Tilde_Cy'+Tilde_W'*D_y', -Tilde_X ]<=...
                    -eps*eye(n_y+n_p,n_y+n_p) ]        
        
        % 15.10.11) Make LMI for Hinf-Control
          gamma = sdpvar(1,1,'symmetric');   
            set_LMIs = [ set_LMIs,...
                gamma>=eps,...
                [   Tilde_Lambda_p*Tilde_X+Tilde_X*Tilde_Lambda_p'+...
                    Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu',...
                    Tilde_Bd, Tilde_X*Tilde_Cz'+Tilde_W'*D_z';...
                    Tilde_Bd', -gamma*eye(n_d), E_z';...
                    Tilde_Cz*Tilde_X+D_z*Tilde_W,...
                    E_z, -gamma*eye(n_z)...
                ]<=-eps*eye(n_p+n_d+n_z)...
            ];
        

        % 15.10.12) Execute Optimization 
        tic
        optimize(set_LMIs,...
                 [],Yalmip_sdpsettings); 
        time_theorem_1=toc;
        
        % 15.10.13) Extract optimal variables 
        optimal_Tilde_W = value(Tilde_W);
        optimal_Tilde_X = value(Tilde_X);
        
        % 15.10.14) Compute and update feedback matrix
        
        Kpf_parcial=optimal_Tilde_W/(optimal_Tilde_X);
        Kpf_theorem_1 = Kpf_theorem_1 + Kpf_parcial * (Q * L_p') / 2;     
                
        % 16.10.12) Execute Optimization   
        tic
        optimize(set_LMIs,...
                 c_H2*rho+c_Hinf*gamma,Yalmip_sdpsettings); 
        time_theorem_7=toc;
        
        % 15.10.13) Extract optimal variables 
        optimal_Tilde_W = value(Tilde_W);
        optimal_Tilde_X = value(Tilde_X);
        
        % 15.10.14) Compute and update feedback matrix
        Kpf_parcial=optimal_Tilde_W/(optimal_Tilde_X);
        Kpf_theorem_7 = Kpf_theorem_7 + Kpf_parcial * (Q * L_p') / 2;     
                 optimize(set_LMIs,[],Yalmip_sdpsettings);
        

    else
        disp('uncontrollable pole detected')
    end

% 15.20) End multi-step optimization time


% 16) Make closed loop system d(t) to y(t) via theorem 2 
system_y_theorem_1=...
ss(A_x+B_u*Kpf_theorem_1,B_d,C_y+D_y*Kpf_theorem_1,E_y);

% 17) Make closed loop system d(t) to z(t) via theorem 2 
system_z_theorem_1=...
ss(A_x+B_u*Kpf_theorem_1,B_d,C_z+D_z*Kpf_theorem_1,E_z); 

% 18) Make closed loop system d(t) to y(t) via multistep theorem 6
system_y_theorem_7=...
ss(A_x+B_u*Kpf_theorem_7,B_d,C_y+D_y*Kpf_theorem_7,E_y);

% 19) Make closed loop system d(t) to z(t) via multistep theorem 6
system_z_theorem_7=...
ss(A_x+B_u*Kpf_theorem_7,B_d,C_z+D_z*Kpf_theorem_7,E_z); 

% 20) Closed-Loop system poles via theorem 2
Poles_system_theorem_1=cplxpair(pole(system_y_theorem_1));
    
% 21) Closed-Loop system poles via theorem 6   
Poles_system_theorem_7=cplxpair(pole(system_y_theorem_7));
       
    
    Table_Poles=table(Open_Loop_Poles,...
               Poles_system_theorem_1,...
               Poles_system_theorem_7);
           
% 22) Display table                   
    Table_Poles,
    
% 23) Compute open-loop system norms
%       23.1) H2 norm for Gdy in open loop
    norm_H2(1,1)=norm(Canon_Godys,2);
    
%       23.2) variation of H2 norm for Gdy in open loop            
    discrepance_norm_H2(1,1)=0;
    
%       23.3) H-infinity norm for Gdy in open loop            
    norm_Hinf(1,1)=norm(Canon_Godzs,'inf');
    
%       23.4) variation of H-infinity norm for Gdinf in open loop            
    discrepance_norm_Hinf(1,1)=0;
    
%       23.5) Cost function in open loop            
    Cost_Function(1,1)=c_H2*norm_H2(1,1)+c_Hinf*norm_Hinf(1,1);
    
%       23.6) variation of cost function for Gdinf in open loop            
    Discrepance_cost_function(1,1)=0;
    
%   24) Compute closed-loop norms for free allocation
%       24.1) H2 norm for Gdy in closed loop
    norm_H2(2,1)=norm(system_y_theorem_1,2);

%       24.2) variation of H2 norm for Gdy in closed loop            
discrepance_norm_H2(2,1)=...
    100*(norm_H2(2,1)-norm_H2(1,1))/norm_H2(1,1);

%       24.3) H-infinity norm for Gdy in closed loop            
norm_Hinf(2,1)=norm(system_z_theorem_1,'inf');

%       24.4) variation of H-infinity norm for Gdinf in closed loop            
discrepance_norm_Hinf(2,1)=...
    100*(norm_Hinf(2,1)-norm_Hinf(1,1))/norm_Hinf(1,1);

%       24.5) Cost function for free allocation            
    Cost_Function(2,1)=c_H2*norm_H2(2,1)+c_Hinf*norm_Hinf(2,1);
    
%       24.6) variation of cost function for free allocation            
    Discrepance_cost_function(2,1)=0;
    
%   25) Compute closed-loop norms for partial allocation
%       25.1) H2 norm for Gdy in closed loop
norm_H2(3,1)=norm(system_y_theorem_7,2);

%       25.2) variation of H2 norm for Gdy in closed loop            
discrepance_norm_H2(3,1)=...
    100*(norm_H2(3,1)-norm_H2(1,1))/norm_H2(1,1);

%       25.3) H-infinity norm for Gdy in closed loop            
norm_Hinf(3,1)=norm(system_z_theorem_7,'inf');

%       25.4) variation of H-infinity norm for Gdinf in closed loop            
discrepance_norm_Hinf(3,1)=...
    100*(norm_Hinf(3,1)-norm_Hinf(1,1))/norm_Hinf(1,1);

%       25.5) Cost function for partial allocation            
Cost_Function(3,1)=c_H2*norm_H2(3,1)+c_Hinf*norm_Hinf(3,1);
    
%       25.6) variation of cost function for partial allocation            
Discrepance_cost_function(3,1)=0;
    
%   26) Compute Optimization Time
%       26.1) Optimization time for open loop
Optimization_Time(1,1)=0;
    
%       26.2) Optimization time for theorem 2
Optimization_Time(2,1)=...
    time_theorem_1; 
    
%       26.3) Optimization time for multistep theorem 6
Optimization_Time(3,1)=...
    time_theorem_7; 

%   27) Table to compare norms
Transfer_Functions=["Open Loop";"Theorem 2";...
                    "Multistep Theorem 6"];
Table_Norms = table(Transfer_Functions,...
norm_H2,discrepance_norm_H2,...
norm_Hinf,discrepance_norm_Hinf,...
Optimization_Time),

%   28) Table - Optimization Time
percentual_discrepance_optimization_time=...
    100*(time_theorem_7-time_theorem_1)/time_theorem_1;

Table_Time_Optimization=table(time_theorem_1,...
    time_theorem_7,...
    percentual_discrepance_optimization_time),


% 29) Plot for Poles and Position

%   29.1) Eigenvalues (Poles) to be inserted in the plot
eigenvalues = [ Open_Loop_Poles',...
               Poles_system_theorem_1',...
               Poles_system_theorem_7'];

%   29.2) Computing boundary curves of the regions
% Inicializando arrays        
x_horizontal_strip=[];
x_sector=[] ;
y_inferior_sector=[] ;
y_superior_sector=[] ;
x_disk=[] ;
y_disk=[] ;
x_alphav=[] ;
y_alphav=[] ;
x_betav=[] ;
y_betav=[] ;            
            
% 29.3) Definindo maior e menor parte real dos polos do sistema
min_real_aut=min(real(eigenvalues)); 
max_real_aut=max(real(eigenvalues));                

% 29.4) Computar maior e menor parte real da regiao disco
min_real_disk=-q_d-r_d;
max_real_disk=-q_d+r_d;

% 29.5) Computar maior e menor parte real da regiao dos pontos coordenados
min_real=min([min_real_aut,min_real_disk,-beta_v]);

max_real=max([max_real_aut,max_real_disk,-alpha_v]);
max_real=min([max_real,0]);

% 29.6) Computar maior parte imaginaria dos polos do sistema
max_imag_aut=max(imag(eigenvalues));

% 29.7) Computar maior parte imaginaria da regiao setor
max_imag_sector=-min_real*tan(theta_s);

% 29.8) Computar maior parte imaginaria dos pontos coordenados
max_imag=max([max_imag_aut, max_imag_sector, r_d]);

% 29.9) Construir pontos coordenados das curvas de fronteira
for i=1:100
    xt=min_real+(i-1)*(max_real-min_real)/99;
    
    if ~isempty(theta_s)
        x_sector(i,1)=xt;
        y_superior_sector(i,1)=xt*tan(theta_s);
        y_inferior_sector(i,1)=xt*tan(-theta_s);
    end
    
    if ~isempty(r_d)
        if isempty(q_d)
            x_disk(i,1) = r_d*cos(2*pi*i/100);
            y_disk(i,1) = r_d*sin(2*pi*i/100);
        else
            x_disk(i,1) = -q_d+r_d*cos(2*pi*i/100);
            y_disk(i,1) = -q_d+r_d*sin(2*pi*i/100);
        end
    end
    
    if ~isempty(alpha_v)
        x_alphav(i,1) = -alpha_v;
        y_alphav(i,1) = -max_imag+2*max_imag*(i-1)/99;
    end
    
    if ~isempty(beta_v)
        x_betav(i,1) = -beta_v;
        y_betav(i,1) = -max_imag+2*max_imag*(i-1)/99;
    end
end

%   29.10) Graphs           
    figure
    subplot(211)
    plot(real(Open_Loop_Poles),imag(Open_Loop_Poles),'ok',...
         real(Poles_system_theorem_1),imag(Poles_system_theorem_1),'xr',...
         real(Poles_system_theorem_7),imag(Poles_system_theorem_7),'xb',...
         x_alphav,y_alphav,'-c',...
         x_betav,y_betav,'-c',...
         x_sector,y_inferior_sector,'-c',...
         x_sector,y_superior_sector,'-c',...
         x_disk,y_disk,'-c',...
         'LineWidth',2,'MarkerSize',8);

title('Complex Plane','FontSize',12)
legend('Open','Theorem 2','Theorem 6','Frontier','FontSize',10)
    subplot(212)
    plot(real(Open_Loop_Poles),imag(Open_Loop_Poles),'ok',...
         real(Poles_system_theorem_1),imag(Poles_system_theorem_1),'xr',...
         real(Poles_system_theorem_7),imag(Poles_system_theorem_7),'xb',...
         x_alphav,y_alphav,'-c',...
         x_betav,y_betav,'-c',...
         x_sector,y_inferior_sector,'-c',...
         x_sector,y_superior_sector,'-c',...
         x_disk,y_disk,'-c',...
         'LineWidth',2,'MarkerSize',8);
xlim([-3 -0.35]);
ylim([-3 3]);
title('Complex Plane','FontSize',12)
legend('Open','Theorem 1','Theorem 7','Frontier','FontSize',10)

%   30) Numerical Integral Table of the Absolute Value of the Control Action
%      30.1) Transient Response - free allocation
[y_theorem_1,t_theorem_1,x_theorem_1] = ...
    step(system_y_theorem_1,0:0.01:14);

%      30.2) Control action - free allocation         
u_theorem_1(1:size(x_theorem_1,1),1:n_u,1)=...
    x_theorem_1(:,:,1)*Kpf_theorem_1';

%      30.3) Transient Response - partial multistep allocation                                 
[y_theorem_7,t_theorem_7,x_theorem_7] = ...
    step(system_y_theorem_7,0:0.01:14);

%      30.4) Control action - partial multistep allocation 
u_theorem_7(1:size(x_theorem_1,1),1:n_u,1)=...
    x_theorem_7(:,:,1)*Kpf_theorem_7';

%   31) Figure - control actions
figure
%       31.1) Response u1(t) to step input d1(t)    
subplot(211)
plot(t_theorem_1,u_theorem_1(1:end,1),'-r',...
     t_theorem_7,u_theorem_7(1:end,1),'-b')
title('u_1(t) before step d(t)','FontSize',12)
legend( 'theorem 1','theorem 7','FontSize',10)

%       31.2) Response u2(t) to step input d1(t)            
subplot(212)
plot(t_theorem_1,u_theorem_1(1:end,2),'-r',...
     t_theorem_7,u_theorem_7(1:end,2),'-b')
title('u_2(t) before step d(t)','FontSize',12)
legend('theorem 1','theorem 7','FontSize',10)
    