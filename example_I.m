clear all;
close all;
clc;

%1) Damping coefficient, Natural Frequency, and Gain of the partial functions
    %1.1) Input Damping coefficient
    zeta_c(1,1) = 0.1;  % Function 1
    zeta_c(1,2) = 0.12; % Function 2
    zeta_c(1,3) = 0.15; % Function 3
    zeta_c(1,4) = 0.2;  % Function 4
    zeta_c(1,5) = 0.3;  % Function 5

    %1.2) Input Natural frequency
    w_c(1,1) = 0.1;  % Function 1
    w_c(1,2) = 0.12; % Function 2
    w_c(1,3) = 0.15; % Function 3
    w_c(1,4) = 0.2;  % Function 4
    w_c(1,5) = 0.3;  % Function 5

    %1.3) Gain
    KYc(1:5)=ones(1,5);
    KZc(1:5)=ones(1,5);
    KZc(1,3)=-1;

%2) Make G_c(s), G_ody(s) and G_odz(s)  

    Godys=tf(0,1); %Make G_ody(s)
    Godzs=tf(0.01,1); %Make G_odz(s)
    for c=1:5 %loop for add G_c(s) in G_ody(s) and G_odz(s)
        num = w_c(1,c)^2; % Numerator G_c(s)
        den = [1, 2*zeta_c(1,c)*w_c(1,c), w_c(1,c)^2]; % Denominator G_c(s) 
        G_c = tf(num, den); % Make G_c(s)
        Godys=Godys+KYc(1,c)*G_c; % Add G_c(s) in G_ody(s)
        Godzs=Godzs+KZc(1,c)*G_c; % Add G_c(s) in G_odz(s)
    end
    
%3) Compute Matrices on LCTI-MIMO-PSTF

    %3.1) Compute canonical controlable representation of G_ody(s) 
    Canon_Godys = ss(Godys);
    
    %3.2) Compute canonical controlable representation of G_odz(s)
    Canon_Godzs = ss(Godzs);
    
    %3.3) Compute matrices of open-loop system configuration
    A_x=Canon_Godys.A;
    B_d=Canon_Godys.B;
    C_y=Canon_Godys.C;
    E_y=Canon_Godys.D;
    C_z=Canon_Godzs.C;
    E_z=Canon_Godzs.D;
    
    %3.4) Vector's dimensions

    n_x=size(A_x,1);
    n_y=size(C_y,1);
    n_z=size(C_z,1);
    n_d=size(B_d,2);    
    n_u=2;
    n_ola=n_u;

    %3.6) Compute Matrices of closed-loop system configuration 
    
    B_u=zeros(n_x,n_u);
    B_u(2,1)=1;
    B_u(3,1)=1;
    B_u(2,2)=-1;
    B_u(3,2)=-1;
    
    D_y=0.001*ones(n_y,n_u);
    D_z=D_y;

%5) Input for cost function 

c_H2=1;
c_Hinf=2;    

%6) Input parameter for D-stability

theta_s=1.0974;
w_H=1.2566;
r_d=0.6999;
q_d=0;
alpha_v=0.45799;
beta_v=0.8;
    
%7)  Configure optimization settings for YALMIP
Yalmip_sdpsettings =...
     sdpsettings('verbose',1,'solver','lmilab','debug',1);

%8)  Make Decision Variables for D-stability LMIs
X  = sdpvar(n_x,n_x,'symmetric');
W  = sdpvar(n_u,n_x,'full');
gamma = sdpvar(1,1,'symmetric'); 
Z = sdpvar(n_y,n_y,'symmetric');
rho = sdpvar(1,1,'symmetric');   
        
%9) Make LMIs for D-stability
set_LMIs_theorem_2 = X>=eps*eye(n_x);

    %9.1) LMI for Re(s)<=-alpha_v
    if ~isempty(alpha_v)
        Ax_Valpha = A_x+alpha_v*eye(n_x);
        set_LMIs_theorem_2 =...
            [ set_LMIs_theorem_2, Ax_Valpha*X+X*Ax_Valpha'+...
                 B_u*W+W'*B_u'<=-eps*eye(n_x)  ];
    end

    %9.2) LMI for Re(s)=>-beta_v
    if ~isempty(beta_v)
        Ax_Vbeta = -A_x-beta_v*eye(n_x);
        set_LMIs_theorem_2 = ...
            [ set_LMIs_theorem_2, Ax_Vbeta*X+X*Ax_Vbeta'-...
                 B_u*W-W'*B_u'<=-eps*eye(n_x) ];
    end

    %9.3) LMI for ABS(s)<=r_d
    if ~isempty(r_d)
        Ax_Dqr = A_x+q_d*eye(n_x); 
        set_LMIs_theorem_2 = [ set_LMIs_theorem_2,...
                [ -r_d*X, Ax_Dqr*X+B_u*W;...
                  X*Ax_Dqr'+W'*B_u', -r_d*X ]...
                  <=-eps*eye(2*n_x)    ];
    end

    %9.4) LMI for -w_H<=Imag(s)<=w_H
    if ~isempty(w_H)
        set_LMIs_theorem_2 =[ set_LMIs_theorem_2,...
                [ -w_H*X,...
                  X*A_x'-A_x*X+...
                  W'*B_u'-B_u*W;...
                  A_x*X-X*A_x'+...
                  B_u*W-W'*B_u',...
                  -w_H*X...
                ]<=-eps*eye(2*n_x) ];
    end

    %9.5) LMI for -theta_s<=angle(s)<=theta_s
    if ~isempty(theta_s)
        Ax_sin = A_x*sin(theta_s);
        Ax_cos = A_x*cos(theta_s);
        Bu_sin = B_u*sin(theta_s);
        Bu_cos = B_u*cos(theta_s);
        set_LMIs_theorem_2 =[ set_LMIs_theorem_2,...
                [ Ax_sin*X+X*Ax_sin'+...
                  Bu_sin*W+W'*Bu_sin',...
                  Ax_cos*X-X*Ax_cos'+...
                  Bu_cos*W-W'*Bu_cos';...
                  X*Ax_cos'-Ax_cos*X+...
                  W'*Bu_cos'-Bu_cos*W,...
                  Ax_sin*X+X*Ax_sin'+...
                  Bu_sin*W+W'*Bu_sin'...
                ]<=-eps*eye(2*n_x) ];
    end

    %9.6) LMI for H2-Control
    if ~isempty(c_H2)
        
        set_LMIs_theorem_2 = [ set_LMIs_theorem_2, A_x*X+X*A_x'+...
                 B_u*W+W'*B_u'+...
                 B_d*B_d'<=-eps*eye(n_x) ];
        set_LMIs_theorem_2 = [ set_LMIs_theorem_2,...
                    Z>=eps*eye(n_y),...
                    trace(Z)<=rho,...
                    rho>=eps ];
        set_LMIs_theorem_2 = [set_LMIs_theorem_2,...
            [ -Z, C_y*X+D_y*W; X'*C_y'+W'*D_y', -X ]<=...
                -eps*eye(n_y+n_x,n_y+n_x) ];        
    end

    %9.7) LMIs for Hinfinity-Control
    if ~isempty(c_Hinf)
          
        set_LMIs_theorem_2 = [ set_LMIs_theorem_2,...
                    gamma>=eps,...
                    [   A_x*X+X*A_x'+B_u*W+W'*B_u',...
                        B_d, X*C_z'+W'*D_y';...
                        B_d', -gamma*eye(n_d), E_z';...
                        C_z*X+D_z*W,...
                        E_z, -gamma*eye(n_z)...
                    ]<=-eps*eye(n_x+n_d+n_z)...
                ];
    end

%10) Do optimization
tic
optimize(set_LMIs_theorem_2,...
            c_H2*rho+c_Hinf*gamma,...
            Yalmip_sdpsettings);
time_theorem_2=toc;

%11) Extract optimal solution variables for Theorem_2
optimal_W = value(W);
optimal_X = value(X);

%12) Compute feedback matrix via classical mixed H2/Hinf control
%D-stability
Kpf_theorem_2=optimal_W/optimal_X;   

%13) Input Q_p
Q = [1 1;1i -1i];

Open_Loop_Poles=cplxpair(pole(Canon_Godys));

% 14) Start optimization time count for partial allocation
Kpf_theorem_7=zeros(n_u,n_x);

% 15) Counting structure for multi-step partial allocation
tic
for j = 1:5
    % 15.1) Obtain the Eigenvalue and Eigenvector Matrices
    [Right_Eigenvectors, Eigenvalues, Left_Eigenvectors] = ...
        eig(A_x + B_u * Kpf_theorem_7);

    % 15.2) Define the poles x
    pole1 = Open_Loop_Poles(2 * j - 1, 1);
    pole2 = Open_Loop_Poles(2 * j, 1);

    % 15.3) Calculate the absolute difference between each element of
    % Lambda and poles
    dif_pole1 = abs(Eigenvalues - pole1);
    dif_pole2 = abs(Eigenvalues - pole2);

    % 15.4) Find the index of the element with the smallest difference
    [~, minIndex_1] = min(dif_pole1(:));
    [~, minIndex_2] = min(dif_pole2(:));

    % 15.5) Convert the linear index to matrix indices
    [row_1, col_1] = ind2sub(size(Eigenvalues), minIndex_1);
    [row_2, col_2] = ind2sub(size(Eigenvalues), minIndex_2);

    % 15.6) Group the indices
    row_poles = min(row_1, row_2):max(row_1, row_2);
    col_poles = min(col_1, col_2):max(col_1, col_2);

    % 15.7) Find the matrices for partial allocation
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

        if ~isempty(C_y)
            Tilde_Cy = C_y*L_p*inv(L_p'*L_p)*Q'/2;
        else
            Tilde_Cy=[];
        end
        if ~isempty(C_z)
            Tilde_Cz = C_z*L_p*inv(L_p'*L_p)*Q'/2;
        else
            Tilde_Cz=[];
        end

        % 15.10.2) Make Decision Variables 
        
        Tilde_X  = sdpvar(n_ola,n_ola,'symmetric');
        Tilde_W  = sdpvar(n_u,n_ola,'full');
        
        % 15.10.3) Make LMIs for theorem 6
        set_LMIs_theorem_7 = Tilde_X>=eps*eye(n_ola);

        % 15.10.4) Make LMI for Re(s)<=-alpha_v
        if ~isempty(alpha_v)
            Tilde_Ax_Valpha = Tilde_Lambda_p+alpha_v*eye(n_ola);
            set_LMIs_theorem_7 = [ set_LMIs_theorem_7,...
                Tilde_Ax_Valpha*Tilde_X+Tilde_X*Tilde_Ax_Valpha'+...
                Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu'<=-eps*eye(n_ola)  ];
        end

        % 15.10.5) Make LMI for Re(s)>=-beta_v
        if ~isempty(beta_v)
            Tilde_Ax_Vbeta = -Tilde_Lambda_p-beta_v*eye(n_ola);
            set_LMIs_theorem_7 = [ set_LMIs_theorem_7,...
                Tilde_Ax_Vbeta*Tilde_X+Tilde_X*Tilde_Ax_Vbeta'-...
                Tilde_Bu*Tilde_W-Tilde_W'*Tilde_Bu'<=-eps*eye(n_ola) ];
        end

        % 15.10.6) Make LMI for abs(s)<=-r_d
        if ~isempty(r_d)
            Tilde_Ax_Dqr = Tilde_Lambda_p+q_d*eye(n_ola); 
            set_LMIs_theorem_7 = [ set_LMIs_theorem_7,...
                [ -r_d*Tilde_X,...
                  Tilde_Ax_Dqr*Tilde_X+Tilde_Bu*Tilde_W;...
                  Tilde_X*Tilde_Ax_Dqr'+Tilde_W'*Tilde_Bu',...
                  -r_d*Tilde_X ]...
                  <=-eps*eye(2*n_ola)    ];
        end

        % 15.10.7) Make LMI for -theta_s<=angle(s)<=theta_s
        if ~isempty(theta_s)
            Tilde_Axsin = Tilde_Lambda_p*sin(theta_s);
            Tilde_Axcos = Tilde_Lambda_p*cos(theta_s);
            Tilde_Busin = Tilde_Bu*sin(theta_s);
            Tilde_Bucos = Tilde_Bu*cos(theta_s);
            set_LMIs_theorem_7 =[ set_LMIs_theorem_7,...
                [ Tilde_Axsin*Tilde_X+Tilde_X*Tilde_Axsin'+...
                  Tilde_Busin*Tilde_W+Tilde_W'*Tilde_Busin',...
                  Tilde_Axcos*Tilde_X-Tilde_X*Tilde_Axcos'+...
                  Tilde_Bucos*Tilde_W-Tilde_W'*Tilde_Bucos';...
                  Tilde_X*Tilde_Axcos'-Tilde_Axcos*Tilde_X+...
                  Tilde_W'*Tilde_Bucos'-Tilde_Bucos*Tilde_W,...
                  Tilde_Axsin*Tilde_X+Tilde_X*Tilde_Axsin'+...
                  Tilde_Busin*Tilde_W+Tilde_W'*Tilde_Busin'...
                ]<=-eps*eye(2*n_ola) ];
        end

        % 15.10.9) Make LMI for -w_H<=imag(s)<=w_H
        if ~isempty(w_H)
            LMIs_CMHD =[ set_LMIs_theorem_7,...
                [ -w_H*Tilde_X,...
                  Tilde_X*Tilde_Lambda_p'-Tilde_Lambda_p*Tilde_X+...
                  Tilde_W'*Tilde_Bu'-Tilde_Bu*Tilde_W;...
                  Tilde_Lambda_p*Tilde_X-Tilde_X*Tilde_Lambda_p'+...
                  Tilde_Bu*Tilde_W-Tilde_W'*Tilde_Bu',...
                  -w_H*Tilde_X...
                ]<=-eps*eye(2*n_ola) ];
        end 

        % 15.10.10) Make LMI for H2-Control
        if ~isempty(c_H2)

            Z = sdpvar(n_y,n_y,'symmetric');
            rho = sdpvar(1,1,'symmetric');   
            
            set_LMIs_theorem_7 = [ set_LMIs_theorem_7, ...
                Tilde_Lambda_p*Tilde_X+Tilde_X*Tilde_Lambda_p'+...
                Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu'+...
                Tilde_Bd*Tilde_Bd'<=-eps*eye(n_ola) ];
            
            set_LMIs_theorem_7 = [ set_LMIs_theorem_7,...
                        Z>=eps*eye(n_y),...
                        trace(Z)<=rho,...
                        rho>=eps ];
                    
            set_LMIs_theorem_7 = [set_LMIs_theorem_7,...
                [ -Z, Tilde_Cy*Tilde_X+D_y*Tilde_W;...
                    Tilde_X'*Tilde_Cy'+Tilde_W'*D_y', -Tilde_X ]<=...
                    -eps*eye(n_y+n_ola,n_y+n_ola) ]        
        end

        % 15.10.11) Make LMI for Hinf-Control
        if ~isempty(c_Hinf)
            gamma = sdpvar(1,1,'symmetric');   
            set_LMIs_theorem_7 = [ set_LMIs_theorem_7,...
                gamma>=eps,...
                [   Tilde_Lambda_p*Tilde_X+Tilde_X*Tilde_Lambda_p'+...
                    Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu',...
                    Tilde_Bd, Tilde_X*Tilde_Cz'+Tilde_W'*D_z';...
                    Tilde_Bd', -gamma*eye(n_d), E_z';...
                    Tilde_Cz*Tilde_X+D_z*Tilde_W,...
                    E_z, -gamma*eye(n_z)...
                ]<=-eps*eye(n_ola+n_d+n_z)...
            ];
        end

        % 15.10.12) Execute Optimization   
        if ~isempty(c_H2)
             if ~isempty(c_Hinf)
                optimize(set_LMIs_theorem_7,...
                    c_H2*rho+c_Hinf*gamma,Yalmip_sdpsettings); 
             else
                 optimize(set_LMIs_theorem_7,...
                     c_H2*rho,Yalmip_sdpsettings);
             end
        else
            if ~isempty(c_Hinf)
                optimize(set_LMIs_theorem_7,...
                    c_Hinf*gamma,Yalmip_sdpsettings); 
             else
                 optimize(set_LMIs_theorem_7,[],Yalmip_sdpsettings);
             end
        end

        % 15.10.13) Extract optimal variables 
        optimal_Tilde_W = value(Tilde_W);
        optimal_Tilde_X = value(Tilde_X);
        
        % 15.10.14) Compute and update feedback matrix
        Kpf_parcial=optimal_Tilde_W*inv(optimal_Tilde_X);
        Kpf_theorem_7 = Kpf_theorem_7 + Kpf_parcial * (Q * L_p') / 2;

    else
        disp('uncontrollable pole detected')
    end

end

% 15.20) End multi-step optimization time
time_multistep_theorem_7 = toc;

% 16) Make closed loop system d(t) to y(t) via theorem 2 
system_y_theorem_2=...
ss(A_x+B_u*Kpf_theorem_2,B_d,C_y+D_y*Kpf_theorem_2,E_y);

% 17) Make closed loop system d(t) to z(t) via theorem 2 
system_z_theorem_2=...
ss(A_x+B_u*Kpf_theorem_2,B_d,C_z+D_z*Kpf_theorem_2,E_z); 

% 18) Make closed loop system d(t) to y(t) via multistep theorem 6
system_y_theorem_7=...
ss(A_x+B_u*Kpf_theorem_7,B_d,C_y+D_y*Kpf_theorem_7,E_y);

% 19) Make closed loop system d(t) to z(t) via multistep theorem 6
system_z_theorem_7=...
ss(A_x+B_u*Kpf_theorem_7,B_d,C_z+D_z*Kpf_theorem_7,E_z); 

% 20) Closed-Loop system poles via theorem 2
Poles_system_theorem_2=cplxpair(pole(system_y_theorem_2));
    
% 21) Closed-Loop system poles via theorem 6   
Poles_system_theorem_7=cplxpair(pole(system_y_theorem_7));
       
    
    Table_Poles=table(Open_Loop_Poles,...
               Poles_system_theorem_2,...
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
    norm_H2(2,1)=norm(system_y_theorem_2,2);

%       24.2) variation of H2 norm for Gdy in closed loop            
discrepance_norm_H2(2,1)=...
    100*(norm_H2(2,1)-norm_H2(1,1))/norm_H2(1,1);

%       24.3) H-infinity norm for Gdy in closed loop            
norm_Hinf(2,1)=norm(system_z_theorem_2,'inf');

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
    time_theorem_2; 
    
%       26.3) Optimization time for multistep theorem 6
Optimization_Time(3,1)=...
    time_multistep_theorem_7; 

%   27) Table to compare norms
Transfer_Functions=["Open Loop";"Theorem 2";...
                    "Multistep Theorem 6"];
Table_Norms = table(Transfer_Functions,...
norm_H2,discrepance_norm_H2,...
norm_Hinf,discrepance_norm_Hinf,...
Optimization_Time),

%   28) Table - Optimization Time
percentual_discrepance_optimization_time=...
    100*(time_multistep_theorem_7-time_theorem_2)/time_theorem_2;

Table_Time_Optimization=table(time_theorem_2,...
    time_multistep_theorem_7,...
    percentual_discrepance_optimization_time),


% 29) Plot for Poles and Position

%   29.1) Eigenvalues (Poles) to be inserted in the plot
eigenvalues = [ Open_Loop_Poles',...
               Poles_system_theorem_2',...
               Poles_system_theorem_7'];

%   29.2) Computing boundary curves of the regions
% Inicializando arrays        
x_horizontal_strip=[];
y_inferior_horizontal_strip=[];
y_superior_horizontal_strip=[];
x_sector=[] ;
y_inferior_sector=[] ;
y_superior_sector=[] ;
x_parabola=[] ;
y_inferior_parabola=[] ;
y_superior_parabola=[] ;
x_disk=[] ;
y_disk=[] ;
x_alphav=[] ;
y_alphav=[] ;
x_betav=[] ;
y_betav=[] ;            
            
% 29.3) Defining the largest and smallest real part of the poles of the system
min_real_aut=min(real(eigenvalues)); 
max_real_aut=max(real(eigenvalues));                

% 29.4) Compute the largest and smallest real part of the disk region
min_real_disk=-q_d-r_d;
max_real_disk=-q_d+r_d;

29.5) Compute the largest and smallest real parts of the region of the coordinate points
min_real=min([min_real_aut,min_real_disk,-beta_v]);

max_real=max([max_real_aut,max_real_disk,-alpha_v]);
max_real=min([max_real,0]);

% 29.6) Compute the largest imaginary part of the poles of the system
max_imag_aut=max(imag(eigenvalues));

% 29.7) Compute the largest imaginary part of the sector region
max_imag_sector=-min_real*tan(theta_s);

% 29.8) Compute the largest imaginary part of the coordinate points
max_imag=max([max_imag_aut, max_imag_sector, w_H, r_d]);

% 29.9) Determine the coordinate points corresponding to the boundary curves
for i=1:100
    xt=min_real+(i-1)*(max_real-min_real)/99;
    if ~isempty(w_H)
        x_horizontal_strip(i,1)=xt;
        y_inferior_horizontal_strip(i,1)=-w_H;
        y_superior_horizontal_strip(i,1)=w_H;
    end
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
         real(Poles_system_theorem_2),imag(Poles_system_theorem_2),'xr',...
         real(Poles_system_theorem_7),imag(Poles_system_theorem_7),'xb',...
         x_alphav,y_alphav,'-c',...
         x_betav,y_betav,'-c',...
         x_sector,y_inferior_sector,'-c',...
         x_sector,y_superior_sector,'-c',...
         x_disk,y_disk,'-c',...
         x_horizontal_strip,y_inferior_horizontal_strip,'-c',...
         x_horizontal_strip,y_superior_horizontal_strip,'-c',...
         'LineWidth',2,'MarkerSize',8);
xlim([-0.81 0]);
ylim([-1.3 1.3]);
title('Complex Plane','FontSize',12)
legend('Open','Theorem 2','Theorem 7','Frontier','FontSize',10)
    subplot(212)
    plot(real(Open_Loop_Poles),imag(Open_Loop_Poles),'ok',...
         real(Poles_system_theorem_2),imag(Poles_system_theorem_2),'xr',...
         real(Poles_system_theorem_7),imag(Poles_system_theorem_7),'xb',...
         x_alphav,y_alphav,'-c',...
         x_betav,y_betav,'-c',...
         x_sector,y_inferior_sector,'-c',...
         x_sector,y_superior_sector,'-c',...
         x_disk,y_disk,'-c',...
         x_horizontal_strip,y_inferior_horizontal_strip,'-c',...
         x_horizontal_strip,y_superior_horizontal_strip,'-c',...
         'LineWidth',2,'MarkerSize',8);
xlim([-0.46 -0.455]);
ylim([-0.5295 -0.527]);
title('Complex Plane','FontSize',12)
legend('Open','Theorem 2','Theorem 7','Frontier','FontSize',10)

%   30) Numerical Integral Table of the Absolute Value of the Control Action
%      30.1) Transient Response - free allocation
[y_theorem_2,t_theorem_2,x_theorem_2] = ...
    step(system_y_theorem_2,0:0.01:14);

%      30.2) Control action - free allocation         
u_theorem_2(1:size(x_theorem_2,1),1:n_u,1)=...
    x_theorem_2(:,:,1)*Kpf_theorem_2';

%      30.3) Transient Response - partial multistep allocation                                 
[y_theorem_7_multistep,t_theorem_7_multistep,x_theorem_7_multistep] = ...
    step(system_y_theorem_7,0:0.01:14);

%      30.4) Control action - partial multistep allocation 
u_theorem_7_multistep(1:size(x_theorem_2,1),1:n_u,1)=...
    x_theorem_7_multistep(:,:,1)*Kpf_theorem_7';

%   31) Figure - control actions
figure
%       31.1) Response u1(t) to step input d1(t)    
subplot(211)
plot(t_theorem_2,u_theorem_2(1:end,1),'-r',...
     t_theorem_7_multistep,u_theorem_7_multistep(1:end,1),'-b')
title('u_1(t) before step d(t)','FontSize',12)
legend( 'theorem 2','theorem 7','FontSize',10)

%       31.2) Response u2(t) to step input d1(t)            
subplot(212)
plot(t_theorem_2,u_theorem_2(1:end,2),'-r',...
     t_theorem_7_multistep,u_theorem_7_multistep(1:end,2),'-b')
title('u_2(t) before step d(t)','FontSize',12)
legend('theorem 2','theorem 7','FontSize',10)

