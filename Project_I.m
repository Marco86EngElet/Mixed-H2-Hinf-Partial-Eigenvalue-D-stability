clc, close, clear

format SHORT
%% 1) Inputs and Outputs

    % 1.1) Number of dynamics states
    n_x=10;
    
    % 1.2) Number of outputs without disturbs
    n_y=1;
    
    % 1.3) Number of outputs with disturs
    n_z=1;

    % 1.4) Number of disturbs
    n_d=1;    

    % 1.5) Number of input control 
    n_u=2;

    % 1.6) Number of allocable poles by feedback
    n_ola=n_u;
 
%% 2) Specify type of allocable poles

    %2.1) Number of real poles
    n_r=0;
    
    %2.2) Number of pair complex conjugate poles
    n_pcc=5;

%3) Parameters of D-stability    
    alpha_v=0.4578;
    beta_v=0.8;
    e_P=[];
    q_d=0;
    r_d=0.6999;
    w_H=1.2566;
    theta_s=1.0974;
    
    
%% 4) Compute transfer open-loop transfer functions

    %4.1) Inputs: damping ration
    
    zeta_c(1,1) = 0.1;  % Partial function 1
    zeta_c(2,1) = 0.12; % Partial function 2
    zeta_c(3,1) = 0.15; % Partial function 3
    zeta_c(4,1) = 0.2;  % Partial function 4
    zeta_c(5,1) = 0.3;  % Partial function 5

    %4.2) Inputs: natural frequency
    
    w_c(1,1) = 0.1;  % Partial function 1
    w_c(2,1) = 0.12; % Partial function 2
    w_c(3,1) = 0.15; % Partial function 3
    w_c(4,1) = 0.2;  % Partial function 4
    w_c(5,1) = 0.3;  % Partial function 5

    %4.3) Inputs: Gain 
    
    KYc=ones(5,1); %for "Gody(s)" 
    
    KZc=ones(5,1); %for "Godz(s)"
    KZc(3,1)=-1; 

    %4.4) Make "Gody(s)" e "Godz(s)" objects 

        %4.4.1) Initializing
    
        Godys=tf(0,1); % Gody(s)=0;
        Godzs=tf(0.01,1); % Godz(s)=0.01;
        
        %4.4.2) "for" loop to make sum of partial functions 
        
        for c=1:5 
            num = w_c(c,1)^2; % Numerator 
            den = [1, 2*zeta_c(c,1)*w_c(c,1), w_c(c,1)^2]; % Denominator  
            G_c = tf(num, den); % Partial fraction without gain
            Godys=Godys+KYc(c,1)*G_c; % Add Partial fraction in "Gody(s)"
            Godzs=Godzs+KZc(c,1)*G_c; % Add Partial fraction in "Godz(s)"
        end

%% 5) Canoncial representation of Open-Loop System 
    
    %5.1) Gody(s) 
    
    Canon_Godys = ss(Godys);
    
    %5.2) Godz(s)
    
    Canon_Godzs = ss(Godzs);
    
    %5.3) Extract matrices
    
    A_x=Canon_Godys.A;
    B_d=Canon_Godys.B;
    C_y=Canon_Godys.C;
    E_y=Canon_Godys.D;
    C_z=Canon_Godzs.C;
    E_z=Canon_Godzs.D; 

%% 6) Compute Norms H2 and Hinf of open-loop system
    
    Compute_Norms_and_Poles_Open_Loop
    
%% 7) Compute natural frequency and damp ration of open-loop system

    [wn_God,zeta_God,p_God]=damp(Canon_Godys);
        
%% 8) Compute unit impulse response of Gody(s) 

    impulse_simulation_time=0:0.01:14;
    
    Compute_Impulse_Simulation_Data
    
%% 9) Input control matrices 'B_u', 'D_y', 'D_z'.
    
    %9.1) B_u
    
    B_u=zeros(n_x,n_u);
    B_u(2,1)=1;
    B_u(3,1)=1;
    B_u(2,2)=-1;
    B_u(3,2)=-1;
    
    %9.2) D_y e D_z
    
    D_y=0.001*ones(n_y,n_u);
    D_z=D_y;
    
%% 10) Yalmip Settings for SDP optimization

    %10.1) Base configuration
    
    SDP_settings =...
         sdpsettings('verbose',1,'solver','lmilab','debug',1);

    %10.2) Cost Function Weights 
    
    c_H2=1;
    c_Hinf=2;

%% 11) Compute number of uncontrollabel systems states
    
    Compute_number_of_unallocatable_states

%% 12)  Make decision variables for SDP-LMI optimization

X  = sdpvar(n_x,n_x,'symmetric');
W  = sdpvar(n_u,n_x,'full');
gamma = sdpvar(1,1,'symmetric'); 
Z = sdpvar(n_y,n_y,'symmetric');
rho = sdpvar(1,1,'symmetric');   
        
%% 13) Make LMIs for optimal mixed H2-Hinf control with D-stability

LMIs_mixed_control_D_stable

%% 14) Compute Theorem 2

tic %Start time counter
optimize(set_LMIs_classic,...
            c_H2*rho+c_Hinf*gamma,...
            SDP_settings);
classic_optimization_time=toc; %End time counter

%% 15) Extract optimal solution by Theorem 2

optimal_W = value(W);
optimal_X = value(X);

%% 16) Compute Feedback matrix by Theorem 2

Kpf_classic=optimal_W/optimal_X;

%% 17) State Space for Closed Loop System by Theorem 2

    % 17.1) State Space Gcdy(s)
    
    Gcdys_classic=...
        ss(A_x+B_u*Kpf_classic,B_d,C_y+D_y*Kpf_classic,E_y);

    % 17.2) State Space Gcdz(s)  
    
    Gcdzs_classic=...
        ss(A_x+B_u*Kpf_classic,B_d,C_z+D_z*Kpf_classic,E_z);

    
%% 18) Compute norms H_2 e H_inf of Closed Loop System by Theorem 2 

    %18.1) H_2 norm of Gcdy(s)
    
    H2_norm_Gcdys_classic=norm(Gcdys_classic,2);
    
    %18.2) H_inf norm of Gcdz(s)
    
    Hinf_norm_Gcdzs_classic=norm(Gcdzs_classic,'inf');
    
%% 19) Poles of Closed Loop System by Theorem 2

    Poles_Gcd_classic=cplxpair(pole(Gcdys_classic));
    
%% 20) Response Gcdy(s) by unit impulse disturb

    Compute_yt_ut_impulse_response_classic_method
    
%% 21) Computer Theorem 7 multistep  

    %21.1) Initial Feedback matrix
    
    Kpf_partial=zeros(n_u,n_x);
    
    %21.2) Input Q_p  
    
    Q_p = [1 1;1i -1i];

    %21.3) Start optimization time counter 
    tic
    
    %21.4) Counter for multistep feedback
    
    for j = 1:n_pcc % counter depends of pair-conjugate complex poles
        
        %21.4.1) Computer Eigenvectors and Eigenvalues
        
        [Right_Eigenvectors, Eigenvalues, Left_Eigenvectors] = ...
            eig(A_x + B_u * Kpf_partial);

        %21.4.2) Chose j-th pair-conjugate complex poles    
        
        pole1 = Open_Loop_Poles(2 * j - 1, 1);
        pole2 = Open_Loop_Poles(2 * j, 1);

        %21.4.3) Compute the differences between eigenvalues and poles

        dif_pole1 = abs(Eigenvalues - pole1);
        dif_pole2 = abs(Eigenvalues - pole2);

        %21.4.4) Find the index of the elements with the smallest difference between the eigenvalues and the selected poles for allocation

        
        [~, minIndex_1] = min(dif_pole1(:));
        [~, minIndex_2] = min(dif_pole2(:));

        %21.4.5) Convert linear indices "minIndex\_1" and "minIndex\_2" 
        % of the smallest value to subscripts 
        % (row and column, stored in "row\_1", "row\_2", "col\_1", "col\_2") 
        % of the matrix "Eigenvalues".
        
        [row_1, col_1] = ind2sub(size(Eigenvalues), minIndex_1);
        [row_2, col_2] = ind2sub(size(Eigenvalues), minIndex_2);

        %21.4.6) Arrange the row and column intervals within "row_poles" and "col_poles."   
        
        row_poles = min(row_1, row_2):max(row_1, row_2);
        col_poles = min(col_1, col_2):max(col_1, col_2);

        %21.4.7) Find submatrices "Lambda_p" and "L_p"
        
        Lambda_p = Eigenvalues(row_poles, col_poles);
        L_p = Left_Eigenvectors(1:end, col_poles);

        %21.5) Compute mixed D-stable optimal control with partial allocation
        LMIs_mixed_control_partial_D_stable
    
            %21.6) Compute optimization    
            if ~isempty(c_H2) 
                 if ~isempty(c_Hinf) 
                     
                    optimize(set_LMIs_partial,...
                        c_H2*rho+c_Hinf*gamma,SDP_settings); 
                    
                 else 
                     
                     optimize(set_LMIs_partial,...
                         c_H2*rho,SDP_settings);
                     
                 end
            else 
                if ~isempty(c_Hinf) 
                    
                    optimize(set_LMIs_partial,...
                        c_Hinf*gamma,SDP_settings); 
                    
                else 
                    
                     optimize(set_LMIs_partial,[],SDP_settings);
                     
                end 
            end

            %21.7) Extract decision variable by Theorem 7
            
            optimal_Tilde_W = value(Tilde_W);
            optimal_Tilde_X = value(Tilde_X);

            %21.8) Update Feedback matrix
            
            KD_partial=optimal_Tilde_W*inv(optimal_Tilde_X);
            Kpf_partial = Kpf_partial + KD_partial * (Q_p * L_p') / 2;


    end %End Multistep Process
    
    % 21.9) End Optimization time counter
    
    partial_optimization_time = toc;
    
%% 22) State-space representation of the closed-loop system designed via multi-step Theorem 7.

    % 22.1) State-Space for Gcdy(s) 
    
    Gcdys_partial=...
        ss(A_x+B_u*Kpf_partial,B_d,C_y+D_y*Kpf_partial,E_y);

    % 22.2) State-Space for Gcdz(s) 
    
    Gcdzs_partial=...
        ss(A_x+B_u*Kpf_partial,B_d,C_z+D_z*Kpf_partial,E_z); 

%% 23) Compute H_2 and H_inf norms of closed-loop system designed via multi-step Theorem 7.

    %23.1) H_2 norm of Gcdy(s)
    
    H2_norm_Gcdys_partial=norm(Gcdys_partial,2);
    
    %23.2) H_inf norm of Gcdz(s)
    
    Hinf_norm_Gcdzs_partial=norm(Gcdzs_partial,'inf');
    
%% 24) Compute poles of closed-loop system designed via multi-step Theorem 7.

    Poles_Gcd_Partial=cplxpair(pole(Gcdys_partial));
    
%% 25) Obtain Damping Coefficients and Natural Frequency.

[wn_theo7,zeta_theo7,p_theo7]=damp(Gcdys_partial);
cont=1;

%% 26) Compute Gcdy(s) response of unit impulse response

    Compute_yt_ut_impulse_response_partial_method
    
%% 27) Make Tables

    Compute_Tables
          
%% 28) Construct Bode plots for open-loop and closed-loop LTI-MIMO systems.

    figure 
    sigmaplot(Gcdzs_classic,Gcdzs_partial,Canon_Godzs)
    title('Singular Value G_d_z(s)')
    legend('Gcdz(s) Classic','Gcdz(s) Partial',...
           'Godz(s)','FontSize',10)
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); % Set thickness for all lines
    
%% 29) Construct transient response plots of g(t)^2 to a unit impulse disturbance for open-loop and closed-loop LTI-MIMO systems.

    figure
    subplot(121)
    plot(t_Gcdys_impulse_classic,y2_Gcdys_impulse_classic,...
         t_Gcdys_impulse_partial,y2_Gcdys_impulse_partial,...
         t_Canon_Godys_impulse,y2_Canon_Godys_impulse);
    xlim([0 14])
    title('Impulse Response')
    xlabel('t')
    ylabel('y^2(t)')
    legend('Gcdy(s) classic','Gcdy(s) partial','Gody(s)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(122)
    plot(t_Gcdys_impulse_classic,y2_Gcdys_impulse_classic,...
         t_Gcdys_impulse_partial,y2_Gcdys_impulse_partial,...
         t_Canon_Godys_impulse,y2_Canon_Godys_impulse);
    xlim([0 5])
    ylim([0 0.008])
    title('Classic Vs Partial')
    xlabel('t')
    ylabel('y^2(t)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('Response of y(t) to the impulse')
    
%% 30) Construct a plot of u(t) for open-loop and closed-loop LTI-MIMO systems in response to a unit impulse disturbance.
    
    figure
    subplot(121)
    semilogy(t_Gcdys_impulse_classic,u2_Gcdys_impulse_classic(:,1),...
         t_Gcdys_impulse_partial,u2_Gcdys_impulse_partial(:,1));
    xlim([0,11]) 
    title('u_1(t)')
    xlabel('t')
    ylabel('u^2_1(t)')
    legend('G_c_d_y(s) classic','G_c_d_y(s) partial')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(122)
    semilogy(t_Gcdys_impulse_classic,u2_Gcdys_impulse_classic(:,2),...
         t_Gcdys_impulse_partial,u2_Gcdys_impulse_partial(:,2));
    xlim([0,11]) 
    title('u_2(t)')
    xlabel('t')
    ylabel('u^2_2(t)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('Response of u(t) to impulse')
    
%% 31) Pole Location Plot

    %31.1) Eigenvalues (Poles) to be plotted.
    All_poles=[ Open_Loop_Poles',...
                Poles_Gcd_classic',...  
                Poles_Gcd_Partial'];

    %31.2) Computing region boundary curves.    
    [  x_horizontal_strip,y_lower_horizontal_strip,...
        y_upper_horizontal_strip,...
        x_sector,y_lower_sector,y_upper_sector,...
        x_disk,y_disk,...
        x_alphav,y_alphav,...
        x_betav,y_betav ] = coordinates_D_regions(...
        All_poles,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P); 

    %31.3) Graph           
    figure
    subplot(121)
    plot(real(Open_Loop_Poles),imag(Open_Loop_Poles),'ok',...
         real(Poles_Gcd_classic),...
         imag(Poles_Gcd_classic),'+r',...
         real(Poles_Gcd_Partial),...
         imag(Poles_Gcd_Partial),'xb',...
         x_alphav,y_alphav,'-k',...
         x_betav,y_betav,'-k',...
         x_sector,y_lower_sector,'-k',...
         x_sector,y_upper_sector,'-k',...
         x_disk,y_disk,'-k',...
         x_horizontal_strip,y_lower_horizontal_strip,'-k',...
         x_horizontal_strip,y_upper_horizontal_strip,'-k',...
         'LineWidth',2,'MarkerSize',8);
   legend('Open Loop','Classic','Partial','Boundary','FontSize',10) 
   grid on
   
   subplot(122)
   plot(real(Open_Loop_Poles),imag(Open_Loop_Poles),'ok',...
         real(Poles_Gcd_classic),...
         imag(Poles_Gcd_classic),'+r',...
         real(Poles_Gcd_Partial),...
         imag(Poles_Gcd_Partial),'xb',...
         x_alphav,y_alphav,'-k',...
         x_betav,y_betav,'-k',...
         x_sector,y_lower_sector,'-k',...
         x_sector,y_upper_sector,'-k',...
         x_disk,y_disk,'-k',...
         x_horizontal_strip,y_lower_horizontal_strip,'-k',...
         x_horizontal_strip,y_upper_horizontal_strip,'-k',...
         'LineWidth',2,'MarkerSize',8);
   ylim([-0.6 0.6])
   xlim([-0.71 0])
   grid on
   sgtitle('Poles and Location in Complex Plane','FontSize',12)
