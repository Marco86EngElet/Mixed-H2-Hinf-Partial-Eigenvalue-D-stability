clc, close, clear

format SHORT

%% 1) Inputs and Outputs

    % 1.1) Dinamics states
    A_x=toeplitz([0.1:0.01:0.5,0.01]);
    n_x=size(A_x,1);
    
    % 1.2) Number of outputs not directly influenced by disturbances.
    n_y=1;
    
    % 1.3) Number of outputs influenced by disturbances.
    n_z=1;

    % 1.) Number of disturbs
    n_d=5;    

    % 1.) Number of inputs control 
    n_u=2;

    % ) Number of allocable poles by feedback
    n_ola=n_u;
 
%% 2) Enter the number of real poles and complex conjugate pole pairs 

    %2.1)  Number of real poles
    n_r=n_x;
    
    %2.2) Number of complex conjugate pole pairs 
    n_pcc=(n_x-2)/2;

%% 3) Parameters for D-stability LMIs
    alpha_v=1.2321;
    beta_v=1.6667;
    e_P=[];
    q_d=0;
    r_d=1.9216;
    w_H=2.4993;
    theta_s=0.8092;
    
%% 4) Canonical Representation of Open-Loop System Transfer Functions
    
    %4.1) B_d
    B_d=0.1*eye(n_x,n_d);
    
    %4.2) C_y
    C_y=zeros(n_y,n_x);
    C_y(1,end)=1;
    C_y(1,end-1)=-0.5;
    C_y(1,end-2)=0.25;
    
    %4.3) C_z
    n_z=1;
    C_z=zeros(n_z,n_x);
    C_z(1,1)=1;
    C_z(1,end)=-1;
    
    %4.4) E_y
    E_y=zeros(n_y,n_d);
        
    %5.4) E_z
    E_z=zeros(n_z,n_d);
    E_z(1,1)=0.01;
    E_z(1,2)=-0.01;
    E_z(1,3)=0.01;
    E_z(1,4)=-0.01;
    E_z(1,5)=0.01;
    
    %4.5) Canonical Representation of Gody(s)
    
    Canon_Godys = ss(A_x,B_d,C_y,E_y);
    
    %4.6) Canonical Representation of Godz(s)
    
    Canon_Godzs = ss(A_x,B_d,C_z,E_z);

%% 5) Compute Norms and Poles of Open-Loop System
    
    Compute_Norms_and_Poles_Open_Loop    
    
%% 6) Compute Natural Frequency and Damp Ratio

    [wn_God,zeta_God,p_God]=damp(Canon_Godys);

%% 7) Gody(s) Response of unit impulse disturb

    impulse_simulation_time=0:0.03:30;
    
    Compute_Impulse_Simulation_Data
    
%% 8) Input Control Matrices: 'B_u', 'D_y', 'D_z'.
    
    %8.1) B_u
    
    B_u=zeros(n_x,n_u);
    B_u(2,1)=1;
    B_u(3,1)=1;
    B_u(2,2)=-1;
    B_u(3,2)=-1;
    
    %8.2) D_y e D_z
    
    D_y=0.001*ones(n_y,n_u);
    D_z=D_y;
    
%% 9) Configure SDP optimization

    %9.1) Base yalmip settings for SDP
    
    SDP_settings =...
         sdpsettings('verbose',1,'solver','lmilab','debug',1);

    %9.2) Cost Function Weights 
    
    c_H2=[];
    c_Hinf=[];

%% 10) Execute Theorem 1

    %11.1) Instantiate feedback matrix
    
    Kpf_partial=zeros(n_u,n_x);

    %11.2) Input Q_p
    
    Q_p = sqrt(2)*eye(2);

    %11.3) Compute Eigenvalues and Eigenvectors
        
    [Right_Eigenvectors, Eigenvalues, Left_Eigenvectors] = ...
        eig(A_x + B_u * Kpf_partial);
        
    %11.4) Find poles to make matrices for allocation
    [row_poles,col_poles]=find(real(Eigenvalues)>0);
    
    %11.5) Compute matrices for allocation
    Lambda_p = Eigenvalues(row_poles, col_poles);
    L_p = Left_Eigenvectors(1:end, col_poles);    
    
    %11.6) Compute D-stable control via partial pole placement
    LMIs_mixed_control_partial_D_stable

        %11.7) Start time count for optimization via Theorem 1
        tic
    
        %11.8) Compute theorem 1   
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
    
		%11.9) End time for Theorem 1
		classic_optimization_time=toc;
    
        %11.10) Extrair variaveis de decisao para computar K_pf 

        optimal_Tilde_W = value(Tilde_W);
        optimal_Tilde_X = value(Tilde_X);    

%% 12) Compute Feedback Matrix via Theorem 1

Kpf_classic=optimal_Tilde_W/optimal_Tilde_X;
Kpf_classic=Kpf_classic* (Q_p * L_p')/2;

%% 13) State-space representation of the closed-loop system designed via Theorem 1

    % 13.1) Gcdy(s) 
    
    Gcdys_classic=...
        ss(A_x+B_u*Kpf_classic,B_d,C_y+D_y*Kpf_classic,E_y);

    % 13.2) Gcdz(s)  
    
    Gcdzs_classic=...
        ss(A_x+B_u*Kpf_classic,B_d,C_z+D_z*Kpf_classic,E_z);
    
%% 14) Compute H_2 and H_inf Norms for the system designed via Theorem 1

    %14.1) H_2 norm of Gcdy(s)
    
    H2_norm_Gcdys_classic=norm(Gcdys_classic,2);
    
    %14.2) H_inf norm of Gcdz(s)
    
    Hinf_norm_Gcdzs_classic=norm(Gcdzs_classic,'inf');
    
%% 15) Compute the poles of the system designed via Theorem 1

    Poles_Gcd_classic=cplxpair(pole(Gcdys_classic));
    
%% 16) Compute Gcdy(s) responses to a unit impulse disturbance.

    Compute_yt_ut_impulse_response_classic_method

%% 17) Execute theorem 7 
        
        %17.1) Cost Function Wheights
        c_H2=1;
        c_Hinf=2;
        
        %17.2) Start optimization time
        tic
    
        %17.3) Compute theorem 7   
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

    % 17.4) End optimization time
    partial_optimization_time = toc;
    
    %17.5) Extract optimal solution        
        optimal_Tilde_W = value(Tilde_W);
        optimal_Tilde_X = value(Tilde_X); 
        
    %17.6) Compute feedback matrix    
        KD_partial=optimal_Tilde_W*inv(optimal_Tilde_X);
        Kpf_partial = Kpf_partial + KD_partial * (Q_p * L_p') / 2;

        
%% 18) State-space representation of the closed-loop system designed via Theorem 7

    %18.1) Gcdy(s) 
    
    Gcdys_partial=...
        ss(A_x+B_u*Kpf_partial,B_d,C_y+D_y*Kpf_partial,E_y);

    %18.2) Gcdz(s) 
    
    Gcdzs_partial=...
        ss(A_x+B_u*Kpf_partial,B_d,C_z+D_z*Kpf_partial,E_z); 

%% 19) Compute H_2 and H-inf norms for the system designed via Theorem 7

    %19.1) H_2 Norm of Gcdy(s)
    
    H2_norm_Gcdys_partial=norm(Gcdys_partial,2);
    
    %19.2) H_inf Norm of Gcdz(s)
    
    Hinf_norm_Gcdzs_partial=norm(Gcdzs_partial,'inf');
    
%% 20) Compute the poles of the system designed via Theorem 7

    Poles_Gcd_Partial=cplxpair(pole(Gcdys_partial));
    
%% 21) Obtain Damping Ratios and Natural Frequencies

[wn_theo7,zeta_theo7,p_theo7]=damp(Gcdys_partial);
cont=1;
    
%% 22) Compute responses of Gcdy(s) to unit impulse disturbance

    Compute_yt_ut_impulse_response_partial_method

%% 23) Make Tables

    Compute_Tables
          
%% 24) Bode plots for open-loop and closed-loop LTI-MIMO system

    figure 
    subplot(121)
    sigmaplot(Gcdzs_classic,Gcdzs_partial,Canon_Godzs)
    title('')
    legend('Gcdz(s) Classic','Gcdz(s) Partial','Godz(s)','FontSize',10)
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); % Set thickness for all lines

    subplot(122)
    sigmaplot(Gcdzs_classic,Gcdzs_partial)
    title('')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); % Set thickness for all lines

    sgtitle('Singular Values G_d_z(s)')
    
%% 25) Transient response plots of g(t)^2 to unit impulse disturbance of the open-loop and closed-loop LTI-MIMO system

    figure
    
    subplot(151)
    plot(t_Gcdys_impulse_classic,y2_Gcdys_impulse_classic(1:end,1,1),...
         t_Gcdys_impulse_partial,y2_Gcdys_impulse_partial(1:end,1,1));
    title('y^2(t) X d_1(t)')
    xlabel('t')
    ylabel('y^2(t)')
    xlim([0 7])
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(152)
    plot(t_Gcdys_impulse_classic,y2_Gcdys_impulse_classic(1:end,1,2),...
         t_Gcdys_impulse_partial,y2_Gcdys_impulse_partial(1:end,1,2));
    title('y^2(t) X d_2(t)')
    xlabel('t')
    ylabel('y^2(t)')
    xlim([0 7])
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(153)
    plot(t_Gcdys_impulse_classic,y2_Gcdys_impulse_classic(1:end,1,3),...
         t_Gcdys_impulse_partial,y2_Gcdys_impulse_partial(1:end,1,3));
    title('y^2(t) X d_3(t)')
    xlim([0 7])
    xlabel('t')
    ylabel('y^2(t)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(154)
    plot(t_Gcdys_impulse_classic,y2_Gcdys_impulse_classic(1:end,1,4),...
         t_Gcdys_impulse_partial,y2_Gcdys_impulse_partial(1:end,1,4));
    title('y^2(t) X d_4(t)')
    xlabel('t')
    xlim([0 7])
    ylabel('y^2(t)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(155)
    plot(t_Gcdys_impulse_classic,y2_Gcdys_impulse_classic(1:end,1,5),...
         t_Gcdys_impulse_partial,y2_Gcdys_impulse_partial(1:end,1,5));
    title('y^2(t) X d_5(t)')
    xlabel('t')
    ylabel('y^2(t)')
    xlim([0 7])
    legend('Gcdy(s) classic','Gcdy(s) partial')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('Response y^2(t) to the impulse')
    
%% 26) Plot of u(t) for the open-loop and closed-loop LTI-MIMO system in response to a unit step disturbance

    
    figure
    subplot(121)
    semilogy(t_Gcdys_impulse_classic,u2_Gcdys_impulse_classic(:,1),...
         t_Gcdys_impulse_partial,u2_Gcdys_impulse_partial(:,1));
    title('Response u_1(t) to impulse ')
    xlabel('t')
    ylabel('u^2_1(t)')
    xlim([0,20]);
    legend('G_c_d_y(s) classic','G_c_d_y(s) partial')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(122)
    semilogy(t_Gcdys_impulse_classic,u2_Gcdys_impulse_classic(:,2),...
         t_Gcdys_impulse_partial,u2_Gcdys_impulse_partial(:,2));
    title('u^2_2(t) to impulse')
    xlabel('t')
    ylabel('u_2(t)')
    xlim([0,20]);
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('Response u(t) to impulse')

%% 27) Plot for Pole Location

    %27.1) Eigenvalues (Poles) to be inserted into the plot

    All_poles=[ Open_Loop_Poles',...
                  Poles_Gcd_classic',...  
                  Poles_Gcd_Partial'];

    %27.2) Computing boundary curves of the regions
   
    [  x_horizontal_strip,y_lower_horizontal_strip,...
        y_upper_horizontal_strip,...
        x_sector,y_lower_sector,y_upper_sector,...
        x_disk,y_disk,...
        x_alphav,y_alphav,...
        x_betav,y_betav ] = coordinates_D_regions(...
        All_poles,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P); 

    %27.3) Graph            
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
   ylim([-1.1 1.1])
   xlim([-2 0])
   grid on
   sgtitle('Poles and Location in the Complex Plane','FontSize',12)
   
