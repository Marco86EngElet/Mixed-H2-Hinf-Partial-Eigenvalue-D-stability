%A) Obtain the Partial Controllability Matrix
 

    controlability_matrix = ctrb(Lambda_p, (Q_p * L_p') * B_u / 2);

%B) Calculate the number of uncontrollable states

    number_of_non_controlable_states = 2-rank(controlability_matrix);

    %C) Check if it is controllable
    if number_of_non_controlable_states == 0 

        %D) Compute constant matrices resulting from similarity transformations

        Tilde_Lambda_p = Q_p*Lambda_p*Q_p'/2;
        Tilde_Bu = Q_p*L_p'*B_u/2;
        Tilde_Bd = Q_p*L_p'*B_d/2;

        if ~isempty(C_y)

            Tilde_Cy = C_y*L_p*inv(L_p'*L_p)*Q_p'/2;

        else
            Tilde_Cy=[];
        end
        if ~isempty(C_z)

            Tilde_Cz = C_z*L_p*inv(L_p'*L_p)*Q_p'/2;

        else

            Tilde_Cz=[];

        end

        %E) Make decision variables

        Tilde_X  = sdpvar(n_ola,n_ola,'symmetric');
        Tilde_W  = sdpvar(n_u,n_ola,'full');
        rho= sdpvar(1,1,'full');
        gamma=sdpvar(1,1,'full');
        
        %F) Initiate LMIs for Algorithm 7 to allocate the j-th pair of complex poles
        
        set_LMIs_partial = Tilde_X>=eps*eye(n_ola);
        set_LMIs_partial = [set_LMIs_partial,rho>=eps,gamma>=eps];
            
        %G) LMI for "Re(s) >= -alpha_v" (Vertical Strip Region)
        if ~isempty(alpha_v)

            Tilde_Ax_Valpha = Tilde_Lambda_p+alpha_v*eye(n_ola);

            set_LMIs_partial = [ set_LMIs_partial,...
                Tilde_Ax_Valpha*Tilde_X+Tilde_X*Tilde_Ax_Valpha'+...
                Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu'<=-eps*eye(n_ola)  ];
        end

        %H) LMI for "Re(s) >= -beta_v" (Vertical Strip Region)

        if ~isempty(beta_v)

            Tilde_Ax_Vbeta = -Tilde_Lambda_p-beta_v*eye(n_ola);

            set_LMIs_partial = [ set_LMIs_partial,...
                Tilde_Ax_Vbeta*Tilde_X+Tilde_X*Tilde_Ax_Vbeta'-...
                Tilde_Bu*Tilde_W-Tilde_W'*Tilde_Bu'<=-eps*eye(n_ola) ];
        end

        %I) LMI for "abs(s) <= r_d" (Disk Region)
        if ~isempty(r_d)

            Tilde_Ax_Dqr = Tilde_Lambda_p+q_d*eye(n_ola); 

            set_LMIs_partial = [ set_LMIs_partial,...
                [ -r_d*Tilde_X,...
                  Tilde_Ax_Dqr*Tilde_X+Tilde_Bu*Tilde_W;...
                  Tilde_X*Tilde_Ax_Dqr'+Tilde_W'*Tilde_Bu',...
                  -r_d*Tilde_X ]...
                  <=-eps*eye(2*n_ola)    ];
        end

        %J) LMI for "-theta_s ≤ angle(s) ≤ theta_s" (Sector Region)

        if ~isempty(theta_s)

            Tilde_Axsin = Tilde_Lambda_p*sin(theta_s);

            Tilde_Axcos = Tilde_Lambda_p*cos(theta_s);

            Tilde_Busin = Tilde_Bu*sin(theta_s);

            Tilde_Bucos = Tilde_Bu*cos(theta_s);

            set_LMIs_partial =[ set_LMIs_partial,...
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

        %K) LMI for "-w_H ≤ imag(s) ≤ w_H" (Horizontal Strip Region)
        if ~isempty(w_H)
            set_LMIs_partial =[ set_LMIs_partial,...
                [ -w_H*Tilde_X,...
                  Tilde_X*Tilde_Lambda_p'-Tilde_Lambda_p*Tilde_X+...
                  Tilde_W'*Tilde_Bu'-Tilde_Bu*Tilde_W;...
                  Tilde_Lambda_p*Tilde_X-Tilde_X*Tilde_Lambda_p'+...
                  Tilde_Bu*Tilde_W-Tilde_W'*Tilde_Bu',...
                  -w_H*Tilde_X...
                ]<=-eps*eye(2*n_ola) ];
        end 

        %M) LMI for H_2 Optimal Control

        if ~isempty(c_H2)

            Z = sdpvar(n_y,n_y,'symmetric');

            rho = sdpvar(1,1,'symmetric');   

            set_LMIs_partial = [ set_LMIs_partial, ...
                Tilde_Lambda_p*Tilde_X+Tilde_X*Tilde_Lambda_p'+...
                Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu'+...
                Tilde_Bd*Tilde_Bd'<=-eps*eye(n_ola) ];

            set_LMIs_partial = [ set_LMIs_partial,...
                        Z>=eps*eye(n_y),...
                        trace(Z)<=rho,...
                        rho>=eps ];

            set_LMIs_partial = [set_LMIs_partial,...
                [ -Z, Tilde_Cy*Tilde_X+D_y*Tilde_W;...
                    Tilde_X'*Tilde_Cy'+Tilde_W'*D_y', -Tilde_X ]<=...
                    -eps*eye(n_y+n_ola,n_y+n_ola) ]        
        end

        %N) LMI for H-infinity Optimal Control

        if ~isempty(c_Hinf)

            gamma = sdpvar(1,1,'symmetric');   

            set_LMIs_partial = [ set_LMIs_partial,...
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

        else 
            disp('SREAPE incontrolavel')
        end        
