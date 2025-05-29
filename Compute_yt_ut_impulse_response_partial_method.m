%A) Signal y^2(t) generated.
    
    [y_Gcdys_impulse_partial,t_Gcdys_impulse_partial,x_Gcdys_impulse_partial]=...
        impulse(Gcdys_partial,impulse_simulation_time);
    
    y2_Gcdys_impulse_partial=y_Gcdys_impulse_partial.*y_Gcdys_impulse_partial;
    
%B) Signal u(t) generated.
    
    u_Gcdys_impulse_partial(1:size(x_Gcdys_impulse_partial,1),1:n_u,1)=...
        x_Gcdys_impulse_partial(:,:,1)*Kpf_partial';

    u2_Gcdys_impulse_partial=u_Gcdys_impulse_partial.*u_Gcdys_impulse_partial;
