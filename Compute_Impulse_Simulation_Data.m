%B) Signal y(t), x(t), and the 't' associated with impulse d(t)

[y_Canon_Godys_impulse,t_Canon_Godys_impulse,x_Canon_Godys_impulse]=...
	impulse(Canon_Godys,impulse_simulation_time);
    
%C) Signal y(t)^2 associated with impulse d(t)
    
    y2_Canon_Godys_impulse=y_Canon_Godys_impulse.*y_Canon_Godys_impulse;
