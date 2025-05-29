%A) Total Controllability Matrix
    
    total_controllability_matrix=ctrb(A_x,B_u);

    %B) Number of uncontrollable states

    estados_nao_controlaveis_sistema_original=...
        n_x-rank(total_controllability_matrix);
