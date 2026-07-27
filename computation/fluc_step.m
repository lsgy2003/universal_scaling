function state = fluc_step(state, A_init, avg_A_init, A_now)

    %% --- C_AA: temporal correlation of field_A ---
    avg_A_now = sqrt(mean(conj(A_now).* A_now));
    corr_field_AA = mean(conj(A_init).* A_now);

    state.C_AA (state.save_counter)  = abs(corr_field_AA)/avg_A_now/avg_A_init;

    state.save_counter = state.save_counter + 1;

    
 
end