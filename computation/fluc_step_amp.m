function state = fluc_step_amp(state, B_init, avg_B_init, B_now, A_init, avg_A_init, A_now)

    %% --- C_AA: temporal correlation of field_A ---
    avg_A_now = sqrt(mean(conj(A_now).* A_now));
    corr_field_AA = mean(conj(A_init).* A_now);

    state.C_AA (state.save_counter)  = abs(corr_field_AA)/avg_A_now/avg_A_init;

    %% --- C_BB: temporal correlation of field_B ---
    avg_B_now = sqrt(mean(conj(B_now).* B_now));
    corr_field_BB = mean(conj(B_init).* B_now);

    state.C_BB (state.save_counter)  = abs(corr_field_BB)/avg_B_now/avg_B_init;

    %% --- amplitude overlap for field_A ---

    amp_corr_AA = mean(abs(A_init).*abs(A_now));

    state.amp_AA(state.save_counter) = amp_corr_AA /avg_A_now / avg_A_init;


    %% --- phase-only correlation for field_A ---

    epsA = 1e-12;

    A_init_hat = A_init ./ max(abs(A_init), epsA);

    A_now_hat  = A_now  ./ max(abs(A_now),  epsA);

    phase_corr_AA = mean(conj(A_init_hat).*A_now_hat);

    state.phase_AA(state.save_counter) = abs(phase_corr_AA);

    %% -----Update time step-----
    state.save_counter = state.save_counter + 1;
 
end