function state = fluc_step2(state,A_init, avg_A_init, A_now, phase_unw_ini_A, phase_unw)

%% --- C_AA: temporal correlation of field_A ---
    avg_A_now = sqrt(mean(conj(A_now).* A_now));

    corr_field_AA = mean(conj(A_init).* A_now);

    state.C_AA (state.save_counter)  = abs(corr_field_AA)/avg_A_now/avg_A_init;   

  % --- phase-only correlation for field_A ---

    epsA = 1e-12;

    A_init_hat = A_init ./ max(abs(A_init), epsA);

    A_now_hat  = A_now  ./ max(abs(A_now),  epsA);

    phase_corr_AA = mean(conj(A_init_hat).*A_now_hat);

    state.phase_AA(state.save_counter) = abs(phase_corr_AA);

    %% --- C'_AA ---

    state.C1_AA (state.save_counter)  = corr_field_AA/avg_A_now/avg_A_init;

    % --- phase-only correlation for field_A ---

    state.phase1_AA(state.save_counter) = phase_corr_AA;


   %% --- C''_AA ---
    corr2_field_AA = conj(A_init).* A_now;

    state.C2_AA (state.save_counter,:)  = (corr2_field_AA/avg_A_now/avg_A_init)';
  
    % --- phase-only correlation for field_A ---
    phase2_corr_AA = conj(A_init_hat).*A_now_hat;

    state.phase2_AA(state.save_counter,:) = phase2_corr_AA;

    
    %% -- Building blocks of Var(theta): spatial variance of unwinding phase ---
    Delta_phase_A = phase_unw{1} - phase_unw_ini_A;

    %size(Nt,1)
    state.Var1(state.save_counter) = mean(Delta_phase_A.^2); 
    state.Var2(state.save_counter) = mean(Delta_phase_A).^2;
    state.Var3(state.save_counter) = mean(Delta_phase_A);
    %size(Nt,Nx)
    state.Var4(state.save_counter,:) = Delta_phase_A.^2'; 
    state.Var5(state.save_counter,:) = Delta_phase_A';

                               
    %% -----Update time step-----
    state.save_counter = state.save_counter + 1;
 
end