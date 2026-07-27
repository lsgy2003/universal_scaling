function state = freq_step(state,freqA,phase_unw)

    % --- average frequency ---
    state.freq_A (state.counter) = mean(abs(freqA));

    uhat_A = fft(freqA - mean(freqA))';
    
    state.Sk_A = state.Sk_A + abs(uhat_A);

    % -- Var(Omega): spatial variance of frequency ---

    state.Var_freq = state.Var_freq + ((freqA - mean(freqA)).^2)';

    % -- width of phase ---
     state.width_A(state.counter) = mean((phase_unw{1} - mean(phase_unw{1})).^2);

    %% Time step update
     state.counter = state.counter + 1;

end