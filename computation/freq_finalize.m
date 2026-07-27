function result = freq_finalize(state, Nx, Nt)

    %time average
    Sk = state.Sk_A / Nt / Nx; %Nx is the normalization factor of FFT
    Var   = state.Var_freq / Nt;  

    % shift zero mode to center
    Sk = fftshift(Sk);

    result.Sk_A = Sk;
    result.freq_A = state.freq_A;
    result.Var_freq = Var;
    result.width_A = state.width_A;

end
