function result = fluc_finalize(state, Nt_settle)

    %Calculate the power spectral density
    C_AA_temp = state.C_AA - mean(state.C_AA);
    %[PSD,f] = pwelch(C_AA_temp, [], [], [], fs);   % default Welch

    YA = fft(C_AA_temp); 
    P2 = abs(YA/Nt_settle);
    P1 = P2(1:floor(Nt_settle/2+1)); %single-sided spectrum
    P1(2:end-1) = 1*P1(2:end-1);

    %result.f = f;
    result.P_AA = P1;
    result.C_AA = state.C_AA;
end
