function result = fluc_finalize_amp(state, fs)

    %Calculate the power spectral density
    C_AA_temp = state.C_AA - mean(state.C_AA);
    [PSD,f] = pwelch(C_AA_temp, [], [], [], fs);   % default Welch

    result.f = f;
    result.P_AA = PSD;
    result.C_AA = state.C_AA;
    result.C_BB = state.C_BB;
    result.amp_AA = state.amp_AA;
    result.phase_AA = state.phase_AA;
end
