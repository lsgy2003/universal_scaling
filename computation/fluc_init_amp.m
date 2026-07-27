function state = fluc_init_amp(Nt_settle)

    state.C_AA = zeros(Nt_settle,1);
    state.C_BB = zeros(Nt_settle,1);
    state.amp_AA = zeros(Nt_settle,1);
    state.phase_AA = zeros(Nt_settle,1);

    state.save_counter = 1;
    
end
