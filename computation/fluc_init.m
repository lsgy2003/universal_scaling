function state = fluc_init(Nt_settle)

    %temporal correlation C_AA
    state.C_AA = zeros(Nt_settle,1);

    state.save_counter = 1;
    
end
