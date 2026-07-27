function state = fluc_init2(Nt_settle, Nx)

    %temporal correlations
    state.C_AA = zeros(Nt_settle,1);
    state.C1_AA = zeros(Nt_settle,1);
    state.C2_AA = zeros(Nt_settle,Nx);
    state.phase_AA = zeros(Nt_settle,1);
    state.phase1_AA = zeros(Nt_settle,1);
    state.phase2_AA = zeros(Nt_settle,Nx);

    %Var(theta)
    state.Var1 = zeros(Nt_settle,1);
    state.Var2 = zeros(Nt_settle,1);
    state.Var3 = zeros(Nt_settle,1);
    state.Var4  = zeros(Nt_settle,Nx);
    state.Var5  = zeros(Nt_settle,Nx);

    state.save_counter = 1;
    
end
