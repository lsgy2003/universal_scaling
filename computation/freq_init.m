function state = freq_init(Nx, Nt)
    
    state.Sk_A  = zeros(1,Nx);
    
    %frequency & its spatial variance % width of phase
    state.freq_A = zeros(Nt,1);
    state.Var_freq = zeros(1, Nx);
    state.counter = 1;
    state.width_A = zeros(Nt,1);
    
end
