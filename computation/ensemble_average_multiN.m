function avg = ensemble_average_multiN(results, N_list, save_raw_trials)

if nargin < 3
    save_raw_trials = false;
end

nN = numel(N_list);

for k = 1:nN
    N = N_list(k);

    % take first N trials
    tmp = [results{1:N}];

    %% --- ensemble averages of raw building blocks ---
    avg_fluc.C_AA = mean(cat(2,tmp.C_AA),2);
    avg_fluc.C1_AA     = mean(cat(2,tmp.C1_AA),2);
    avg_fluc.C2_AA     = mean(cat(3,tmp.C2_AA),3);

    avg_fluc.phase_AA = mean(cat(2,tmp.phase_AA),2);
    avg_fluc.phase1_AA = mean(cat(2,tmp.phase1_AA),2);
    avg_fluc.phase2_AA = mean(cat(3,tmp.phase2_AA),3);

    avg_fluc.Var1 = mean(cat(2,tmp.Var1),2);
    avg_fluc.Var2 = mean(cat(2,tmp.Var2),2);
    avg_fluc.Var3 = mean(cat(2,tmp.Var3),2); % Var1-3, Nt x N

    avg_fluc.Var4 = mean(cat(3,tmp.Var4),3);   % Var4-5, Nt x Nx x N
    avg_fluc.Var5 = mean(cat(3,tmp.Var5),3);   

    %% --- Calculate desired quantities ---
    avg.C_AA(:,k)     = avg_fluc.C_AA;
    avg.C1_AA(:,k) = abs(avg_fluc.C1_AA);
    avg.C2_AA(:,k) = mean(abs(avg_fluc.C2_AA),2);
    
    avg.phase_AA(:,k) = avg_fluc.phase_AA;
    avg.phase1_AA(:,k) = abs(avg_fluc.phase1_AA);
    avg.phase2_AA(:,k) = mean(abs(avg_fluc.phase2_AA),2);

    avg.Var_A(:,k)  = avg_fluc.Var1 - avg_fluc.Var2;
    avg.Var1_A(:,k) = avg_fluc.Var1 - avg_fluc.Var3.^2;
    avg.Var2_A(:,k) = mean(avg_fluc.Var4,2) - mean(avg_fluc.Var5.^2,2);

    %% --- optionally save raw trial data ---
    if save_raw_trials
        avg.raw{k} = tmp;
        avg.raw_N(k) = N;
    end
end

avg.N_list = N_list;

end