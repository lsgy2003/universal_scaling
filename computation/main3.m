clear all; 

%% ---------------- Geometry ----------------
Lx = 2^8;
dx = 1/2;
x  = (-Lx/2:dx:Lx/2)';
Nx = Lx/dx + 1;

iLeft   = [Nx, 1:Nx-1]';
iCenter = (1:Nx)';
iRight  = [2:Nx, 1]';

%% ---------------- Time control ----------------
dt = 0.001;
Evo = 1000;
T = 100; %total time
T_settle = 20; %relaxation time
Nt_total = T/dt + 1;
dt_save = dt*Evo;
Nt = T/dt_save + 1; %total number of output slices 
Nt_settle = (T- T_settle)/dt_save + 1;
fs = 1/dt_save; %sampling frequency in power spectrum
%% ---------------- Sweep ----------------

DBB = 1;
DAA = 100;
DAB = 0;
DBA = 0;

%% ---------------- Physical parameters ----------------
jp=0.0097;

ini = load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/data/jm=-0.25.mat');
jp_all = ini.jp_all;
idx = find(abs(jp_all - jp) < 1e-12, 1);
A0 = ini.A_ini(idx);
B0 = ini.B_ini(idx);

jAA = 1; jBB = 1;
jm  = -0.25;
jAB = jp + jm;
jBA = jp - jm;
gA = 0;
gB = 0;


rhoA = 1; rhoB = 1;

eta_c = 1/2*(jAA + jBB)*rhoA;
eta   = 0.5*eta_c;

xi_A = 0.01; xi_B = 0.01;
Eta_A = xi_A / sqrt(dt*dx);
Eta_B = xi_B / sqrt(dt*dx);


%% Output root
outroot = '/Users/Phantom/Desktop/Code/test';
if ~exist(outroot,'dir'); mkdir(outroot); end

outfile  = fullfile(outroot, sprintf('jp=%.4f.mat', jp));


%% Trials
nTrials = 2^3;   % maximum ensemble size for one-run N scan
N_list  = 2.^(1:3);
save_raw_trials = false; % set true only if you want to save all 1024 trials

results = cell(nTrials,1);

parfor trial = 1:nTrials
    % Record Cputime
    Cputime_ini=cputime;

    % Each worker uses its own RNG stream
    rng(trial,'twister');

    % ---------- Initialize ONCE per trial ----------
    A_now = A0*ones(Nx,1);
    B_now = B0*ones(Nx,1);
   
    % --- initialize analysis states ---
    fluc_state = fluc_init2(Nt_settle, Nx);

    A_init = [];
    avg_A_init = NaN;
    phase_init_A = [];
    have_init = false;

    A_prev_save = [];
    B_prev_save = [];
    have_prev_save = false;

    phase_unw_ini_A = [];

    % ---- unwrap bookkeeping RESET per DA ----iii
    phase_prev = cell(1,2);
    phase_unw  = cell(1,2);
    phase_ref  = cell(1,2);
    Indicator  = {zeros(Nx,1), zeros(Nx,1)};

    for it = 1:Nt_total
        tnow = (it-1)*dt;

        A = A_now;
        B = B_now;
        
        delA= (A(iRight)-A(iLeft))/(2*dx);
        delB= (B(iRight)-B(iLeft))/(2*dx);
        del2A = (A(iLeft) + A(iRight) - 2*A(iCenter)) / dx^2;
        del2B = (B(iLeft) + B(iRight) - 2*B(iCenter)) / dx^2;

        Coupled_AA = -(eta - jAA*rhoA);
        Coupled_AB =  (jAB*rhoA);
        Coupled_BA =  (jBA*rhoB);
        Coupled_BB = -(eta - jBB*rhoB);

        NL_AA = -abs(jAA*A + jAB*B).^2/(2*eta) - gA*(conj(delA).*delA);
        NL_BB = -abs(jBA*A + jBB*B).^2/(2*eta) - gB*(conj(delB).*delB);

        dAdt = Coupled_AA*A + NL_AA.*A + Coupled_AB*B + DAA*del2A + DAB*del2B ...
            + Eta_A*(randn(Nx,1) + 1i*randn(Nx,1));
        dBdt = Coupled_BB*B + NL_BB.*B + Coupled_BA*A + DBB*del2B + DBA*del2A ...
            + Eta_B*(randn(Nx,1) + 1i*randn(Nx,1));

        A = A + dAdt*dt;
        B = B + dBdt*dt;

        A_next = A;
        B_next = B;

        phase_now = {angle(A), angle(B)};

        for n = 1:2
            if it == 1
                phase_prev{n} = phase_now{n};
                phase_unw{n}  = phase_now{n};
                phase_ref{n}  = phase_unw{n};
            else
                dp = phase_now{n} - phase_prev{n};
                Indicator{n} = Indicator{n} + (dp < -pi) - (dp > pi);
                phase_unw{n} = phase_now{n} + 2*pi*Indicator{n};
                phase_prev{n} = phase_now{n};
            end
        end
        % ===================== Save Functions ====================
        
        if mod(it-1, Evo) == 0    
            freqA = (phase_unw{1} - phase_ref{1}) / (dt*Evo);
            freqB = (phase_unw{2} - phase_ref{2}) / (dt*Evo);
            phase_ref{1} = phase_unw{1};
            phase_ref{2} = phase_unw{2};
   
            % current saved snapshot (at time tnow)
            A_save = A_now;
            B_save = B_now;
          
            %
            % ---------- Temporal Correlation -----------
            if ~have_init && it >= T_settle/dt
                 % capture the first saved frame at/after settle time
                 % Use A_next because phase_unw was just updated from A_next.
                 A_init = A_next;
                 avg_A_init = sqrt(mean(abs(A_init).^2));   % RMS amplitude
                 phase_unw_ini_A = phase_unw{1};               % fixed unwrapped phase reference
                 have_init = true;
            end

            if have_init && it > T_settle/dt
                fluc_state = fluc_step2(fluc_state,A_init, avg_A_init, A_now, phase_unw_ini_A, phase_unw);
            end
            %}
        end

        % shift
        A_now = A_next;
        B_now = B_next;
        
    end

    % --- finalize ---
    fluc_result = fluc_finalize2(fluc_state);

    %----accumulate per trial----
    trial_result = struct();
    trial_result.fluc    = fluc_result;

    results{trial} = trial_result;
    
    Cputime = cputime - Cputime_ini;
    fprintf('CPU time %d:  %8.2f\n',trial,Cputime);
    
end
%% ----- Ensembel Average --------
R = [results{:}];   % struct array

avg_fluc_N  = ensemble_average_multiN({R.fluc}, N_list, save_raw_trials);

%% -- save results --
save(outfile, 'avg_fluc_N', 'N_list', 'Lx', 'jp', 'nTrials','-v7.3');