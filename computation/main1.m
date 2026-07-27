clear all; 

%% ---------------- Geometry ----------------
Lx = 2^10;
dx = 1/2;
x  = (-Lx/2:dx:Lx/2)';
Nx = Lx/dx + 1;

iLeft   = [Nx, 1:Nx-1]';
iCenter = (1:Nx)';
iRight  = [2:Nx, 1]';

%% ---------------- Time control ----------------
dt = 0.001;
Evo = 1000;
T = 10000; %total time
T_settle = 7000; %relaxation time
Nt_total = T/dt + 1;
dt_save = dt*Evo;
Nt = T/dt_save + 1; %total number of output slices 
Nt_settle = (T- T_settle)/dt_save + 1;
fs = 1/dt_save; %sampling frequency in power spectrum
%% ---------------- Sweep ----------------

DBB = 1;
DAA = 100;

%% ---------------- Physical parameters ----------------
jp=0.0020;

ini = load('/Users/Phantom/Desktop/Code/data/jm=-0.25_ini.mat'); %initial are mean-field steady states
jp_all = ini.jp_all;
A0 = ini.A(jp_all == jp);
B0 = ini.B(jp_all == jp);

jAA = 1; jBB = 1;
jm  = -0.25;
jAB = jp + jm;
jBA = jp - jm;


rhoA = 1; rhoB = 1;

eta_c = 1/2*(jAA + jBB)*rhoA;
eta   = 0.5*eta_c;

xi_A = 0.002; xi_B = 0.002;
Eta_A = xi_A / sqrt(dt*dx);
Eta_B = xi_B / sqrt(dt*dx);


%% Trials
nTrials = 96;   % set as needed

results = cell(nTrials,1);


%% Output root
outroot = '/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/';
if ~exist(outroot,'dir'); mkdir(outroot); end

outfile  = fullfile(outroot, sprintf('xi=%.3f_jp=%.4f.mat', xi_A, jp));


parfor trial = 1:nTrials

    % Record Cputime
    Cputime_ini=cputime;

    % Each worker uses its own RNG stream
    rng(trial,'twister');

    % ---------- Initialize ONCE per trial ----------
    A_now = A0*ones(Nx,1);
    B_now = B0*ones(Nx,1);

    % --- initialize analysis states ---
    freq_state   = freq_init(Nx, Nt);
    fluc_state = fluc_init(Nt_settle);

    A_init = [];
    avg_A_init = NaN;
    have_init = false;

    A_prev_save = [];
    B_prev_save = [];
    have_prev_save = false;

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

        NL_AA = -abs(jAA*A + jAB*B).^2/(2*eta);
        NL_BB = -abs(jBA*A + jBB*B).^2/(2*eta);

        dAdt = Coupled_AA*A + NL_AA.*A + Coupled_AB*B + DAA*del2A ...
            + Eta_A*(randn(Nx,1) + 1i*randn(Nx,1));
        dBdt = Coupled_BB*B + NL_BB.*B + Coupled_BA*A + DBB*del2B ...
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

            % ---------- Frequency -----------
            freq_state = freq_step(freq_state, freqA, phase_unw);

            
            % ---------- Temporal Correlation -----------
            if ~have_init && it >= T_settle/dt
                 % capture the first frame at/after settle time
                 A_init = A_now;
                 avg_A_init = sqrt(mean(abs(A_init).^2));   % same as /Nx
                 have_init = true;
            end

            if have_init && it > T_settle/dt
                fluc_state = fluc_step(fluc_state, A_init, avg_A_init, A_now);
            end
            %
        end

        % shift
        A_now = A_next;
        B_now = B_next;
        
    end

    % --- finalize ---
    freq_result   = freq_finalize(freq_state, Nt, Nx);
    fluc_result = fluc_finalize(fluc_state, fs);

    %----accumulate per trial----
    trial_result = struct();
    trial_result.freq      = freq_result;
    trial_result.fluc    = fluc_result;

    results{trial} = trial_result;
    
    Cputime = cputime - Cputime_ini;
    fprintf('CPU time %d:  %8.2f\n',trial,Cputime);
    
end
%% ----- Ensembel Average --------
R = [results{:}];   % struct array

avg_freq      = ensemble_average({R.freq});
avg_fluc    = ensemble_average({R.fluc});

%save averaged results
save(outfile, ...
    'avg_fluc', 'avg_freq',  ...
    'Lx', 'jp', 'nTrials', ...
    '-v7.3'); 
