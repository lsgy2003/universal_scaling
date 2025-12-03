clear all;
%temporal correlation functions
Lx=2^10;
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.01;
T=10000;
Evo=100;
Nt = T/dt;
tw=[1000 4000 7000];
Dt=1:dt*Evo:(T-tw(1))+1;
n = 240; %number of trials
jp =[0.025, 0.030, 0.040, 0.050, 0.070, 0.100, 0.150, 0.200, 0.250, 0.300, 0.400];

for l = 1:length(jp)
	
	for i = 1:length(tw)
	  ts=tw(i):dt*Evo:T;    
	  avg_corr_psi_AA = zeros(1,length(ts));
	  corr_psi_AA = zeros(1,length(ts));
	  avg_corr_psi_BB = zeros(1,length(ts));
          corr_psi_BB = zeros(1,length(ts)); %full field correlation initialization
	  
	  avg_corr_amp_AA = zeros(1,length(ts));
      corr_amp_AA = zeros(1,length(ts));
      avg_corr_amp_BB = zeros(1,length(ts));
	  corr_amp_BB = zeros(1,length(ts)); %amplitude correlation initialization
	  
	  avg_corr_phase_AA = zeros(1,length(ts));
      corr_phase_AA = zeros(1,length(ts));
      avg_corr_phase_BB = zeros(1,length(ts));
      corr_phase_BB = zeros(1,length(ts)); %phase correlation initialization
	 
	  for k=1:n
    		m=load(sprintf('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=%.3f/test/%d.mat',jp(l),k));
    		Psi_ini_A = m.Output.Solution{1,tw(i)}{1,1};
		    Psi_ini_B = m.Output.Solution{1,tw(i)}{1,2};
		
		    amp_ini_A = abs(m.Output.Solution{1,tw(i)}{1,1});
		    amp_ini_B = abs(m.Output.Solution{1,tw(i)}{1,2});
		
		    phase_ini_A = Psi_ini_A./amp_ini_A;
            phase_ini_B = Psi_ini_B./amp_ini_B;
    		
		    avg_Psi_ini_AA = sum(conj(Psi_ini_A).*Psi_ini_A)/Nx;
		    avg_Psi_ini_BB = sum(conj(Psi_ini_B).*Psi_ini_B)/Nx;
                
		    avg_amp_ini_AA = sum(amp_ini_A.^2)/Nx;
            avg_amp_ini_BB = sum(amp_ini_B.^2)/Nx;

            avg_phase_ini_AA = sum(conj(phase_ini_A).*phase_ini_A)/Nx;
            avg_phase_ini_BB = sum(conj(phase_ini_B).*phase_ini_B)/Nx;		
    		
		    for j=1:length(ts)
        		Psi_A = m.Output.Solution{1,tw(i)+j}{1,1};
			    Psi_B = m.Output.Solution{1,tw(i)+j}{1,2};
			
			    amp_A = abs(m.Output.Solution{1,tw(i)+j}{1,1});
                amp_B = abs(m.Output.Solution{1,tw(i)+j}{1,2});
                        
			    phase_A = Psi_A./amp_A;
                phase_B = Psi_B./amp_B;

        		corr_psi_AA(j) = sum(conj(Psi_ini_A).*Psi_A)/Nx; %full field correlation
			    corr_psi_BB(j) = sum(conj(Psi_ini_B).*Psi_B)/Nx;
			
			    corr_amp_AA(j) = sum(conj(amp_ini_A).*amp_A)/Nx;%amplitude correlation
                corr_amp_BB(j) = sum(conj(amp_ini_B).*amp_B)/Nx;

                corr_phase_AA(j) = sum(conj(phase_ini_A).*phase_A)/Nx;%phase correlation
                corr_phase_BB(j) = sum(conj(phase_ini_B).*phase_B)/Nx;
    		end
    		avg_corr_psi_AA = avg_corr_psi_AA + corr_psi_AA/avg_Psi_ini_AA;
		    avg_corr_psi_BB = avg_corr_psi_BB + corr_psi_BB/avg_Psi_ini_BB;
		
		    avg_corr_amp_AA = avg_corr_amp_AA + corr_amp_AA/avg_amp_ini_AA;
            avg_corr_amp_BB = avg_corr_amp_BB + corr_amp_BB/avg_amp_ini_BB;

            avg_corr_phase_AA = avg_corr_phase_AA + corr_phase_AA/avg_phase_ini_AA;
		    avg_corr_phase_BB = avg_corr_phase_BB + corr_phase_BB/avg_phase_ini_BB;
	  end
	  avg_corr_psi_AA = abs(avg_corr_psi_AA/n);
	  avg_corr_psi_BB = abs(avg_corr_psi_BB/n);

	  avg_corr_amp_AA = abs(avg_corr_amp_AA/n);
      avg_corr_amp_BB = abs(avg_corr_amp_BB/n);

	  avg_corr_phase_AA = abs(avg_corr_phase_AA/n);
      avg_corr_phase_BB = abs(avg_corr_phase_BB/n);


	  filename=sprintf('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=%.3f/test/2^10_2_tw=%d.mat',jp(l),tw(i));
	  save(filename,'avg_corr_psi_AA','avg_corr_psi_BB','avg_corr_amp_AA','avg_corr_amp_BB','avg_corr_phase_AA','avg_corr_phase_BB');
	end
end
