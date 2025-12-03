clear all;
clf;

%%Setting the geometry
Lx=2^8;
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.01;
T=800;
Evo=100;
Nt = T/dt;
t=0:dt:T;
ts=0:dt*Evo:T;

%% Setting the physical problem

jp=0.002;
ini=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/uni_ini.mat'); %uniform steady states computed from mean-field theory
jp_all=ini.jp_all;
k=find(jp_all==jp);
A0=ini.A(k);
B0=ini.B(k);

jAA=0.5;
jBB=0.5;
jm=-0.25;
jAB=jp+jm;
jBA=jp-jm;

DA=1;
DB=1;

%Stochastic Noise Amplitude
xi_A=0.1; %\sigma
xi_B=0.1;
Eta_A=xi_A/sqrt(dt*dx); % Noise amplitude
Eta_B=xi_B/sqrt(dt*dx);

a=1; %Gaussian perturbation amplitude
sig=10; %Gaussian perturbation width
%
filename2='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/lin_stability/2^8_phase_dif_noise_Var.mat';
save(filename2,'jp','jm','xi_A','DA','rhoA','Lx','T','a','sig');

%
u= VideoWriter('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/topo_vortex_noise3','MPEG-4');
open(u); 


n=240;
parfor i=1:n
%% Zero dim linear equation

iLeft=[Nx, 1:Nx-1]'; %perodic boundary condition
iCenter = (1:Nx)';
iRight = [2:Nx, 1]';


 Cputime_ini=cputime;

 %
 filename1='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/lin_stability/2^8_phase_dif_noise.mat';
 m=matfile(filename1,'Writable',true);
 %}


    %%Initialization of variables    
    phase=cell(1,2); %winding phases
    amp=cell(1,2);
    phase1=cell(1,2);
    phase2=cell(1,2);
    phase3=cell(1,2); %unwinding phases
    phase4=cell(1,2); %unwinding phases in the out-of-phase basis
    phase5=cell(1,2); %for the calculation of the frequency    
    corr_phase=cell(1,2);
    corr_phase1=cell(1,2); %correlation of the unwinded phase
    density=cell(1,2);
    %Output variables
    Output = struct();
    Output.corr_phase=cell(1,2);
    Output.density =cell(1,2);
    Output.freq =cell(1,2);
    Output.Solution=cell(1,2);
    Output.Phase = cell(1,2); %Output phase
    Indicator=cell(1,2);
    down=cell(1,2);
    up=cell(1,2);
    for n = 1:2
        Indicator{n}=zeros(Nx,1);
        down{n}=zeros(Nx,1);
        up{n}=zeros(Nx,1);
    end
    freq=cell(1,2);
    site_freq=cell(1,2);

    N_steps=0;

    %%Initial Condition
        % Computing a random normal matrix
    rng(n,'twister');
    A = a*(randn(Nx,1) + 1i*rand(Nx,1));  %random initial
    B = a*(randn(Nx,1) + 1i*rand(Nx,1));
    %A =ones(Nx,1)*A0; %uniform initial
    %B =ones(Nx,1)*B0;  
    %A = A0*ones(Nx,1)+(a*exp(-x.^2/sig)); %uniform initial + a small real Gaussian perturbation
    %B = B0*ones(Nx,1)+(a*exp(-x.^2/sig));


    dAdt=zeros(Nx,1);
    dBdt=zeros(Nx,1);


    for j=1:length(t)

        %finite difference method + perodic boundary conditions
        delA= (A(iRight)-A(iLeft))/(2*dx);
        delB= (B(iRight)-B(iLeft))/(2*dx);
        del2A= (A(iLeft)+A(iRight)-2*A(iCenter))/dx^2;
        del2B= (B(iLeft)+B(iRight)-2*B(iCenter))/dx^2;

        %%Setting the Physical Conditions
        %Potential
        Coupled_AA =  jAA;
        Coupled_AB =  jAB;
        Coupled_BA =  jBA;
        Coupled_BB =  jBB;

        %Nonlinearity
        NL_AA = -abs(jAA*A+jAB*B).^2/(2*eta);
        NL_BB = -abs(jBA*A+jBB*B).^2/(2*eta);

        %Full equations
        dAdt= Coupled_AA*A+NL_AA.*A+Coupled_AB*B+DA*del2A+Eta_A*(randn(Nx,1)+1i*randn(Nx,1));
        dBdt= Coupled_BB*B+NL_BB.*B+Coupled_BA*A+DB*del2B+Eta_B*(randn(Nx,1)+1i*randn(Nx,1)); %with Nonlinearity term

        A=A+dAdt*dt;
        B=B+dBdt*dt;
        phase{1} = angle(A);
        phase{2} = angle(B);
        amp{1} = abs(A);
        amp{2} = abs(B);

        %%Unwind the phase in time and space
        for n = 1:2
            % initialize the first time step
            if (j==1)
                phase1{n} = phase{n};
                phase3{n} = phase1{n};
                phase5{n} = phase{n};
                corr_phase{n} = sum((phase3{n}-sum(phase3{n})/Nx).^2)/Nx; %Calculate the autocorrelation of the phase (width)
                density{n} = sum(amp{n}.^2)/Nx;

            % for all the following time steps
            elseif (j > 1)
                phase2{n} = phase{n};
            %Modify the phase in time
                down{n} = (phase2{n}-phase1{n}<-pi);
                up{n} = (phase2{n}-phase1{n}> pi);
                Indicator{n}=Indicator{n}+down{n}-up{n};
                phase3{n}=phase2{n}+2*pi*Indicator{n}; %phase in time
            %Calculate the autocorrelation of the phase (width)
                corr_phase{n} = sum((phase3{n}-sum(phase3{n})/Nx).^2)/Nx;
            %Calculate the density of the profile
                density{n} = sum(amp{n}.^2)/Nx;
            %Repeat the iteration
                phase1{n} = phase2{n};
            end
        end
                                                                                                                                                                                          
        %phase_p = (phase3{1}+phase3{2})/sqrt(2); % in-phase mode
        %phase_m = (phase3{1}-phase3{2})/sqrt(2); % out-of-phase mode
        %corr_phase{1} = sum((phase_p-sum(phase_p)/Nx).^2)/Nx;
        %corr_phase{2} = sum((phase_m-sum(phase_m)/Nx).^2)/Nx;
        phase_dif = phase{1}-phase{2}; %phase difference
        phase_dif_winding = atan2(sin(phase_dif),cos(phase_dif));


        %}
        % Set output variables
        if mod(j-1,Evo)==0
            N_steps = N_steps+1;
            Output.Solution{N_steps}{1} = A; 
            Output.Solution{N_steps}{2} = B;
           for n=1:2
            freq{n}=(phase3{n}-phase5{n})/(dt*Evo);%calculate the angular frequency;
            phase5{n}=phase3{n};
            Output.corr_phase{n}((j-1)/Evo+1) = corr_phase{n};
            Output.density{n}((j-1)/Evo+1) = density{n};
            Output.freq{n}((j-1)/Evo+1) = sum(abs(freq{n}))/Nx;
            Output.phase_dif{n}((j-1)/Evo+1)= sum(abs(phase_dif_winding))/Nx;
            Output.Phase{N_steps}{n} = phase3{n};
            Output.Freq{N_steps}{n} = freq{n};
           end
            
           %%Plot 
           f = figure(1);
           f.Position = [300 300 800 900];
           
           subplot(3,1,1)
           plot(x,freq{1},x,freq{2},'LineWidth',1); %x,freq{2}
           axis([-Lx/2 Lx/2 -2 2]);
           xlabel('x')
           ylabel('\Omega_A(x), \Omega_B(x)')
           title('Frequency $\Omega_{A}$ and $\Omega_{B}$','Interpreter','latex','FontSize',16)
           legend('\Omega_{A}','\Omega_{B}','Location','southeast','FontSize',18);
           ax = gca;
           ax.FontSize=18;
           drawnow;
           %
           subplot(3,1,2)
           plot(x,abs(A),x,abs(B),'LineWidth',1);
           axis([-Lx/2 Lx/2 0 1]);
           title('Amplitude $|\vec{P}_A|$ and $|\vec{P}_B|$','Interpreter','latex','FontSize',16)
           xlabel('x')
           ylabel('|P_A(x)|, |P_B(x)|')
           legend('|P_A|','|P_B|','Location','southeast','FontSize',18);
           ax = gca;
           ax.FontSize=18;
           drawnow;
            
           %
           subplot(3,1,3)
           phase_dif = phase{1}-phase{2};
           phase_dif_winding = atan2(sin(phase_dif),cos(phase_dif));

           plot(x,phase_dif_winding,'Color','k','LineWidth',1.5);
           hold off;
           axis([-Lx/2 Lx/2 -2*pi 2*pi]);
           title('Phase $\Delta \theta = \theta_{A} - \theta_{B}$','Interpreter','latex','FontSize',16)
           xlabel('x')
           ylabel('\Delta \theta(x)')
           ax = gca;
           ax.FontSize=18;
           drawnow;
           %}
         

           figure(2)
           plot(x,freq{1},'Color',[0 0 0]/255,'LineWidth',2); %x,freq{2}
           set(gca,'color', [211 211 211]/255);
           axis([-Lx/2 Lx/2 -1 1]);
           xlabel('x')
           ylabel('\Omega_A')
           ax = gca;
           ax.FontSize=18;
           drawnow;
 
           fprintf('Time Step %d\n',(j-1)/Evo); 
           %}

         
           frame=getframe(gcf); 
           writeVideo(u,frame);  
             
 
        end
                    
    end

    Cputime = cputime - Cputime_ini;
    fprintf('CPU time %d:  %8.2f\n',i,Cputime);
    m.Output=Output;
end
close(u)