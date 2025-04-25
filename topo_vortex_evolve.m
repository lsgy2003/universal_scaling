clear all;
%%Pattern formation in phase difference and amplitude
Lx=2^8;
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.01;
T=800;
Evo=100;
Nt = T/dt;
ts=0:dt*Evo:T;

m= load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/lin_stability/2^8_phase_dif_noise6.mat');

T0=0;%transient period
phaseA_matrix=zeros(Nx,T+1-T0); %create 2D matrix of phaseA
phaseAB_matrix=zeros(Nx,T+1-T0); %create 2D matrix of phaseA-phaseB
p_matrix=zeros(Nx,T+1-T0); %create 2D matrix of |P(x,t)|^2
f_matrix=zeros(Nx,T+1-T0); %create 2D matrix of \Omega(x,t)

for t=1:T+1-T0
    phaseA_matrix(:,t) = angle(m.Output.Solution{1,T0+t}{1,1});
    phase_dif = angle(m.Output.Solution{1,T0+t}{1,1})-angle(m.Output.Solution{1,t}{1,2});
    phaseAB_matrix(:,t) = atan2(sin(phase_dif),cos(phase_dif));
    p_matrix(:,t)= abs(m.Output.Solution{1,T0+t}{1,1});
    f_matrix(:,t)=m.Output.Freq{1,T0+t}{1,1};      
end

clf;
%
f=figure(1);
f.Position = [300 300 800 300];
imagesc(x, T0:T, phaseAB_matrix');
set(gca, 'YDir', 'normal');
ax = gca;
ax.FontSize=18;
xlabel('x');
ylabel('t');
title('Phase $\Delta\theta$','FontSize',20,'Interpreter','latex','FontWeight','bold');
colormap(parula);
cmap = colormap;  % Get the hsv colormap with 256 colors
cmap = cmap * 0.5;  % Reduce the saturation/brightness for a softer appearance
colorbar;  % To show the color scale
c=colorbar;
clim([-pi pi]);
c.Label.String='\Delta\theta';
c.Label.FontSize=18;
c.Label.Rotation = 0;   % Set the label orientation (e.g., horizontal if vertical colorbar)
c.Label.Position = [0.8, -3.2, 0];  % Shift label (adjust values as needed)
set(c, 'Ticks', [-pi 0 pi], 'TickLabels', {'-\pi', '0', '\pi'});  % Define ticks and labels
saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_phase_noise6.fig')
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_phase_noise6.eps','ContentType','vector','BackgroundColor','none')
%}
%
f=figure(2);
f.Position = [300 300 800 300];
imagesc(x, T0:T, log(p_matrix'));
set(gca, 'YDir', 'normal');
ax = gca;
ax.FontSize=18;

xlabel('x');
ylabel('t');
title('Amplitude of $\vec{P}_A$','FontSize',20,'Interpreter','latex','FontWeight','bold');
colorbar;
caxis([-3 2]);
c=colorbar;
set(c, 'Ticks', [-2 0 2], 'TickLabels', {'10^{-2}','10^0','10^2'});  % Customize ticks for log scale, '10^{-3}',,'10^2'
c.Label.String='$|\vec{P}_A|$';
c.Label.FontSize=14;
c.Label.Rotation = 0;
c.Label.Position = [0.6, -3, 0]; 
c.Label.Interpreter = 'latex'; 
%c.Ticks=[0 0.1 0.2 0.3 0.4];
saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_amp_noise6.fig')
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_amp_noise6.eps','ContentType','vector','BackgroundColor','none')
%}

%{
subplot(2,1,3)
imagesc(x, T0:T, phaseAB_matrix');
set(gca, 'YDir', 'normal');
ax = gca;
ax.FontSize=12;
colorbar;  % To show the color scale
xlabel('x');
ylabel('t');
title('Phase Difference of P_A and P_B','FontSize',10);
c=colorbar;
c.Label.String='\theta_A-\theta_B';
c.Label.FontSize=12;
c.Ticks=[-3 -2 -1 0 1 2 3];

subplot(2,2,4)
imagesc(x, T0:T, f_matrix');
set(gca, 'YDir', 'normal');
ax = gca;
ax.FontSize=12;
colorbar;  % To show the color scale
xlabel('x');
ylabel('t');
title('Frequency of P_A','FontSize',10);
c=colorbar;
c.Label.String='\Omega_A';
c.Label.FontSize=12;
c.Ticks=[-2 -1 0 1 2];
%}
 