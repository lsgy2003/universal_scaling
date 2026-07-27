clear all;
%%Pattern formation in phase difference and amplitude
Lx=2^11;
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.001;
T=1000;
Evo=1000;
Nt = T/dt;
ts=0:dt*Evo:T;

m= load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral/jp=0.010_2^11_2.mat');

T0=0;%transient period
phaseA_matrix=zeros(Nx,T+1-T0); %create 2D matrix of phaseA
phaseAB_matrix=zeros(Nx,T+1-T0); %create 2D matrix of phaseA-phaseB
p_matrix=zeros(Nx,T+1-T0); %create 2D matrix of |P(x,t)|^2
%f_matrix=zeros(Nx,T+1-T0); %create 2D matrix of \Omega(x,t)

for t=1:T+1-T0
    phaseA_matrix(:,t) = angle(m.Output.Solution{1,T0+t}{1,1});
    phase_dif = angle(m.Output.Solution{1,T0+t}{1,1})-angle(m.Output.Solution{1,t}{1,2});
    phaseAB_matrix(:,t) = atan2(sin(phase_dif),cos(phase_dif));
    p_matrix(:,t)= abs(m.Output.Solution{1,T0+t}{1,1});
end

clf;
%
f=figure(1);
f.Position = [300 300 800 300];
imagesc(x, T0:T, phaseA_matrix');
set(gca, 'YDir', 'normal');
ax = gca;
ax.FontSize=18;
axis([-1024 1024 400 1000])
xlabel('x');
ylabel('t');
title('Phase $\theta_A$','FontSize',20,'Interpreter','latex','FontWeight','bold');
load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral/myColor3.mat');
colormap(Custom3);
colorbar;  % To show the color scale
c=colorbar;
clim([-pi pi]);
c.Label.String='$\theta_A$';
c.Label.FontSize=18;
c.Label.Interpreter = 'latex'; 
c.Label.Rotation = 0;   % Set the label orientation (e.g., horizontal if vertical colorbar)
c.Label.Position = [0.8, -3.2, 0];  % Shift label (adjust values as needed)
set(c, 'Ticks', [-pi 0 pi], 'TickLabels', {'-\pi', '0', '\pi'});  % Define ticks and labels

%
saveas(gcf,'/Users/Phantom/Desktop/Code/data/phaseA_evo_jp=0.010_2.fig')
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/phaseA_evo_jp=0.010_2.eps','ContentType','vector','BackgroundColor','none')
%}
%
f=figure(2);
f.Position = [300 300 800 300];
imagesc(x, T0:T, log(p_matrix'));
load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral/myColor2.mat');
colormap(Custom2);
set(gca, 'YDir', 'normal');
ax = gca;
ax.FontSize=18;
axis([-1024 1024 400 1000]);
xlabel('x');
ylabel('t');
title('Amplitude of $\vec{P}_A$','FontSize',20,'Interpreter','latex','FontWeight','bold');
colorbar;
clim([-3 2]);
c=colorbar;
set(c, 'Ticks', [-2 0 2], 'TickLabels', {'10^{-2}','10^0','10^2'});  % Customize ticks for log scale, '10^{-3}',,'10^2'
c.Label.String='$|\vec{P}_A|$';
c.Label.FontSize=14;
c.Label.Rotation = 0;
c.Label.Position = [0.6, -3, 0]; 
c.Label.Interpreter = 'latex'; 

%
saveas(gcf,'/Users/Phantom/Desktop/Code/data/jm=-0.25/phaseA_amp_jp=0.010_2.fig')
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/jm=-0.25/phaseA_amp_jp=0.010_2.eps','ContentType','vector','BackgroundColor','none')
%}