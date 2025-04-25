clear all;

sigma=0.5;
jp=0.002;

T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
ts=0:dt*Evo:T;
Lx=[2^8 2^9 2^10 2^11 2^12 2^13 2^14]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^8.mat');
freq0=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^9.mat');
freq1=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^10.mat');
freq2=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^11.mat');
freq3=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^12.mat');
freq4=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^13.mat');
freq5=m.omega1(:);

%m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^14.mat');
%freq6=m.omega1(:);

clf;
figure(1)
loglog(ts,freq0,'Color',[100, 149, 237]/255);
hold on;
loglog(ts,freq1,'Color',[255, 191, 0]/255);
hold on;
loglog(ts,freq2,'Color',[255, 127, 80]/255);
hold on;
loglog(ts,freq3,'Color',[204, 204, 255]/255);
hold on;
loglog(ts,freq4,'Color',[159, 226, 191]/255);
hold on;
loglog(ts,freq5,'Color',[0 0 0]);
hold on;
xlabel('$t$','Interpreter','latex')
ylabel('$\overline{\Omega}_{A}(L,t)$','Interpreter','latex')
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','L=2^{13}','Location','southwest'); %,'L=2^{14}'
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;


saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_L_sigma=0.5.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_L_sigma=0.5.eps','ContentType','vector','BackgroundColor','none');





