clear all;

jp=0.002;

T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
ts=0:dt*Evo:T;
sigma=[0.001 0.01 0.1 0.2 0.3 0.4 0.5];


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.001/Size/2^13.mat');
freq0=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.01/Size/2^13.mat');
freq1=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.1/Size/2^13.mat');
freq2=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.2/Size/2^13.mat');
freq3=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.3/Size/2^13.mat');
freq4=m.omega1(:);
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.4/Size/2^13.mat');
freq5=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^13.mat');
freq6=m.omega1(:);
%}
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
loglog(ts,freq5,'Color',[0 0.4470 0.7410]);
hold on;
loglog(ts,freq6,'Color',[0 0 0]);
hold on;
xlabel('$t$','Interpreter','latex')
ylabel('$\overline{\Omega}_{A}(\sigma,t)$','Interpreter','latex')
lgd=legend('\sigma=0.001','\sigma=0.01','\sigma=0.1','\sigma=0.2','\sigma=0.3','\sigma=0.4','\sigma=0.5','Location','southwest'); %,'L=2^{14}'
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;

saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_xi_L=2^13_1.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_xi_L=2^13_1.eps','ContentType','vector','BackgroundColor','none');


figure(2)
avg_freq0=mean(freq0(1000:end));
avg_freq1=mean(freq1(1000:end));
avg_freq2=mean(freq2(1000:end));
avg_freq3=mean(freq3(1000:end));
avg_freq4=mean(freq4(8000:end));
avg_freq5=mean(freq5(1000:end));
avg_freq6=mean(freq6(1000:end));
avg_freq=[avg_freq0 avg_freq1 avg_freq2 avg_freq3 avg_freq4 avg_freq5 avg_freq6];
loglog(sigma,avg_freq,'--o','MarkerSize',10);
xlabel('\sigma')
ylabel('\Omega(\sigma)')
ax = gca;
ax.FontSize=18;
%saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_xi_L=2^13_2.fig');





