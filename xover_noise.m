%aging effect
clear all;

jp=0.011;

T=5000;
T1=10000;
T2=2000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 800 1000 4000 8000];
Dt=1:dt*Evo:(T-tw(3))+1;
Dt1=1:dt*Evo:(T1-tw(3))+1;
Dt2=1:dt*Evo:(T2-tw(2))+1;
Lx=2^12;
sigma=[0.01 0.02 0.03 0.1];

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.01/Size/jp=0.011_tw=1000.mat');
corr1=m.corr(:);
Dcorr1 = -log(corr1);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.02/Size/jp=0.011_tw=1000.mat');
corr2=m.corr(:);
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.03/Size/jp=0.011_tw=1000.mat');
corr3=m.corr(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.1/Size/jp=0.011_tw=800.mat');
corr4=m.corr(:);
Dcorr4 = -log(corr4);

clf;

% create a default colormap 
length = 3;
dark_b = [255 51 51]/255;
light_b = [255 182 193]/255;
colors = [linspace(light_b(1),dark_b(1),length)', linspace(light_b(2),dark_b(2),length)', linspace(light_b(3),dark_b(3),length)'];

figure(1)
loglog(Dt1, Dcorr1,'LineWidth',1,'Color',colors(1,:));
hold on;
loglog(Dt, Dcorr2,'LineWidth',1,'Color',colors(2,:));
hold on;
loglog(Dt, Dcorr3,'LineWidth',1,'Color',colors(3,:));
hold on;
loglog(Dt2, Dcorr4,'LineWidth',1,'Color',[204,0,25]/255);
hold on;
%}
xlabel('t','FontSize',20)
ylabel('\chi_{AA}(t)','FontSize',20)
ax = gca;
ax.FontSize=16;
lgd=legend('\sigma=0.01','\sigma=0.02','\sigma=0.03','\sigma=0.10','Location','southeast'); 
lgd.FontSize=18;
axis([10^0 4*10^3 3*10^(-5) 6*10^0])
%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xovr_noise.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_noise.eps','ContentType','Vector','BackgroundColor','None');
%}