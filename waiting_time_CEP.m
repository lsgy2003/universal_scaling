%aging effect
clear all;

sigma=0.005;
jp=0.0097;

T=10000;
T2=5000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 1000 3000 7000];
Lx=2^12;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^11_tw=100.mat');
corr1=m.corr(:);
Dcorr1 = -log(corr1);
Dt1=1:dt*Evo:(T2-tw(1))+1;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^11_tw=1000.mat');
corr2=m.corr(:);
Dcorr2 = -log(corr2);
Dt2=1:dt*Evo:(T2-tw(2))+1;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^11_tw=3000.mat');
corr3=m.corr(:);
Dcorr3 = -log(corr3);
Dt3=1:dt*Evo:(T2-tw(3))+1;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^11_tw=7000.mat');
corr4=m.corr(:);
Dcorr4 = -log(corr4);
Dt4=1:dt*Evo:(T-tw(4))+1;

clf;

% create a default colormap 
length = 3;
dark_b = [0 102 204]/255;
light_b = [204 229 255]/255;
dark=[0 0 0];
colors = [linspace(light_b(1),dark_b(1),length)', linspace(light_b(2),dark_b(2),length)', linspace(light_b(3),dark_b(3),length)'];

figure(1)
loglog(Dt1, Dcorr1,'LineWidth',1,'Color',colors(1,:));
hold on;
loglog(Dt2, Dcorr2,'LineWidth',1,'Color',colors(2,:));
hold on;
loglog(Dt3, Dcorr3,'LineWidth',1,'Color',colors(3,:));
hold on;
loglog(Dt4, Dcorr4,'LineWidth',1,'Color',dark);
hold on;

xlabel('t','FontSize',20)
ylabel('\chi_{AA}(t_0,t_0+t)','FontSize',20)
ax = gca;
ax.FontSize=16;
axis([1 4900 10^(-6) 10^(-1)])
lgd=legend('t_0=100','t_0=1000','t_0=3000','t_0=7000','Location','northwest'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=60000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;

%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/waiting_time_CEP_2^11.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/waiting_time_CEP_2^11.eps','ContentType','Vector','BackgroundColor','None');
%}