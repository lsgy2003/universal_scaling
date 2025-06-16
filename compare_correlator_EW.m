%compare different correlation

clear all;

sigma=0.005;
jp=0.0400;
T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = 1000;
Dt=1:dt*Evo:(T-tw(1))+1;
Lx=2^11;
lgLx =log(Lx);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Size/2^11_tw=1000.mat');
corr1=m.corr_AA(:);
Dcorr1 = -log(corr1);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^11_tw=1000_2.mat');
Dcorr2=m.avg_corr_psi_AA(:);
%Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^11_tw=1000_3.mat');
corr3=m.avg_corr_psi_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^11_tw=1000_4.mat');
corr4=m.avg_corr_psi_AA(:);
Dcorr4 = -log(corr4);

clf;

figure(1)

loglog(Dt,Dcorr1,'Color','k','LineWidth',2);
hold on;
loglog(Dt,Dcorr2,'Color','r','LineWidth',2);
hold on;
loglog(Dt,Dcorr3,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2);
hold on;
loglog(Dt,Dcorr3,':','Color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on;
x1=Dt(20:300).^0.5*10^(-4)*0.3;
loglog(Dt(20:300),x1,'--k','LineWidth',1);
txt = {'$\beta=\frac{1}{4}$'};%,P_g(1)/2
text(4*10^(1),0.5*10^(-3),txt,'Interpreter','latex','FontSize',20);
hold on;
x2=Dt(1000:8000)*10^(-6);
loglog(Dt(1000:8000),x2,'--k','LineWidth',1);
txt = {'$\beta=\frac{1}{2}$'};%,P_g(1)/2
text(0.8*10^(3),5*10^(-3),txt,'Interpreter','latex','FontSize',20);
hold off;
xlabel('t','FontSize',20)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$-\log C_{AA}$','$-\log H_{AA}$','$-\log C^{\prime}_{AA}$','$-\log H^{\prime}_{AA}$','Location','northwest','interpreter','latex'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
%axis([10^0 3*10^3 10^(-6) 10^1])

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_compare.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_compare.eps','ContentType','vector','BackgroundColor','none');

