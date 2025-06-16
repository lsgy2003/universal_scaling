%compare different correlation

clear all;

sigma=0.005;
jp=0.00970;
T=10000;
T1=5000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw=7000;
Dt=1:dt*Evo:(T-tw(1))+1;
Dt1=1:dt*Evo:(T1-tw(1))+1;
Lx=2^11;
lgLx =log(Lx);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^11_tw=7000.mat');
corr1=m.corr_AA(:);
Dcorr1 = -log(corr1);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^11_tw=7000_2.mat');
corr2=m.avg_corr_psi_AA(:);
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^11_tw=7000_3.mat');
corr3=m.avg_corr_psi_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^11_tw=7000_4.mat');
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
xlabel('t','FontSize',18)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$-\log C_{AA}$','$-\log H_{AA}$','$-\log C^{\prime}_{AA}$','$-\log H^{\prime}_{AA}$','Location','northwest','interpreter','latex'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
axis([10^0 3*10^3 10^(-6) 10^2])

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_compare3.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_compare3.eps','ContentType','vector','BackgroundColor','none');

figure(2)

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/fit_error_tw=7000.mat');
width1 = m.width_s;
std1=m.std_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/ws_jp=0.0097_2.mat');
width2 = m.width_s;
std2 =m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/ws_jp=0.0097_3.mat');
width3 = m.width_s;
std3 =m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/ws_jp=0.0097_4.mat');
width4 = m.width_s;
std4 =m.std_width_s;


Lx=[2^8 2^9 2^10 2^11 2^12];
errorbar(Lx, width1, std1, 'square','MarkerSize',5,'Color','k');
hold on;
errorbar(Lx(1:end-1), width2, std2, 'diamond','MarkerSize',5,'Color','r');
hold on;
errorbar(Lx(1:end-1), width3, std3, '^','MarkerSize',5,'Color',[0.9290 0.6940 0.1250]);
hold on;
errorbar(Lx(1:end-1), width4, std4, 'o','MarkerSize',5,'Color',[0.4660 0.6740 0.1880]);
hold on;
set(gca, 'XScale','log', 'YScale','log')


P_g=polyfit(log(Lx(:)),log(width1(:)),1);
x = 2^7:2^13;
y = x.^P_g(1)*exp(1)^P_g(2);
loglog(x,y,'Color','k','LineWidth',1.5);
txt = {'$\alpha=1.35\pm0.05$'};
text(10^3,3*10^(-4),txt,'Color','k','Interpreter','latex','FontSize',25);
txt = {'$\propto L^{2\alpha}$'};
text(0.1*10^3,0.9*10^(0),txt,'Interpreter','latex','FontSize',30);

P_g=polyfit(log(Lx(1:end-1)),log(width2),1);
x = 2^7:2^12;
y = x.^P_g(1)*exp(1)^P_g(2);
loglog(x,y,'Color','r','LineWidth',1.5);
txt = {'$\alpha=0.01\pm0.07$'};%,P_g(1)/2
text(0.5*10^3,0.9*10^0,txt,'Color','r','Interpreter','latex','FontSize',25);
hold on;

P_g=polyfit(log(Lx(1:end-1)),log(width3),1);
x = 2^7:2^12;
y = x.^P_g(1)*exp(1)^P_g(2);
loglog(x,y,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
txt = {'$\alpha=0.01\pm0.07$'};%,P_g(1)/2
text(0.5*10^3,0.3*10^0,txt,'Color',[0.9290 0.6940 0.1250],'Interpreter','latex','FontSize',25);
hold on;

P_g=polyfit(log(Lx(1:end-1)),log(width4),1);
x = 2^7:2^12;
y = x.^P_g(1)*exp(1)^P_g(2);
loglog(x,y,':','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
txt = {'$\alpha=0.01\pm0.07$'};%,P_g(1)/2
text(0.5*10^3,0.1*10^0,txt,'Color',[0.4660 0.6740 0.1880],'Interpreter','latex','FontSize',25);
hold on;



xlabel('L')
%ylabel('\chi_{AA}^s(L)')
ax = gca;
ax.FontSize=25;

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_compare3.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_compare3.eps','ContentType','vector','BackgroundColor','none');
