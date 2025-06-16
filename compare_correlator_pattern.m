%compare different correlation

clear all;

sigma=0.1;
jp=0.0400;
T=10000;
%T=2000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [1000];
Dt=1:dt*Evo:(T-tw(1))+1;
Lx=2^11;
lgLx =log(Lx);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/2^11_tw=1000.mat');
corr1=m.avg_corr_psi_AA(:);
Dcorr1 = -log(corr1);

%corr3=m.avg_corr_phase_AA(:);
%Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/New/2^11_tw=1000_2.mat');
Dcorr2=m.avg_corr_psi_AA(:);
%Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/New/2^11_tw=1000_3.mat');
corr3=m.avg_corr_psi_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/New/2^11_tw=1000_4.mat');
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
hold off;

xlabel('t','FontSize',20)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$-\log C_{AA}$','$-\log H_{AA}$','$-\log C^{\prime}_{AA}$','$-\log H^{\prime}_{AA}$','Location','northwest','interpreter','latex'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
%axis([10^0 3*10^3 10^(-6) 10^1])

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_compare.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_compare.eps','ContentType','vector','BackgroundColor','none');

%{
figure(2)

loglog(Dt,Dcorr1,Dt,Dcorr3,'LineWidth',2);
xlabel('t','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$-\log C_{AA}$','$-\log \langle |\overline{e^{i\Delta_{\theta_A}}}|\rangle$','Location','southeast','interpreter','latex');
lgd.FontSize=25;
%axis([10^0 3*10^3 10^(-6) 10^1])

%ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_phase_only.fig';
%saveas(gcf,ImageID);
%exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_phase_only.eps','ContentType','vector','BackgroundColor','none');
%}

figure(2)

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/ws_jp=0.040.mat');
width1 = m.width_s;
std1=m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/New/ws_jp=0.040_2.mat');
width2 = m.width_s;
std2 =m.std_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/New/ws_jp=0.040_3.mat');
width3 = m.width_s;
std3 =m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/New/ws_jp=0.040_4.mat');
width4 = m.width_s;
std4 =m.std_width_s;


Lx=[2^7 2^8 2^9 2^10 2^11 2^12 2^13];
errorbar(log(Lx), width1, std1, 'square','MarkerSize',5,'Color','k');
hold on;
errorbar(log(Lx(2:end-1)), width2, std2, 'diamond','MarkerSize',5,'Color','r');
hold on;
errorbar(log(Lx(2:end-1)), width3, std3, '^','MarkerSize',5,'Color',[0.9290 0.6940 0.1250]);
hold on;
errorbar(log(Lx(2:end-1)), width4, std4, 'o','MarkerSize',5,'Color',[0.4660 0.6740 0.1880]);
hold on;

P_g=polyfit(log(Lx(:)),width1(:),1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,'Color','k','LineWidth',1.5);
hold on;
txt = {'$\gamma=0.25 \pm 0.01$'};%,P_g(1)/2
text(7,1.5,txt,'Color','k','Interpreter','latex','FontSize',25);
hold on;

P_g=polyfit(log(Lx(2:end-1)),width2,1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,'Color','r','LineWidth',1.5);
txt = {'$\gamma=-0.01 \pm 0.01$'};%,P_g(1)/2
text(5,2.8,txt,'Color','r','Interpreter','latex','FontSize',25);
hold on;


P_g=polyfit(log(Lx(2:end-1)),width3(:),1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
hold on;
txt={'$\propto 2\gamma \log L$'};
text(4.1,5.8,txt,'Interpreter','latex','FontSize',30);%'Interpreter','latex',
txt = {'$\gamma=0.26 \pm 0.01$'};%,P_g(1)/2
text(7,4.1,txt,'Color',[0.9290 0.6940 0.1250],'Interpreter','latex','FontSize',25);
hold on;

P_g=polyfit(log(Lx(2:end-1)),width4(:),1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,':','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
hold on;
txt={'$\propto 2\gamma \log L$'};
text(4.1,5.8,txt,'Interpreter','latex','FontSize',30);%'Interpreter','latex',
txt = {'$\gamma=0.26 \pm 0.01$'};%,P_g(1)/2
text(7,3.6,txt,'Color',[0.4660 0.6740 0.1880],'Interpreter','latex','FontSize',25);
hold off;


xlabel('log(L)')
%ylabel('\chi_{AA}^s(L)')
ax = gca;
ax.FontSize=25;

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_compare2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_compare2.eps','ContentType','vector','BackgroundColor','none');
