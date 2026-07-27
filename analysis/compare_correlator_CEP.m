%compare different correlation

clear all;

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

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/S1A.mat');
corr0=m.avg.C_AA(:); %\chi_aa
Dcorr0 = -log(corr0);

corr1=m.avg.C1_AA(:); %\chi'_aa
Dcorr1 = -log(corr1);

corr2=m.avg.phase1_AA(:);
Dcorr2 = -log(corr2);

corr3=m.avg.Var1_A(:);


corr4=m.avg.C2_AA(:); %\chi''_aa
Dcorr4 = -log(corr4);
corr5=m.avg.phase2_AA(:);
Dcorr5 = -log(corr5);
corr6=m.avg.Var2_A(:);
%}

clf;

figure(1)
loglog(Dt,Dcorr1,'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
hold on;
loglog(Dt,Dcorr2,':','Color',[169, 169, 169]/255,'LineWidth',2);
hold on;
loglog(Dt,corr3/2,'--','Color',[0 0.4470 0.7410],'LineWidth',1);
hold on;
xlabel('t','FontSize',18)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$\chi^{\prime}_{AA}$','$-\log \Big|\Big\langle\overline{ {\rm exp}[{i\Delta_{\theta_a}]}}\Big\rangle\Big|$','$\rm Var^{\prime}[\Delta_{\theta_A}]/2$','Location','northwest','interpreter','latex'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
axis([10^0 3*10^3 10^(-6) 10^2])

%ImageID='/Users/Phantom/Desktop/Code/data/Figure/CEP_compare1.fig';
%saveas(gcf,ImageID);
%exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/Figure/CEP_compare1.eps','ContentType','vector','BackgroundColor','none');

figure(2)
loglog(Dt,Dcorr0,'Color','k','LineWidth',1);
hold on;
loglog(Dt,Dcorr4,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on;
loglog(Dt,Dcorr5,':','Color',[150, 75, 0]/255,'LineWidth',2);
hold on;
loglog(Dt,corr6/2,'--','Color','r','LineWidth',1);
hold on;
x1=Dt(5:100).^2*10^(-4)*0.1;
loglog(Dt(5:100),x1,'--k','LineWidth',1);
txt = {'$\beta=1$'};%,P_g(1)/2
text(0.6*10^(2),5*10^(-1),txt,'Interpreter','latex','FontSize',18);
xlabel('t','FontSize',18)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$\chi_{AA}$','$\chi^{\prime\prime}_{AA}$','$-\log \overline{|\langle {\rm exp}[{i\Delta_{\theta_a}]}\rangle|}$','$\rm Var^{\prime\prime}[\Delta_{\theta_A}]/2$','Location','northwest','interpreter','latex'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
axis([10^0 3*10^3 10^(-6) 10^2])

ImageID='/Users/Phantom/Desktop/Code/data/CEP_compare2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/CEP_compare2.eps','ContentType','vector','BackgroundColor','none');

figure(3)

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fit_error.mat');
width1 = m.width_s;
std1=m.std_s;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/CEP_C1_AA.mat');
width2 = m.width_s;
std2 =m.std_s;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/CEP_C2_AA.mat');
width3 = m.width_s;
std3 =m.std_s;


Lx=[2^8 2^9 2^10 2^11 2^12];
errorbar(Lx, width1, std1, 'square','MarkerSize',5,'Color','k');
hold on;
errorbar(Lx(1:end-1), width2, std2, 'diamond','MarkerSize',5,'Color',[0.9290 0.6940 0.1250]);
hold on;
errorbar(Lx(1:end-1), width3, std3, '^','MarkerSize',5,'Color',[0.4660 0.6740 0.1880]);
hold on;

set(gca, 'XScale','log', 'YScale','log')


P_g=polyfit(log(Lx(:)),log(width1(:)),1);
x = 2^7:2^13;
y = x.^P_g(1)*exp(1)^P_g(2);
loglog(x,y,'Color','k','LineWidth',1.5);
txt = {'$\chi_{AA}:\alpha=1.35\pm0.03$'};
text(0.7*10^3,3*10^(-4),txt,'Color','k','Interpreter','latex','FontSize',20);
txt = {'$\propto L^{2\alpha}$'};
text(1.5*10^2,10^(-1),txt,'Interpreter','latex','FontSize',25);

P_g=polyfit(log(Lx(1:end-1)),log(width2),1);
x = 2^7:2^12;
y = x.^P_g(1)*exp(1)^P_g(2);
loglog(x,y,'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
txt = {'$\chi^{\prime}_{AA}:\alpha=-0.02\pm0.03$'};%,P_g(1)/2
text(0.5*10^3,0.9*10^0,txt,'Color',[0.9290 0.6940 0.1250],'Interpreter','latex','FontSize',20);
hold on;

P_g=polyfit(log(Lx(1:end-1)),log(width3),1);
x = 2^7:2^12;
y = x.^P_g(1)*exp(1)^P_g(2);
loglog(x,y,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
txt = {'$\chi^{\prime\prime}_{AA}:\alpha=-0.02\pm0.03$'};%,P_g(1)/2
text(0.5*10^3,0.2*10^0,txt,'Color',[0.4660 0.6740 0.1880],'Interpreter','latex','FontSize',20);
hold on;
%{
P_g=polyfit(log(Lx(1:end-1)),log(width4),1);
x = 2^7:2^12;
y = x.^P_g(1)*exp(1)^P_g(2);
loglog(x,y,':','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
txt = {'$\alpha=0.01\pm0.07$'};%,P_g(1)/2
text(0.5*10^3,0.1*10^0,txt,'Color',[0.4660 0.6740 0.1880],'Interpreter','latex','FontSize',25);
hold on;
%}


xlabel('L')
%ylabel('\chi_{AA}^s(L)')
ax = gca;
ax.FontSize=20;

ImageID='/Users/Phantom/Desktop/Code/data/CEP_compare3.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/CEP_compare3.eps','ContentType','vector','BackgroundColor','none');
