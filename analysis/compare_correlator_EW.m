%compare different correlation

clear all;


T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = 1000;
Dt=1:dt*Evo:(T-tw(1))+1;
Lx=2^11;
lgLx =log(Lx);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/S1B_1.mat');
corr0_1=m.avg.C_AA(:); %\chi_aa

corr1_1=m.avg.C1_AA(:); %\chi'_aa

corr2_1=m.avg.phase1_AA(:);

corr3_1=m.avg.Var1_A(:);


corr4_1=m.avg.C2_AA(:); %\chi''_aa

corr5_1=m.avg.phase2_AA(:);

corr6_1=m.avg.Var2_A(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/S1B_2.mat');
corr0_2=m.avg.C_AA(:); %\chi_aa
Dcorr0 = -log(corr0_1/2+corr0_2/2);

corr1_2=m.avg.C1_AA(:); %\chi'_aa
Dcorr1 = -log(corr1_1/2+corr1_2/2);

corr2_2=m.avg.phase1_AA(:);
Dcorr2 = -log(corr2_1/2+corr2_2/2);

corr3_2=m.avg.Var1_A(:);
corr3 = corr3_1/2+corr3_2/2;


corr4_2=m.avg.C2_AA(:); %\chi''_aa
Dcorr4 = -log(corr4_1/2+corr4_2/2);

corr5_2=m.avg.phase2_AA(:);
Dcorr5 = -log(corr5_1/2+corr5_2/2);

corr6_2=m.avg.Var2_A(:);
corr6 = corr6_1/2+corr6_2/2;

clf;

figure(1)
loglog(Dt,Dcorr1,'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
hold on;
loglog(Dt,Dcorr2,':','Color',[169, 169, 169]/255,'LineWidth',2);
hold on;
loglog(Dt,corr3/2,'--','Color',[0 0.4470 0.7410],'LineWidth',1);
hold on;
x1=Dt(30:300).^0.5*10^(-4)*0.25;
loglog(Dt(30:300),x1,'--k','LineWidth',1);
txt = {'$\beta=\frac{1}{4}$'};%,P_g(1)/2
text(4*10^(1),0.5*10^(-3),txt,'Interpreter','latex','FontSize',20);
hold on;
x2=Dt(2000:8000)*10^(-6)*0.6;
loglog(Dt(2000:8000),x2,'--k','LineWidth',1);
txt = {'$\beta=\frac{1}{2}$'};%,P_g(1)/2
text(1.5*10^(3),5*10^(-3),txt,'Interpreter','latex','FontSize',20);
hold off;
xlabel('t','FontSize',20)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$\chi^{\prime}_{AA}$','$-\log \Big|\Big\langle\overline{ {\rm exp}[{i\Delta_{\theta_a}]}}\Big\rangle\Big|$','$\rm Var^{\prime}[\Delta_{\theta_A}]/2$','Location','southeast','interpreter','latex'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
%axis([10^0 3*10^3 10^(-6) 10^1])

%ImageID='/Users/Phantom/Desktop/Code/data/Figure/EW_compare1.fig';
%saveas(gcf,ImageID);
%exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/Figure/EW_compare1.eps','ContentType','vector','BackgroundColor','none');

figure(2)
loglog(Dt,Dcorr0,'Color','k','LineWidth',1);
hold on;
loglog(Dt,Dcorr4,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on;
loglog(Dt,Dcorr5,':','Color',[150, 75, 0]/255,'LineWidth',2);
hold on;
loglog(Dt,corr6/2,'--','Color','r','LineWidth',1);
hold on;
x1=Dt(30:300).^0.5*10^(-4)*0.1;
loglog(Dt(30:300),x1,'--k','LineWidth',1);
txt = {'$\beta=\frac{1}{4}$'};%,P_g(1)/2
text(4*10^(1),0.2*10^(-3),txt,'Interpreter','latex','FontSize',20);
hold on;
x2=Dt(2000:8000)*10^(-6)*0.3;
loglog(Dt(2000:8000),x2,'--k','LineWidth',1);
txt = {'$\beta=\frac{1}{2}$'};%,P_g(1)/2
text(1.5*10^(3),2*10^(-3),txt,'Interpreter','latex','FontSize',20);
hold off;
xlabel('t','FontSize',20)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$\chi_{AA}$','$\chi^{\prime\prime}_{AA}$','$-\log \overline{|\langle {\rm exp}[{i\Delta_{\theta_a}]}\rangle|}$','$\rm Var^{\prime\prime}[\Delta_{\theta_A}]/2$','Location','southeast','interpreter','latex'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
axis([2 9*10^3 10^(-6) 10^(-2)])

ImageID='/Users/Phantom/Desktop/Code/data/EW_compare2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/EW_compare2.eps','ContentType','vector','BackgroundColor','none');
