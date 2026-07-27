%%compare different correlation

clear all;

sigma=0.005;
jp=0.0097;
T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = 7000;
Dt=1:dt*Evo:(T-tw(1))+1;
ts = 0:dt*Evo:T;
Lx=2^10;
lgLx =log(Lx);


m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/S1A.mat');
corr1=m.avg.C_AA(:);
Dcorr1 = -log(corr1);

corr2=m.avg.phase_AA(:);
Dcorr2 = -log(corr2);

corr3=m.avg.Var_A(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^10.mat');

corr4 = m.avg_freq.width_A(:);




clf;

figure(1)
loglog(Dt,Dcorr1,'Color','k','LineWidth',2);
hold on;
loglog(Dt,Dcorr2,':','Color',[150, 75, 0]/255,'LineWidth',2);
hold on;
loglog(Dt,corr3/2,'--','Color',[0 0.4470 0.7410],'LineWidth',1);
hold on;
loglog(ts,corr4,'Color','r','LineWidth',1);
hold off;
xlabel('t','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$\chi_{AA}$','$-\log \langle |\overline{e^{i\Delta_{\theta_A}}}| \rangle$','${\rm Var}[\Delta_{\theta_A}]/2$','$w_A$','Location','northwest','interpreter','latex'); 
lgd.FontSize=18;
%axis([10^0 3*10^3 10^(-6) 10^1])

ImageID='/Users/Phantom/Desktop/Code/data/CEP_approx.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/CEP_approx.eps','ContentType','vector','BackgroundColor','none');

