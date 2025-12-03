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
ts = 0:dt*Evo:T;
Lx=2^11;
lgLx =log(Lx);


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New2/2^11_tw=1000_chi2.mat');
corr1=m.avg_corr_psi2_AA(:);
Dcorr1 = -log(corr1);

corr2=m.avg_corr_phase2_AA(:);
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New2/2^11_tw=1000_var1_unwrapped.mat');
corr3=m.avg_Var2_AA(:);
corr4 = m.avg_width_A(:);




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

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_approx.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_approx.eps','ContentType','vector','BackgroundColor','none');

