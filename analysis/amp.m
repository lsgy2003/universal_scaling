clear all;

T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = 7000;
Dt=1:dt*Evo:(T-tw)+1;

Lx=[2^8 2^9 2^10 2^11 2^12]; 
lgLx =log(Lx);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^8_amp.mat');
amp_corr2=m.avg_fluc.amp_AA(:);
Dcorr2 = -log(amp_corr2);
%
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^9_amp.mat');
amp_corr3=m.avg_fluc.amp_AA(:);
Dcorr3 = -log(amp_corr3);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^10_amp.mat');
amp_corr4=m.avg_fluc.amp_AA(:);
Dcorr4 = -log(amp_corr4);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^11_amp.mat');
amp_corr5=m.avg_fluc.amp_AA(:);
Dcorr5 = -log(amp_corr5);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^12_amp.mat');
amp_corr6=m.avg_fluc.amp_AA(:);
Dcorr6 = -log(amp_corr6);

clf;
%
figure(1)

loglog(Dt, amp_corr2,'Linewidth',1,'Color',[100, 149, 237]/255);
hold on;
loglog(Dt, amp_corr3,'Linewidth',1,'Color',[255, 191, 0]/255);
hold on;
loglog(Dt, amp_corr4,'Linewidth',1,'Color',[255, 127, 80]/255);
hold on;
loglog(Dt, amp_corr5,'Linewidth',1,'Color',[204, 204, 255]/255);
hold on;
loglog(Dt, amp_corr6,'Color',[159, 226, 191]/255);
hold on;
%}
xlabel('t')
txt = {'$\langle \frac{|P_A|(t_0,L)|P_A|(t_0+t,L)}{|P_A|^2(t_0,L)}\rangle \sim 1 $'};%{\rm const.}
text(10, 0.999999,txt,'Interpreter','latex','FontSize',22);
hold off;
ax = gca;
ax.FontSize=18;
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','southwest'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=12;
axis([1 3000 0.999997 1]);
%{
ImageID='/Users/Phantom/Desktop/Code/data/Figure/Updated/amp_CEP.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/Figure/Updated/amp_CEP.eps','ContentType','vector','BackgroundColor','none')
%}
%{
phase_corr6 = m.avg_corr_phase_AA(:);
psi_corr6 = m.avg_corr_psi_AA(:);


figure(2)
loglog(Dt,-log(psi_corr6),Dt,-log(phase_corr6),'LineWidth',1.5);
xlabel('t')
ax = gca;
ax.FontSize=18;
lgd=legend('$\chi_{AA}$','$-\log \langle |\overline{e^{i\Delta_{\theta_A}}}| \rangle$','Location','southeast','Interpreter','latex');
lgd.FontSize=18;

ImageID='/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/pattern_phase_only.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/pattern_phase_only.eps','ContentType','Vector','BackgroundColor','None');
%}

