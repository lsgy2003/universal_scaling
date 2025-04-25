clear all;

sigma=0.1;
jp=0.040;
T=10000;
dt = 0.01;
Evo=100;
Nt = T/dt;
tw = 1000;
Dt=1:dt*Evo:(T-tw(1))+1;
%Dt1=1:dt*Evo:(T1-tw(1))+1;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);

%{
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/2^7/tw=1000.mat');
corr1=m.corr(:);
Dcorr1 = -log(corr1);
%}
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/2^8_tw=1000.mat');
amp_corr2=m.avg_corr_amp_AA(:);
Dcorr2 = -log(amp_corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/2^9_tw=1000.mat');
amp_corr3=m.avg_corr_amp_AA(:);
Dcorr3 = -log(amp_corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/2^10_tw=1000.mat');
amp_corr4=m.avg_corr_amp_AA(:);
Dcorr4 = -log(amp_corr4);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/2^11_tw=1000.mat');
amp_corr5=m.avg_corr_amp_AA(:);
Dcorr5 = -log(amp_corr5);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/2^12_tw=1000.mat');
amp_corr6=m.avg_corr_amp_AA(:);
Dcorr6 = -log(amp_corr6);
%{
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/2^13/tw=1000.mat');
corr7=m.corr(:);
Dcorr7 = -log(corr7);
%}
clf;
%
figure(1)

loglog(Dt, amp_corr2,'Color',[100, 149, 237]/255);
hold on;
loglog(Dt, amp_corr3,'Color',[255, 191, 0]/255);
hold on;
loglog(Dt, amp_corr4,'Color',[255, 127, 80]/255);
hold on;
loglog(Dt, amp_corr5,'Color',[204, 204, 255]/255);
hold on;
loglog(Dt, amp_corr6,'Color',[159, 226, 191]/255);
hold on;
%}
xlabel('t')
txt = {'$\langle \frac{|P_A|(t_0,L)|P_A|(t_0+t,L)}{|P_A|^2(t_0,L)}\rangle \sim {\rm const.}$'};%,P_g(1)/2
text(20,0.988,txt,'Interpreter','latex','FontSize',24);

ax = gca;
ax.FontSize=18;
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','southwest'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=14;
%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/amp_scaling_tw=1000.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/amp_scaling_tw=1000.eps','ContentType','vector','BackgroundColor','none')
%}

phase_corr6 = m.avg_corr_phase_AA(:);
psi_corr6 = m.avg_corr_psi_AA(:);

%
figure(2)
loglog(Dt,phase_corr6,Dt,psi_corr6);
xlabel('t')
ax = gca;
ax.FontSize=14;
lgd=legend('$\langle |\overline{e^{i\Delta_{\theta_A}}}| \rangle$','$C_{AA}(L)$','Location','southwest','Interpreter','latex');
lgd.FontSize=16;

figure(3)
loglog(Dt,-log(phase_corr6),Dt,-log(psi_corr6),'LineWidth',1.5);
xlabel('t')
ax = gca;
ax.FontSize=18;
lgd=legend('$-\log(\langle |\overline{e^{i\Delta_{\theta_A}}}| \rangle)$','$-\log C_{AA}$','Location','northwest','Interpreter','latex');
lgd.FontSize=18;

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_phase_only.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_phase_only.eps','ContentType','Vector','BackgroundColor','None');


