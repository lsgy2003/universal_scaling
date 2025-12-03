clear all;

sigma=0.005;
jp=0.0097;

T=5000;
T1=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 1000 3000 7000];
Dt=1:dt*Evo:(T-tw(4))+1;
Dt1=1:dt*Evo:(T1-tw(4))+1;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^8_tw=7000.mat');
corr1=m.corr_AA(:);
Dcorr1 = -log(corr1);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^9_tw=7000.mat');
corr2=m.corr_AA(:);
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^10_tw=7000.mat');
corr3=m.corr_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^11_tw=7000.mat');
corr4=m.corr_AA(:);
Dcorr4 = -log(corr4);
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^12_tw=7000.mat');
corr5=m.corr_AA(:);
Dcorr5 = -log(corr5);


clf;

figure(1);
%P_g=polyfit(log(Dt1(10:100)),log(Dcorr5(10:100)),1);
%x = (10:100);
%y = x.^P_g(1)*exp(1)^P_g(2);
%plot(x,y,'k--','LineWidth',1);
loglog(Dt1,Dcorr1,'LineWidth',2,'Color',[100, 149, 237]/255);
hold on;
loglog(Dt1,Dcorr2,'LineWidth',2,'Color',[255, 191, 0]/255);
hold on;
loglog(Dt1,Dcorr3,'LineWidth',2,'Color',[255, 127, 80]/255);
hold on;
loglog(Dt1,Dcorr4,'LineWidth',2,'Color',[204, 204, 255]/255);
hold on;
loglog(Dt1,Dcorr5,'LineWidth',2,'Color',[159, 226, 191]/255);
hold on;
%}
%txt = {' \beta=',P_g(1)/2};%,P_g(1)/2
%text(0.8*10^(1),10^(-2),txt,'FontSize',20);
xlabel('t')
ylabel('\chi_{AA}(t,L)')
axis([10^0 3000 10^(-7) 10^(0)]);
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northwest'); 
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;
%
saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_AA_scaling1.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_AA_scaling1.eps','ContentType','vector','BackgroundColor','none');
%}
%% create smaller axes in bottom right, and plot on it

%re-scale with L and t
f=figure(2);
f.OuterPosition = [100 100 600 400];

alpha=1.35;
z=2.00;
%
x1 = Dt1*Lx(1)^(-z);
y1 = Dcorr1*Lx(1)^(-2*alpha);
loglog(x1,y1,'.','MarkerSize',3,'Color',[100, 149, 237]/255); 
hold on;
%}
x2 = Dt1*Lx(2)^(-z);
y2 = Dcorr2*Lx(2)^(-2*alpha);
loglog(x2,y2,'.','MarkerSize',3,'Color',[255, 191, 0]/255); 
hold on;
%
x3 = Dt1*Lx(3)^(-z);
y3 = Dcorr3*Lx(3)^(-2*alpha);
loglog(x3,y3,'.','MarkerSize',3,'Color',[255, 127, 80]/255); 
hold on;
%
x4 = Dt1*Lx(4)^(-z);
y4 = Dcorr4*Lx(4)^(-2*alpha);
loglog(x4,y4,'.','MarkerSize',3,'Color',[204, 204, 255]/255); 
hold on;
%
x5 = Dt1*Lx(5)^(-z);
y5 = Dcorr5*Lx(5)^(-2*alpha);
loglog(x5,y5,'.','MarkerSize',3,'Color',[159, 226, 191]/255); 
hold on;

txt = {'$\alpha=1.35$'};%,P_g(1)/2
text(0.4*10^(-3),5*10^(-12),txt,'Interpreter','latex','FontSize',40);
txt = {'$z=2.00$'};%,P_g(1)/2
text(0.4*10^(-3),2*10^(-12),txt,'Interpreter','latex','FontSize',40);

%xlabel('t/L^{z}')
%ylabel('-log|C_{AA}(t,L)|/L^{2\alpha}')
ax = gca;
ax.FontSize=25;
axis([10^(-5) 10^(-2) 10^(-12) 10^(-10)]);
hold off;
%
saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_scaling3.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_scaling3.eps','ContentType','vector','BackgroundColor','none');

