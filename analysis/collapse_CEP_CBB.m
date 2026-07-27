clear all;

sigma=0.005;
jp=0.0097;

T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 1000 3000 7000];
Dt=1:dt*Evo:(T-tw(4))+1;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);


m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^8_amp.mat');
corr1=m.avg_fluc.C_BB(:);
Dcorr1 = -log(corr1);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^9_amp.mat');
corr2=m.avg_fluc.C_BB(:);
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^10_amp.mat');
corr3=m.avg_fluc.C_BB(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^11_amp.mat');
corr4=m.avg_fluc.C_BB(:);
Dcorr4 = -log(corr4);
%
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^12_amp.mat');
corr5=m.avg_fluc.C_BB(:);
Dcorr5 = -log(corr5);



clf;

figure(1);
loglog(Dt,Dcorr1,'LineWidth',2,'Color',[100, 149, 237]/255);
hold on;
loglog(Dt,Dcorr2,'LineWidth',2,'Color',[255, 191, 0]/255);
hold on;
loglog(Dt,Dcorr3,'LineWidth',2,'Color',[255, 127, 80]/255);
hold on;
loglog(Dt,Dcorr4,'LineWidth',2,'Color',[204, 204, 255]/255);
hold on;
loglog(Dt,Dcorr5,'LineWidth',2,'Color',[159, 226, 191]/255);
hold on;
%}

xlabel('t')
ylabel('\chi_{BB}(t_0,t_0+t;L)')
axis([2*10^0 3000 10^(-6) 10^(0)]);
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northwest'); 
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;
%
saveas(gcf,'/Users/Phantom/Desktop/Code/data/CEP_BB_scaling1.fig');
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/CEP_BB_scaling1.eps','ContentType','vector','BackgroundColor','none');


f=figure(2);
f.OuterPosition = [100 100 600 400];

alpha=1.35;
z=2;
%
x1 = Dt*Lx(1)^(-z);
y1 = Dcorr1*Lx(1)^(-2*alpha);
loglog(x1,y1,'.','MarkerSize',5,'Color',[100, 149, 237]/255); 
hold on;
%}
x2 = Dt*Lx(2)^(-z);
y2 = Dcorr2*Lx(2)^(-2*alpha);
loglog(x2,y2,'.','MarkerSize',5,'Color',[255, 191, 0]/255); 
hold on;
%
x3 = Dt*Lx(3)^(-z);
y3 = Dcorr3*Lx(3)^(-2*alpha);
loglog(x3,y3,'.','MarkerSize',5,'Color',[255, 127, 80]/255); 
hold on;
%
x4 = Dt*Lx(4)^(-z);
y4 = Dcorr4*Lx(4)^(-2*alpha);
loglog(x4,y4,'.','MarkerSize',5,'Color',[204, 204, 255]/255); 
hold on;
%
x5 = Dt*Lx(5)^(-z);
y5 = Dcorr5*Lx(5)^(-2*alpha);
loglog(x5,y5,'.','MarkerSize',5,'Color',[159, 226, 191]/255); 
hold on;

txt = {'$\alpha=1.35(2)$'};%,P_g(1)/2
text(10^(-4),10^(-12),txt,'Interpreter','latex','FontSize',40);
txt = {'$z=2$'};%,P_g(1)/2
text(10^(-4),0.3*10^(-12),txt,'Interpreter','latex','FontSize',40);

%xlabel('t/L^{z}')
%ylabel('\chi_{BB}(t_0,t_0+t;L)/L^{2\alpha}')
ax = gca;
ax.FontSize=18;
axis([10^(-6) 10^(-2) 10^(-13) 5*10^(-11)]);
%axis([3*10^(-3) 10^(1) 10^(-13) 10^(-10)]);
hold off;
%
saveas(gcf,'/Users/Phantom/Desktop/Code/data/CEP_BB_scaling3.fig');
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/CEP_BB_scaling3.eps','ContentType','vector','BackgroundColor','none');
%}


