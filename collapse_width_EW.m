clear all;

sigma=0.005;
jp=0.040;

T=10000;
dt = 0.001;
Evo=1000;
t=0:dt*Evo:T;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Width/2^8_AA.mat');
corr1=m.corr1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Width/2^9_AA.mat');
corr2=m.corr1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Width/2^10_AA.mat');
corr3=m.corr1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Width/2^11_AA.mat');
corr4=m.corr1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Width/2^12_AA.mat');
corr5=m.corr1(:);


clf;
figure(1)

P_g1=polyfit(log(t(10:50)),log(corr5(10:50)),1);
fprintf('%d\n',P_g1);

P_g2=polyfit(log(t(50:400)),log(corr5(50:400)),1);
fprintf('%d\n',P_g2);



loglog(t,corr1,'LineWidth',2,'Color',[100, 149, 237]/255);
hold on;
loglog(t,corr2,'LineWidth',2,'Color',[255, 191, 0]/255);
hold on;
loglog(t,corr3,'LineWidth',2,'Color',[255, 127, 80]/255);
hold on;
loglog(t,corr4,'LineWidth',2,'Color',[204, 204, 255]/255);
hold on;
loglog(t,corr5,'LineWidth',2,'Color',[159, 226, 191]/255);
hold on;
%}

x0 = (5:40);
y0 = x0.^P_g1(1)*exp(1)^P_g1(2);
plot(x0,y0,'k--','LineWidth',1);
txt = {' \beta_0=0.53'};%,P_g(1)/2
text(0.2*10^(1),4*10^(-5),txt,'FontSize',20);
hold on;

x1 = (40:400);
y1 = x1.^P_g2(1)*exp(1)^P_g2(2);
plot(x1,y1,'k--','LineWidth',1);
txt = {' \beta_1=0.23'};%,P_g(1)/2
text(0.8*10^(2),4*10^(-4),txt,'FontSize',20);
hold on;


xlabel('t')
ylabel('w_A(t)')
axis([10^0 10000 10^(-7) 2*10^(-1)]);
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northwest'); 
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;

saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Width/EW_scaling_width1.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Width/EW_scaling_width1.eps','ContentType','Vector','BackgroundColor','None');

%% create smaller axes in bottom right, and plot on it
%{
xstart=.58;
xend=.88;
ystart=.18;
yend=.48;
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on;
%}

%scaling
figure(2)
alpha=0.6;
z=2.70;
%{
x1 = t*Lx(1)^(-z);
y1 = corr1*Lx(1)^(-2*alpha);
loglog(x1,y1,'.','MarkerSize',3,'Color',[100, 149, 237]/255); 
hold on;
%}
x2 = t*Lx(2)^(-z);
y2 = corr2*Lx(2)^(-2*alpha);
loglog(x2,y2,'.','MarkerSize',3,'Color',[255, 191, 0]/255); 
hold on;
%
x3 = t*Lx(3)^(-z);
y3 = corr3*Lx(3)^(-2*alpha);
loglog(x3,y3,'.','MarkerSize',3,'Color',[255, 127, 80]/255); 
hold on;
%
x4 = t*Lx(4)^(-z);
y4 = corr4*Lx(4)^(-2*alpha);
loglog(x4,y4,'.','MarkerSize',3,'Color',[204, 204, 255]/255); 
hold on;
%{
x5 = t*Lx(5)^(-z);
y5 = corr5*Lx(5)^(-2*alpha);
loglog(x5,y5,'.','MarkerSize',3,'Color',[159, 226, 191]/255); 
hold on;
%}
txt = {'$\alpha=0.6$'};%,P_g(1)/2
text(10^(-4),10^(-8),txt,'Interpreter','latex','FontSize',25);
txt = {'$z=2.00$'};%,P_g(1)/2
text(10^(-4),0.4*10^(-8),txt,'Interpreter','latex','FontSize',25);


xlabel('t/L^{z}')
ylabel('w(t)/L^{2\alpha}')
ax = gca;
ax.FontSize=20;
%axis([10^(-3) 0.5*10^(2) 10^(-15) 10^(-10)]);
hold off;

%saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Width/EW_scaling_width3.fig');
