clear all;

T=10000;
dt = 0.001;
Evo=1000;
t=0:dt*Evo:T;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);


m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^8_amp.mat');
corr1=m.avg_freq.width_A(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^9_amp.mat');
corr2=m.avg_freq.width_A(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^10_amp.mat');
corr3=m.avg_freq.width_A(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^11_amp.mat');
corr4=m.avg_freq.width_A(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^12_amp.mat');
corr5=m.avg_freq.width_A(:);


clf;
figure(1)

P_g=polyfit(log(t(200:1800)),log(corr5(200:1800)),1);
fprintf('%d\n',P_g);



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


x1 = (200:1800);
y1 = x1.^P_g(1)*exp(1)^P_g(2);
plot(x1,y1,'k--','LineWidth',1);
txt = {' \beta_1=0.22'};%,P_g(1)/2
text(10^(1),10^(-4),txt,'FontSize',20);
hold on;


xlabel('t')
ylabel('w_A(t,L)')
axis([10^0 10000 10^(-7) 2*10^(-1)]);
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northwest'); 
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;

saveas(gcf,'/Users/Phantom/Desktop/Code/data/EW_scaling_width.fig');
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/EW_scaling_width.eps','ContentType','Vector','BackgroundColor','None');

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
alpha=0.5;
z=2.;
%
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
%
x5 = t*Lx(5)^(-z);
y5 = corr5*Lx(5)^(-2*alpha);
loglog(x5,y5,'.','MarkerSize',3,'Color',[159, 226, 191]/255); 
hold on;
%}
txt = {'$\alpha=1/2$'};%,P_g(1)/2
text(10^(-4),10^(-8),txt,'Interpreter','latex','FontSize',25);
txt = {'$z=2$'};%,P_g(1)/2
text(10^(-4),0.4*10^(-8),txt,'Interpreter','latex','FontSize',25);


xlabel('t/L^{z}')
ylabel('w(t)/L^{2\alpha}')
ax = gca;
ax.FontSize=20;
%axis([10^(-3) 0.5*10^(2) 10^(-15) 10^(-10)]);
hold off;

%saveas(gcf,'/Users/Phantom/Desktop/Code/data/Turing/sigma=0.005/jp=0.0400/Width/EW_scaling_width3.fig');
