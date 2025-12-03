clear all;

sigma=0.005;
jp=0.0097;

T=10000;
dt = 0.001;
Evo=1000;
t=0:dt*Evo:T;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Width/2^8.mat');
corr1=m.corr2(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Width/2^9.mat');
corr2=m.corr2(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Width/2^10.mat');
corr3=m.corr2(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Width/2^11.mat');
corr4=m.corr2(:);
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Width/2^12.mat');
corr5=m.corr2(:);


clf;
figure(1)
%P_g=polyfit(log(t(10:100)),log(corr5(10:100)),1);
%fprintf('%d\n',P_g);
%x = (10:100);
%y = x.^P_g(1)*exp(1)^P_g(2);
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
%plot(x,y,'k--','LineWidth',1);
%txt = {' \beta=0.998'};%,P_g(1)/2
%text(10^(1),5*10^(-3),txt,'FontSize',20);
xlabel('$t$','Interpreter','latex')
ylabel('$w_{\theta_{\parallel}}(t,L)$','Interpreter','latex')
%axis([10^0 10000 10^(-6) 2*10^(-1)]);
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northwest'); 
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;

saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_scaling_widthAB1.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_scaling_widthAB1.eps','ContentType','vector','BackgroundColor','none');

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
alpha=0.25;
z=2.00;

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
%{
x = (1:9)*0.005;
y = x.^(2*0.998)*10^(-10)*7;
plot(x,y,'k--','LineWidth',1);
txt = {'$\beta=0.998$'};%,P_g(1)/2
text(2*10^(-2),10^(-13),txt,'FontSize',20,'Interpreter','latex');
%}

txt = {'$\alpha_{\Delta \theta}=0.25$'};%,P_g(1)/2
text(10^(-3),0.7*10^(-6),txt,'Interpreter','latex','FontSize',25);
txt = {'$z=2.00$'};%,P_g(1)/2
text(10^(-3),0.5*10^(-6),txt,'Interpreter','latex','FontSize',25);


xlabel('$t/L^{z}$','Interpreter','latex')
ylabel('$w_{\theta_{\parallel}}(t,L)/L^{2\alpha_{\Delta \theta}}$','Interpreter','latex')
ax = gca;
ax.FontSize=18;
axis([10^(-7) 10^(0) 10^(-7) 2*10^(-6)]);
hold off;

saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_scaling_widthAB3.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_scaling_widthAB3.eps','ContentType','vector','BackgroundColor','none');
