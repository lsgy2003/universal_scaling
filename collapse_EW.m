clear all;

sigma=0.005;
jp=0.0400;

T=10000;
%T1=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 1000 3000 4000 7000];
Dt=1:dt*Evo:(T-tw(2))+1;
%Dt1=1:dt*Evo:(T1-tw(2))+1;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);

%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^8_tw=1000.mat');
corr1=m.avg_corr_psi_AA(:);
Dcorr1 = -log(corr1);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^9_tw=1000.mat');
corr2=m.avg_corr_psi_AA(:);
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^10_tw=1000.mat');
corr3=m.avg_corr_psi_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^11_tw=1000.mat');
corr4=m.avg_corr_psi_AA(:);
Dcorr4 = -log(corr4);
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^12_tw=1000.mat');
corr5=m.avg_corr_psi_AA(:);
Dcorr5 = -log(corr5);
%}

clf;

figure(1);
P_g=polyfit(log(Dt(1000:5000)),log(Dcorr4(1000:5000)),1);
x = (200:1000);
y = x.^P_g(1)*exp(1)^P_g(2);
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
plot(x,y,'k--','LineWidth',1);
txt = {' \beta=',P_g(1)/2};%,P_g(1)/2
text(10^(2),10^(-3),txt,'FontSize',20);
xlabel('t')
ylabel('\chi_{AA}(t,L)')
%axis([10^0 3000 10^(-7) 10^(0)]);
lgd=legend('L=2^9','L=2^10','L=2^{11}','L=2^{11}','L=2^{12}','Location','northwest'); 
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;

saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_BB_scaling1.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_BB_scaling1.eps','ContentType','vector','BackgroundColor','none');

%% create smaller axes in bottom right, and plot on it
%{
xstart=.58;
xend=.88;
ystart=.18;
yend=.48;
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on;
%}
%
%scaling
f=figure(2);
f.Position = [100 100 600 300];

alpha=0.5;
z=2.0;
%
x1 = Dt*Lx(1)^(-z);
y1 = Dcorr1*Lx(1)^(-2*alpha);
h1=loglog(x1,y1,'.','MarkerSize',3,'Color',[100, 149, 237]/255); 
hold on;
%}
x2 = Dt*Lx(2)^(-z);
y2 = Dcorr2*Lx(2)^(-2*alpha);
h2=loglog(x2,y2,'.','MarkerSize',3,'Color',[255, 191, 0]/255); 
hold on;
%
x3 = Dt*Lx(3)^(-z);
y3 = Dcorr3*Lx(3)^(-2*alpha);
h3=loglog(x3,y3,'.','MarkerSize',3,'Color',[255, 127, 80]/255); 
hold on;
%
x4 = Dt*Lx(4)^(-z);
y4 = Dcorr4*Lx(4)^(-2*alpha);
h4=loglog(x4,y4,'.','MarkerSize',3,'Color',[204, 204, 255]/255); 
hold on;
%
x5 = Dt*Lx(5)^(-z);
y5 = Dcorr5*Lx(5)^(-2*alpha);
h5=loglog(x5,y5,'.','MarkerSize',3,'Color',[159, 226, 191]/255); 
hold on;


txt = {'$\alpha=1/2$'};%,P_g(1)/2
text(10^(-4),5*10^(-8),txt,'Interpreter','latex','FontSize',25);
txt = {'$z=2$'};%,P_g(1)/2
text(10^(-4),3*10^(-8),txt,'Interpreter','latex','FontSize',25);
hold off;


xlabel('t/L^{z}')
ylabel('\chi_{BB}(t,L)/L^{2\alpha}')
axis([10^(-7) 10^(-1) 5*10^(-9) 3*10^(-7)]);
ax = gca;
ax.FontSize=18;

saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_BB_scaling2.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_BB_scaling2.eps','ContentType','vector','BackgroundColor','none');

%
% Create a new figure for the legend
legendFig = figure;
axes('Position', [0 0 1 1], 'Visible', 'off'); % Invisible axes

% Plot dummy lines for the legend
hold on;
hLeg(1) = plot(nan, nan, '.', 'MarkerSize', 10, 'Color', [100, 149, 237]/255);
hLeg(2) = plot(nan, nan, '.', 'MarkerSize', 10, 'Color', [255, 191, 0]/255);
hLeg(3) = plot(nan, nan, '.', 'MarkerSize', 10, 'Color', [255, 127, 80]/255);
hLeg(4) = plot(nan, nan, '.', 'MarkerSize', 10, 'Color', [204, 204, 255]/255);
hLeg(5) = plot(nan, nan, '.', 'MarkerSize', 10, 'Color', [159, 226, 191]/255);


% Create the legend
lgd = legend(hLeg, 'L=2^8', 'L=2^9', 'L=2^{10}', 'L=2^{11}', 'L=2^{12}');
set(lgd, 'FontSize', 16, 'Location', 'NorthEast');
lgd.Position = [0.1, 0.1, 0.1, 0.3]; % Adjust size and position of legend within the new figure

% Save or display the legend figure
% Or use savefig for a .fig file:
savefig(legendFig, '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_scaling2_Legend.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_scaling2_Legend.eps','ContentType','vector','BackgroundColor','none');
%}
