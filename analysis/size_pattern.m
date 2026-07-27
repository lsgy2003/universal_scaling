clear all;

T=2000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw=800;
Dt=1:dt*Evo:(T-tw)+1;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/2^8.mat');
corr2=m.avg_fluc.C_AA; %corr(:),corr
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/2^9.mat');
corr3=m.avg_fluc.C_AA;
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/2^10.mat');
corr4=m.avg_fluc.C_AA;
Dcorr4 = -log(corr4);
%
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/2^11.mat');
corr5=m.avg_fluc.C_AA;
Dcorr5 = -log(corr5);
%
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/2^12.mat');
corr6=m.avg_fluc.C_AA;
Dcorr6 = -log(corr6);


clf;

figure(1)
loglog(Dt, Dcorr2,'LineWidth',2,'Color',[100, 149, 237]/255);
hold on;
loglog(Dt, Dcorr3,'LineWidth',2,'Color',[255, 191, 0]/255);
hold on;
loglog(Dt, Dcorr4,'LineWidth',2,'Color',[255, 127, 80]/255);
hold on;
loglog(Dt, Dcorr5,'LineWidth',2,'Color',[204, 204, 255]/255);
hold on;
loglog(Dt, Dcorr6,'LineWidth',2,'Color',[159, 226, 191]/255);
hold on;

xlabel('t')
ylabel('\chi_{AA}(t_0,t_0+t;L)')
ax = gca;
ax.FontSize=18;
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northwest'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=16;
%axis([1 1200 10^(-3) 10^1]);
saveas(gcf,'/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/1.fig');
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/1.eps','ContentType','vector','BackgroundColor','none');

avg_corr2 = mean(Dcorr2(500:end));
std_corr2 = std(Dcorr2(500:end));

avg_corr3 = mean(Dcorr3(500:end));
std_corr3 = std(Dcorr3(500:end));

avg_corr4 = mean(Dcorr4(500:end));
std_corr4 = std(Dcorr4(500:end));

avg_corr5 = mean(Dcorr5(800:end));
std_corr5 = std(Dcorr5(800:end));

avg_corr6 = mean(Dcorr6(1000:end));
std_corr6 = std(Dcorr6(1000:end));
%

width_s = [avg_corr2 avg_corr3 avg_corr4 avg_corr5 avg_corr6];%
std_s = [std_corr2 std_corr3 std_corr4 std_corr5 std_corr6];%

%
%% Log scaling
x = log(Lx(1:end));
y = width_s(1:end);
y_err = std_s(1:end);
[p,S] = polyfit(x,y,1);
[yfit,delta] = polyval(p,x,S);

gamma = p(1)/2;

% standard error of slope
n = length(x);
yres = y - polyval(p,x);
s2 = sum(yres.^2)/(n-2);
Sxx = sum((x-mean(x)).^2);
slope_err = sqrt(s2/Sxx);
standard_errors_gamma = slope_err/2;

fprintf('gamma = %.3f +/- %.3f\n', gamma, standard_errors_gamma);
save('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/fit_error.mat','width_s','std_s','gamma','standard_errors_gamma');
%}
%{
%% Power law fitting
x = log(Lx(2:end));
y = log(width_s(2:end));

[p,S] = polyfit(x,y,1);
[yfit,delta] = polyval(p,x,S);

slope = p(1);
alpha = slope/2;
offset = p(2);

% standard error of slope
n = length(x);
yres = y - polyval(p,x);
s2 = sum(yres.^2)/(n-2);
Sxx = sum((x-mean(x)).^2);
slope_err = sqrt(s2/Sxx);
standard_errors_alpha = slope_err/2;

fprintf('alpha = %.3f +/- %.3f\n', alpha, standard_errors_alpha);
save('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/S11C.mat','width_s','std_s','alpha','standard_errors_alpha','offset');
%}

%Plot the gamma fitting
figure(2)
x1 = linspace(5,10,100);
y1 = p(1)*x1+p(2);
plot(x1,y1,'k','LineWidth',2,'Color','k');
txt = {'$\propto 2\gamma \log L$'};
text(8,1.5,txt,'FontSize',25,'interpreter','latex');
txt = {'$\gamma=0.27 \pm 0.01$'};%,P_g(1)/2
text(8,1.2,txt,'FontSize',25,'interpreter','latex');
hold on;
errorbar(x, y, y_err,'Color','k','LineWidth',2);
%hold on;
%errorbar(log(Lx(1)),width_s(1),std_width_s(1),'square','Color','r','MarkerSize',8,'MarkerFaceColor','r');
%set(gca, 'XScale','log', 'YScale','log')
%hold off;
xlabel('log L')
ylabel('\chi^s_{AA}(L)')
ax = gca;
ax.FontSize=20;

%%Log-fit
%
ImageID='/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/ws_chi1.fig';
saveas(gcf,ImageID);
%exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/Figure/ws_jp=0.0140.eps','ContentType','Vector','BackgroundColor','None');
