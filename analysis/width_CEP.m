clear all;

jp=0.0097;

T=10000;
dt = 0.01;
Evo=100;
Nt = T/dt;
t=0:dt*Evo:T;
Lx=[2^8 2^9 2^10 2^11 2^12];

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^8.mat');
corr2=m.avg_freq.width_A(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^9.mat');
corr3=m.avg_freq.width_A(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^10.mat');
corr4=m.avg_freq.width_A(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^11.mat');
corr5=m.avg_freq.width_A(:);
%
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^12_1.mat');
corr6_1=m.avg_freq.width_A(:);
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^12_2.mat');
corr6_2=m.avg_freq.width_A(:);
corr6 = corr6_1*0.6+corr6_2*0.4;


%%Plot the size dependence of the correlation function
figure(1)

loglog(t, corr2);
hold on;
loglog(t, corr3);
hold on;
loglog(t, corr4);
hold on;
loglog(t, corr5);
hold on;
loglog(t, corr6);
hold on;


xlabel('t','FontSize',20)
ylabel('w','FontSize',20)
ax = gca;
ax.FontSize=16;
lgd=legend('Lx=2^8','Lx=2^9','Lx=2^{10}','Lx=2^{11}','Lx=2^{12}','Location','northwest'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=60000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;

avg_corr2 = mean(corr2(1000:end));
std_corr2 = std(corr2(1000:end));

avg_corr3 = mean(corr3(1000:end));
std_corr3 = std(corr3(1000:end));

avg_corr4 = mean(corr4(2000:end));
std_corr4 = std(corr4(2000:end));

avg_corr5 = mean(corr5(4000:end));
std_corr5 = std(corr5(4000:end));

avg_corr6 = mean(corr6(6000:end));
std_corr6 = std(corr6(6000:end));

width_s = [avg_corr2 avg_corr3 avg_corr4 avg_corr5 avg_corr6];%avg_avg_freq.width_A avg_corr7 
std_s = [std_corr2 std_corr3 std_corr4 std_corr5 std_corr6];%std_avg_freq.width_A std_corr7 

%% fitting
x = log(Lx(:));
y = log(width_s(:));

[p,S] = polyfit(x,y,1);
[yfit,delta] = polyval(p,x,S);

slope = p(1);
alpha = slope/2;

% standard error of slope
n = length(x);
yres = y - polyval(p,x);
s2 = sum(yres.^2)/(n-2);
Sxx = sum((x-mean(x)).^2);
slope_err = sqrt(s2/Sxx);
standard_errors_alpha = slope_err/2;

fprintf('alpha = %.3f +/- %.3f\n', alpha, standard_errors_alpha);
save('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fit_error_width_AA.mat','width_s','std_s','alpha','standard_errors_alpha');

%Plot the alpha fitting
figure(2)
L = 2^7:2^13;
chi = L.^p(1)*exp(1)^p(2);
loglog(L,chi,'k','LineWidth',2);
txt = {'$\alpha_{\theta_A}=1.35\pm0.03$'};
text(10^3,1*10^(-3),txt,'Interpreter','latex','FontSize',20);
txt = {'$w_{\theta_A} \propto L^{2\alpha}$'};
text(10^3,3*10^(-3),txt,'Interpreter','latex','FontSize',20);
hold on;
errorbar(Lx(:),width_s(:),std_s(:),'o','MarkerSize',5,'LineWidth',2,'Color','k');
set(gca, 'XScale','log', 'YScale','log')
%set(gca, 'XScale','log', 'YScale','log')
hold off;
xlabel('L')
ylabel('$w_{\theta_A}(L)_s$','Interpreter','latex')
ax = gca;
ax.FontSize=18;
%title('L^{2\alpha}','FontSize',16)

%
ImageID='/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/alpha_widthAA.fig';
saveas(gcf,ImageID);