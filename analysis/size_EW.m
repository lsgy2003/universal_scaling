clear all;

sigma=0.005;
jp=0.0800;

T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 1000 4000 7000];
Dt=1:dt*Evo:(T-tw(2))+1;
Lx=[2^8 2^9 2^10 2^11 2^12];

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^8.mat');
corr1_1=m.avg_fluc.C_AA(:);
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^8_2.mat');
corr1_2=m.avg_fluc.C_AA(:);
Dcorr1 = -log(corr1_1*0.4+corr1_2*0.6); %240 trials in total

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^9.mat');
corr2_1=m.avg_fluc.C_AA(:);
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^9_2.mat');
corr2_2=m.avg_fluc.C_AA(:);
Dcorr2 = -log(corr2_1*0.4+corr2_2*0.6);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^10.mat');
corr3_1=m.avg_fluc.C_AA(:);
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^10_2.mat');
corr3_2=m.avg_fluc.C_AA(:);
Dcorr3 = -log(corr3_1*0.4+corr3_2*0.6);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^11.mat');
corr4_1=m.avg_fluc.C_AA(:);
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^11_2.mat');
corr4_2=m.avg_fluc.C_AA(:);
Dcorr4 = -log(corr4_1*0.4+corr4_2*0.6);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^12.mat');
corr5_1=m.avg_fluc.C_AA(:);
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/2^12_2.mat');
corr5_2=m.avg_fluc.C_AA(:);
Dcorr5 = -log(corr5_1*0.4+corr5_2*0.6);




clf;
%%Plot the New dependence of the correlation function
figure(1)
loglog(Dt, Dcorr1);
hold on;
loglog(Dt, Dcorr2);
hold on;
loglog(Dt, Dcorr3);
hold on;
loglog(Dt, Dcorr4);
hold on;
loglog(Dt, Dcorr5);
hold on;


xlabel('t','FontSize',20)
ylabel('-log|C(tw,tw+t)|','FontSize',20)
ax = gca;
ax.FontSize=16;
lgd=legend('Lx=2^8','Lx=2^9','Lx=2^{10}','Lx=2^{11}','Lx=2^{12}','Location','northwest'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=60000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;

ImageID='/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/size_scaling.fig';
saveas(gcf,ImageID);

avg_corr1 = mean(Dcorr1(1000:end));
std_corr1 = std(Dcorr1(1000:end));

avg_corr2 = mean(Dcorr2(1000:end));
std_corr2 = std(Dcorr2(1000:end));
%}
avg_corr3 = mean(Dcorr3(1000:end));
std_corr3 = std(Dcorr3(1000:end));

avg_corr4 = mean(Dcorr4(3000:end));
std_corr4 = std(Dcorr4(3000:end));

avg_corr5 = mean(Dcorr5(5000:end));
std_corr5 = std(Dcorr5(5000:end));
%

%}
width_s = [avg_corr1 avg_corr2 avg_corr3 avg_corr4 avg_corr5];%avg_corr1 avg_corr2  avg_corr7 
std_s = [std_corr1 std_corr2 std_corr3 std_corr4 std_corr5];%std_corr1 std_corr2 std_corr7 

%% -- fitting --
x = log(Lx(:));
y = log(width_s(:));

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
save('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/fit_error.mat','width_s','std_s','alpha','standard_errors_alpha','offset');



%Plot the alpha fitting
figure(2)
x = 2^7:2^13;
y = x.^(2*alpha)*exp(1)^p(2);
loglog(x,y,'k--','LineWidth',1);
%plot(x,y,'k--','LineWidth',1);
txt = {' \alpha=', alpha, '\pm', standard_errors_alpha};%,P_g(1)/2
text(10^3,10^(-5),txt,'FontSize',14);
hold on;
errorbar(Lx(:),width_s(:),std_s(:),'o');
%set(gca, 'XScale','log', 'YScale','log')
hold off;
xlabel('L','FontSize',14)
ylabel('-log|C|_s','FontSize',14)
title('L^{2\alpha}','FontSize',16)

%
ImageID='/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/alpha_scaling.fig';
saveas(gcf,ImageID);
%}


