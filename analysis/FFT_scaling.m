clear all;

sigma=0.005;
jp=0.0097;

%temporal correlation functions
Lx=[2^8 2^9 2^10 2^11 2^12];
dx=1/2;
%x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.001;
T=10000;
Evo=1000;
tw=1001;
StopT = T - tw+1;
Nt=StopT/dt;
t=Nt/Evo;


%%z fitting

%peaks=[0.0190 0.0095 0.0048 0.0024 0.00116]; %Psi, sigma=0.01,jp=0.011



peaks=[0.0188 0.0096 0.0048 0.00244 0.00122]; %Psi,sigma=0.005, jp=0.0097
peaks_amp=[3.67*10^(-8) 3.4*10^(-5) 0.0004 0.00407 0.0180];
peaks_err=[0.0005 0.0001 0.0001 0.00004 0.00002];

%fitting
xdata = log(Lx(:));
ydata = log(peaks(:));

[p,S] = polyfit(xdata,ydata,1);
[yfit,delta] = polyval(p,xdata,S);

slope = p(1);
z = slope;

% standard error of slope
n = length(xdata);
yres = ydata - polyval(p,xdata);
s2 = sum(yres.^2)/(n-2);
Sxx = sum((xdata-mean(xdata)).^2);
slope_err = sqrt(s2/Sxx);
standard_errors_z = slope_err;

fprintf('z = %.3f +/- %.3f\n', z, standard_errors_z);


filename='/Users/Phantom/Desktop/Code/data/CEP_fft_peaks.mat';
save(filename,'sigma','jp','Lx','peaks','peaks_amp','peaks_err','z','standard_errors_z');

clf;
figure(1)

x = 2^7:2^13;
y = x.^(p(1))*exp(p(2));

loglog(x,y,'k--','LineWidth',1);
txt = {'$z=0.98(1)$'};%,P_g(1)
text(10^3,0.7*10^(-2),txt,'Interpreter','latex','FontSize',30);
txt2={'$f\propto L^{-z} \propto k^z$'};
text(10^3,10^(-2),txt2,'Interpreter','latex','FontSize',30);%'Interpreter','latex',
hold on;
errorbar(Lx,peaks,peaks_err,'square','MarkerSize',10,'MarkerFaceColor','r');%'MarkerFaceColor','r','MarkerSize',10
hold off;
xlabel('L')
ylabel('f')
ax = gca;
ax.FontSize=30;
%axis([10^2 10^4 2*10^(-4) 5*10^(-2)]);

ImageID = '/Users/Phantom/Desktop/Code/data/CEP_fft_1.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/CEP_fft_1.eps','ContentType','vector','BackgroundColor','none');

figure(2)

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^8_fft.mat');
P1=m.avg_fluc.P_AA(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^9_fft.mat');
P2=m.avg_fluc.P_AA(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^10_fft.mat');
P3=m.avg_fluc.P_AA(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^11_fft.mat');
P4=m.avg_fluc.P_AA(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^12_fft.mat');
P5=m.avg_fluc.P_AA(:);

f= 1/(dt*Evo)*((0:(t/2)))/t;

plot(f(:),P1(:),'LineWidth',2,'Color',[100, 149, 237]/255) ;
hold on;
plot(f(:),P2(:),'LineWidth',2,'Color',[255, 191, 0]/255) ;
hold on;
plot(f(:),P3(:),'LineWidth',2,'Color',[255, 127, 80]/255) ;
hold on;
plot(f(:),P4(:),'LineWidth',2,'Color',[204, 204, 255]/255) ;
hold on;
plot(f(:),P5(:),'LineWidth',2,'Color',[159, 226, 191]/255) ;
hold on;
plot(peaks,peaks_amp,'square','MarkerSize',8,'MarkerFaceColor','r')
hold on;

axis([0. 0.02 0 0.02]);
xlabel('f')
ylabel('Amplitude')
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northeast');
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;

ImageID = '/Users/Phantom/Desktop/Code/data/CEP_fft_2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/CEP_fft_2.eps','ContentType','vector','BackgroundColor','none');



