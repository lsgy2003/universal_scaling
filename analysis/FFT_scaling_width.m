clear all;

sigma=0.005;
jp=0.0097;

%temporal correlation functions
Lx=[2^8 2^9 2^10 2^11 2^12]; 
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.001;
T=10000;
Evo=1000;
Nt=T/dt;
t=Nt/Evo;


%%z fitting

peaks=[0.0360 0.0190 0.0096 0.0048 0.0024]; %width jp=0.0097
peaks_err=[0.001 0.001 0.0001 0.0001 0.0001]; 
peaks_amp = [1.48*10^(-6) 2.21*10^(-5) 0.0003 0.0031 0.019];
%}



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


filename='/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/CEP_fft_peaks_width.mat';
save(filename,'sigma','jp','Lx','peaks','peaks_amp','z','standard_errors_z');

clf;
figure(1)

x = 2^7:2^13;
y = x.^(-params_z(1))*exp(1)^params_z(2);
loglog(x,y,'k--','LineWidth',1);
txt = {'$z=0.975$'};%,P_g(1)
text(10^3,2*10^(-2),txt,'Interpreter','latex','FontSize',30);
txt2={'$f\sim L^{-z} \sim k^z$'};
text(10^3,4*10^(-2),txt2,'Interpreter','latex','FontSize',30);%'Interpreter','latex',
hold on;
errorbar(Lx,peaks,peaks_err,'square','MarkerSize',10,'MarkerFaceColor','r');%'MarkerFaceColor','r','MarkerSize',10
hold off;
xlabel('L')
ylabel('$f_{w_{A}}$','Interpreter','latex')
ax = gca;
ax.FontSize=30;
%axis([10^2 10^4 2*10^(-4) 5*10^(-2)]);

%{
ImageID = '/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fit_fft_width_1.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fit_fft_width_1.eps','ContentType','vector','BackgroundColor','none');
%}
figure(2)


m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^8_width_fft.mat');
P1=m.P(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^9_width_fft.mat');
P2=m.P(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^10_width_fft.mat');
P3=m.P(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^11_width_fft.mat');
P4=m.P(:);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^12_width_fft.mat');
P5=m.P(:);

f= 1/(dt*Evo)*((0:(t/2)))/t;
plot(f(2:end),P1(2:end),'LineWidth',2,'Color',[100, 149, 237]/255) ;
hold on;
plot(f(2:end),P2(2:end),'LineWidth',2,'Color',[255, 191, 0]/255) ;
hold on;
plot(f(2:end),P3(2:end),'LineWidth',2,'Color',[255, 127, 80]/255) ;
hold on;
plot(f(2:end),P4(2:end),'LineWidth',2,'Color',[204, 204, 255]/255) ;
hold on;
plot(f(2:end),P5(2:end),'LineWidth',2,'Color',[159, 226, 191]/255) ;
hold on;
plot(peaks,peaks_amp,'square','MarkerSize',8,'MarkerFaceColor','r')
hold on;

axis([0. 0.04 0 0.025]);
xlabel('$f_{w_{A}}$','Interpreter','latex')
ylabel('Amplitude')
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northeast'); 
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;
%
ImageID = '/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fft_scaling_width_2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fft_scaling_width_2.eps','ContentType','vector','BackgroundColor','none');
%}


