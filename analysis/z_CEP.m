%%Data from sigma=0.005, jp=0.0097

clear all;

Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fft_peaks.mat');
peak1 = m.peaks;
std1=m.peaks_err;
z1=m.z;
std_z1=m.standard_errors_z;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fft_peaks_width.mat');
peak2 = m.peaks;
std2=m.peaks_err;
z2=m.z;
std_z2=m.standard_errors_z;


clf;

%%Fitting \alpha
figure(1);
errorbar(Lx, peak1, std1, 'o','MarkerSize',5,'LineWidth',1,'Color',[0 0 0]);
set(gca, 'XScale','log', 'YScale','log')
hold on;
errorbar(Lx,peak2, std2, 'diamond','MarkerSize',5,'LineWidth',1,'Color',[0.6350 0.0780 0.1840]);
hold on;

%%fit the full correlator
P_g=polyfit(log(Lx(:)),log(peak1(:)),1);
fprintf('%d\n',P_g);
x = 1.5*(2^7:2^12);
y = x.^P_g(1)*exp(1)^P_g(2);

%plot
loglog(x,y,'k','LineWidth',2);
hold on;
txt = {'$z=0.99\pm0.01$'};
text(0.9*10^2,10^(-2),txt,'Interpreter','latex','FontSize',18);
txt = {'$f_{\chi_{AA}} \propto L^{-z}$'};
text(0.9*10^2,0.7*10^(-2),txt,'Interpreter','latex','FontSize',18);
txt={'$f\propto k^z \propto L^{-z}$'};
text(1.3*10^3,3*10^(-2),txt,'Interpreter','latex','FontSize',20,'FontWeight','bold');
hold on;
xlabel('L')
ylabel('f_{\chi_{AA}}(L), f_{w_{A}}(L)')
ax = gca;
ax.FontSize=18;
hold on;

%fit the phase correlator
P_g=polyfit(log(Lx(:)),log(peak2(:)),1);
fprintf('%d\n',P_g);
x = 1.5*(2^7:2^12);
y = x.^P_g(1)*exp(1)^P_g(2);
plot(x,y,'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
txt = {'$z=0.98\pm0.01$'};
text(1.1*10^3,1.5*10^(-2),txt,'Interpreter','latex','FontSize',18,'Color',[0.6350 0.0780 0.1840]);
txt = {'$f_{w_{A}} \propto L^{-z}$'};
text(1.1*10^3,1.1*10^(-2),txt,'Interpreter','latex','FontSize',18,'Color',[0.6350 0.0780 0.1840]);


hold off;


ImageID='/Users/Phantom/Desktop/Code/data/z_CEP1.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/z_CEP1.eps','ContentType','vector','BackgroundColor','none');

