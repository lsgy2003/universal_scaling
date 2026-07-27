%%Data from sigma=0.005, jp=0.0097

clear all;

Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fit_error.mat');
width1 = m.width_s;
std1=m.std_s;
alpha1=m.alpha;
std_alpha1=m.standard_errors_alpha;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fit_error_width_AA.mat');
width2 = m.width_s;
std2=m.std_s;
alpha2=m.alpha;
std_alpha2=m.standard_errors_alpha;

clf;

%%Fitting \alpha
figure(1);
errorbar(Lx, width1, std1, 'o','MarkerSize',5,'LineWidth',1,'Color',[0 0 0]);
set(gca, 'XScale','log', 'YScale','log')
hold on;
errorbar(Lx, width2, std2, 'diamond','MarkerSize',5,'LineWidth',1,'Color',[0.6350 0.0780 0.1840]);
hold on;

%%fit the full correlator
P_g=polyfit(log(Lx(:)),log(width1(:)),1);
fprintf('%d\n',P_g/2);
x = 2^7:2^13;
y = x.^P_g(1)*exp(1)^P_g(2);

%plot
loglog(x,y,'k','LineWidth',3);
hold on;
txt = {'$\alpha=1.35\pm0.02$'};
text(2*10^3,4*10^(-3),txt,'Interpreter','latex','FontSize',18);
txt = {'$\chi_{AA} \propto L^{2\alpha}$'};
text(2*10^3,1.5*10^(-3),txt,'Interpreter','latex','FontSize',18);
hold on;
xlabel('L')
ylabel('\chi_{AA}^s(L), w_{A}^s(L)')
ax = gca;
ax.FontSize=18;
hold on;

%fit the phase correlator
P_g=polyfit(log(Lx(:)),log(width2(:)),1);
fprintf('%d\n',P_g);
x = 2^7:2^13;
y = x.^P_g(1)*exp(1)^P_g(2);
plot(x,y,'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
txt = {'$\alpha=1.35\pm0.0$'};
text(0.6*10^3,10^(-1),txt,'Interpreter','latex','FontSize',18,'Color',[0.6350 0.0780 0.1840]);
txt = {'$w_{A} \propto L^{2\alpha}$'};
text(0.6*10^3,0.4*10^(-1),txt,'Interpreter','latex','FontSize',18,'Color',[0.6350 0.0780 0.1840]);

hold on;

%plot two guidance lines
x=2^7:2^10;
plot(x,x.^3*10^(-11)*0.4,'--','LineWidth',1,'Color',[0 0.4470 0.7410]);
txt={'$\alpha=\frac{3}{2}$'};
text(2*10^2,5*10^(-4),txt,'Interpreter','latex','FontSize',18,'Color',[0 0.4470 0.7410]);%'Interpreter','latex',
hold on;
plot(x,x.*10^(-8)*4,'--','LineWidth',1,'Color',[255, 127, 80]/255);
txt={'$\alpha=\frac{1}{2}$'};
text(4*10^2,5*10^(-5),txt,'Interpreter','latex','FontSize',18,'Color',[255, 127, 80]/255);%'Interpreter','latex',

hold off;


ImageID='/Users/Phantom/Desktop/Code/data/alpha_CEP_AA.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/alpha_CEP_AA.eps','ContentType','vector','BackgroundColor','none');

