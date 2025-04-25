%compare different correlation

clear all;

sigma=0.005;
jp=0.0400;
T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = 1000;
Dt=1:dt*Evo:(T-tw(1))+1;
Lx=2^11;
lgLx =log(Lx);

%% EW
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/Size/2^11_tw=1000.mat');
corr1=m.corr_AA(:);
Dcorr1 = -log(corr1);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^11_tw=1000.mat');
corr2=m.avg_corr_psi_AA(:);
Dcorr2 = -log(corr2);


clf;

figure(1)

loglog(Dt,Dcorr1,Dt,Dcorr2,'LineWidth',2);
hold on;
x1=Dt(20:300).^0.5*10^(-4)*0.3;
loglog(Dt(20:300),x1,'--k','LineWidth',1);
txt = {'$\beta=\frac{1}{4}$'};%,P_g(1)/2
text(4*10^(1),0.5*10^(-3),txt,'Interpreter','latex','FontSize',20);
hold on;
x2=Dt(1000:8000)*10^(-6);
loglog(Dt(1000:8000),x2,'--k','LineWidth',1);
txt = {'$\beta=\frac{1}{2}$'};%,P_g(1)/2
text(0.8*10^(3),5*10^(-3),txt,'Interpreter','latex','FontSize',20);
hold off;
xlabel('t','FontSize',20)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('-log C_{AA}','-log|g^{(1)}_{AA}|','Location','southeast'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
%axis([10^0 3*10^3 10^(-6) 10^1])

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_compare.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/EW_compare.eps','ContentType','vector','BackgroundColor','none');


%% CEP
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/Size/2^11_tw=7000.mat');
corr3=m.corr_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^11_tw=7000.mat');
corr4=m.avg_corr_psi_AA(:);
Dcorr4 = -log(corr4);


figure(2)

loglog(Dt,Dcorr3,Dt,Dcorr4,'LineWidth',2);
xlabel('t','FontSize',20)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('-logC_{AA}','-log|g^{(1)}_{AA}|','Location','southeast'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
axis([10^0 3*10^3 10^(-6) 10^1])

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_compare.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/CEP_compare.eps','ContentType','vector','BackgroundColor','none');

%% pattern


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/2^11_tw=1000.mat');
corr5=m.avg_corr_psi_AA(:);
Dcorr5 = -log(corr5);


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/New/2^11_tw=1000.mat');
corr6=m.avg_corr_psi_AA(:);
Dcorr6 = -log(corr6);

figure(3)

loglog(Dt,Dcorr5,Dt,Dcorr6,'LineWidth',2);
xlabel('t','FontSize',20)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('-logC_{AA}','-log|g^{(1)}_{AA}|','Location','southeast'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
%axis([10^0 3*10^3 10^(-6) 10^1])

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_compare.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_compare.eps','ContentType','vector','BackgroundColor','none');

figure(4)

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/ws_jp=0.040.mat');
width0 = m.width_s;
std0=m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/New/ws_jp=0.040.mat');
width1 = m.width_s;
std1 =m.std_width_s;

Lx=[2^7 2^8 2^9 2^10 2^11 2^12 2^13];
errorbar(log(Lx), width0, std0, 'square','MarkerSize',5);
hold on;
errorbar(log(Lx(2:end-1)), width1, std1, 'diamond','MarkerSize',5);
hold on;

P_g=polyfit(log(Lx(:)),width0(:),1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,'Color',[0 0.4470 0.7410],'LineWidth',1.5);
hold on;
txt = {'$\gamma=0.25 \pm 0.01$'};%,P_g(1)/2
text(7,1.5,txt,'Color',[0 0.4470 0.7410],'Interpreter','latex','FontSize',25);
hold on;


P_g=polyfit(log(Lx(2:end-1)),width1(:),1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
hold on;
txt={'$\propto 2\gamma \log L$'};
text(4.1,5.8,txt,'Interpreter','latex','FontSize',30);%'Interpreter','latex',
txt = {'$\gamma=0.26 \pm 0.01$'};%,P_g(1)/2
text(7,3.7,txt,'Color',[0.8500 0.3250 0.0980],'Interpreter','latex','FontSize',25);
hold off;

xlabel('log(L)')
ax = gca;
ax.FontSize=25;

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_compare2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern_compare2.eps','ContentType','vector','BackgroundColor','none');
