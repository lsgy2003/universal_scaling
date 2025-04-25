clear all;
%temporal correlation functions
Lx=2^12;
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.01;
T=10000;
%T=2000;
Evo=100;
Nt = T/dt;
%{
tw = [400, 800, 1600]; %, 800,1600,3200,6400
Ts=10;
StopT=500; 
DeltaT = StopT/2;
%}
%
tw = [1000 4000 6000];
Ts=10;
StopT=1000; 
DeltaT = StopT/2; 
%}
beta = zeros(1,length(Ts:StopT-DeltaT));
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.15/Psi/jp=0.120_tw=1000.mat');
corr1=m.corr(:);
Dcorr1 = -log(corr1);
Dt1=1:dt*Evo:(T-tw(1))+1;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.15/Psi/jp=0.120_tw=4000.mat');
corr2=m.corr(:);
Dcorr2 = -log(corr2);
Dt2=1:dt*Evo:(T-tw(2))+1;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.15/Psi/jp=0.120_tw=6000.mat');
corr3=m.corr(:);
Dcorr3 = -log(corr3);
Dt3=1:dt*Evo:(T-tw(3))+1;


%
figure(1)
loglog(Dt1,Dcorr1,'LineWidth',2);
hold on;
loglog(Dt2,Dcorr2,'LineWidth',2);
hold on;
loglog(Dt3,Dcorr3,'LineWidth',2);
hold on;
scatter([Ts StopT],[Dcorr3(Ts) Dcorr3(StopT)],80,'filled');
hold on;
x=10:100;
y=0.025*x.^(2*0.32);
loglog(x,y,'--','Color','k','LineWidth',2);
txt = {' $\propto t^{2\beta}$'};
text(10,0.3,txt,'FontSize',25,'Interpreter','latex');

hold on;
xlabel('t')
ylabel('\chi_{AA}(t_0,t_0+t)')
ax = gca;
ax.FontSize=25;
lgd=legend('t_0=1000','t_0=4000','t_0=6000','Location','southeast'); %'tw=1600','tw=3200','tw=6400',
lgd.FontSize=20;
hold off;

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/beta_fit.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/beta_fit.eps','ContentType','vector','BackgroundColor','none')

%
for i = 1:length(beta)
    P_g=polyfit(log(Dt3(Ts+i-1:Ts+i-1+DeltaT)),log(Dcorr3(Ts+i-1:Ts+i-1+DeltaT)),1);
    beta(i)=P_g(1)/2;
end

figure(2)
x=Ts:StopT-DeltaT;
plot(x, beta);
hold on;
y = 0.25*(ones(1,length(x)));
loglog(x,y,'--','LineWidth',1,'Color',[0.9290 0.6940 0.1250])
hold off;
xlabel('t');
ylabel('\beta_3');
axis([x(1) x(end) 0.2 0.8])

beta_stat=datastats(beta');
fprintf('%.3f\n',beta_stat.mean);
fprintf('%.3f\n',beta_stat.std);

filename='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.15/Psi/jp=0.120_psi_beta.mat';
save(filename,'Ts','DeltaT','x','StopT','beta','beta_stat');
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.15/Psi/jp=0.120_psi_beta.fig';
saveas(gcf,ImageID);

