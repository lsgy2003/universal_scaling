%compare different correlation

clear all;

sigma=0.1;
jp=0.0400;
dt = 0.001;
Evo=1000;
T=2000;
tw = 800;
Dt=1:dt*Evo:(T-tw)+1;
Lx=2^11;
lgLx =log(Lx);


m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/2^10.mat');
corr1=m.avg_fluc.C_AA(:); %C_AA
Dcorr1 = -log(corr1);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/S4A.mat');
corr2=m.avg.C1_AA(:); %C1_AA
Dcorr2 = -log(corr2);
corr3=m.avg.phase1_AA(:);
Dcorr3 = -log(corr3);

corr4=m.avg.C2_AA(:); %C2_AA
Dcorr4 = -log(corr4);
corr5=m.avg.phase2_AA(:);
Dcorr5 = -log(corr5);
%}

clf;
figure(1)
loglog(Dt,Dcorr1,'Color','k','LineWidth',1);
hold on;
loglog(Dt,Dcorr2,'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
hold on;
loglog(Dt,Dcorr3,':','Color',[169, 169, 169]/255,'LineWidth',2);
hold on;
loglog(Dt,Dcorr4,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on;
loglog(Dt,Dcorr5,':','Color',[150, 75, 0]/255,'LineWidth',2);
hold off;
%}

xlabel('t','FontSize',20)
%ylabel('-log|C(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=18;
lgd=legend('$\chi_{AA}$','$\chi^{\prime}_{AA}$','$-\log \Big|\Big\langle\overline{ {\rm exp}[{i\Delta_{\theta_a}]}}\Big\rangle\Big|$','$\chi^{\prime\prime}_{AA}$','$-\log \overline{|\langle {\rm exp}[{i\Delta_{\theta_a}]}\rangle|}$','Location','southeast','interpreter','latex'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
axis([2 1200 10^(-2) 10^1])

%
ImageID='/Users/Phantom/Desktop/Code/data/pattern_compare.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/pattern_compare.eps','ContentType','vector','BackgroundColor','none');
%}


figure(2)

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/fit_error.mat');
width1 = m.width_s(1:5);
std1=m.std_width_s(1:5);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/S15B_C1_AA.fit_error.mat');
width2 = m.width_s;
std2 =m.std_width_s;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/S15B_C2_AA.fit_error.mat');
width3 = m.width_s;
std3 =m.std_width_s;


Lx=[2^8 2^9 2^10 2^11 2^12];
errorbar(log(Lx), width1, std1, 'square','MarkerSize',5,'Color','k');
hold on;
errorbar(log(Lx), width2, std2, 'diamond','MarkerSize',5,'Color',[0.9290 0.6940 0.1250]);
hold on;
errorbar(log(Lx), width3, std3, '^','MarkerSize',5,'Color',[0.4660 0.6740 0.1880]);
hold on;
%errorbar(log(Lx(2:end-1)), width4, std4, 'o','MarkerSize',5,'Color',[0.4660 0.6740 0.1880]);
%hold on;

P_g=polyfit(log(Lx(:)),width1(:),1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,'Color','k','LineWidth',1.5);
hold on;
txt = {'$\chi_{AA}(L):\gamma=0.26 \pm 0.01$'};%,P_g(1)/2
text(7,1,txt,'Color','k','Interpreter','latex','FontSize',20);
hold on;

P_g=polyfit(log(Lx(:)),width2,1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
txt = {'$\chi_{AA}^{\prime}(L):\gamma=0.31 \pm 0.03$'};%,P_g(1)/2
text(5,4.8,txt,'Color',[0.9290 0.6940 0.1250],'Interpreter','latex','FontSize',20);
hold on;


P_g=polyfit(log(Lx(:)),width3(:),1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
hold on;
txt={'$\propto 2\gamma \log L$'};
text(4.1,5.5,txt,'Interpreter','latex','FontSize',25);%'Interpreter','latex',
txt = {'$\chi_{AA}^{\prime\prime}(L):\gamma= -0.003 \pm 0.005$'};%,P_g(1)/2
text(4.5,1.8,txt,'Color',[0.4660 0.6740 0.1880],'Interpreter','latex','FontSize',20);
hold on;
%{
P_g=polyfit(log(Lx(2:end-1)),width4(:),1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
plot(x,y,':','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
hold on;
txt={'$\propto 2\gamma \log L$'};
text(4.1,5.8,txt,'Interpreter','latex','FontSize',30);%'Interpreter','latex',
txt = {'$\gamma=0.26 \pm 0.01$'};%,P_g(1)/2
text(7,3.6,txt,'Color',[0.4660 0.6740 0.1880],'Interpreter','latex','FontSize',25);
hold off;
%}

xlabel('log(L)')
%ylabel('\chi_{AA}^s(L)')
ax = gca;
ax.FontSize=18;

ImageID='/Users/Phantom/Desktop/Code/data/pattern_compare2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/pattern_compare2.eps','ContentType','vector','BackgroundColor','none');
