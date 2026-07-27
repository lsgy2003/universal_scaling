clear all;


m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Xover/xi=0.002_ws_all.mat');
sigma0=0.002^2;
j0=m.jp_all;
w0=m.width_s;
std_w0=m.std_width_s;
%
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Xover/xi=0.005_ws_all.mat');
sigma1=0.005^2;
j1=m.jp_all;
w1=m.width_s;
std_w1=m.std_width_s;
%}

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Xover/xi=0.010_ws_all.mat');
sigma2=0.01^2;
j2=m.jp_all;
w2=m.width_s;
std_w2=m.std_width_s;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Xover/xi=0.050_ws_all.mat');
sigma3=0.05^2;
j3=m.jp_all;
w3=m.width_s;
std_w3=m.std_width_s;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Xover/xi=0.080_ws_all.mat');
sigma4=0.08^2;
j4=m.jp_all;
w4=m.width_s;
std_w4=m.std_width_s;


m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Xover/xi=0.100_ws_all.mat');
sigma5=0.1^2;
j5=m.jp_all;
w5=m.width_s;
std_w5=m.std_width_s;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Xover/xi=0.150_ws_all.mat');
sigma6=0.15^2;
j6=m.jp_all;
w6=m.width_s;
std_w6=m.std_width_s;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Xover/xi=0.200_ws_all.mat');
sigma7=0.2^2;
j7=m.jp_all;
w7=m.width_s;
std_w7=m.std_width_s;
%
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Xover/xi=0.500_ws_all.mat');
sigma8=0.5^2;
j8=m.jp_all;
w8=m.width_s;
std_w8=m.std_width_s;
%}

clf;



figure(1)
hold on;
%low noise
h0 = plot(j0,w0/sigma0,'-+','Color',[204, 204, 255]/255,'MarkerSize',5,'LineWidth',1.5);
h1 = plot(j1,w1/sigma1,'-+','Color','k','MarkerSize',5,'LineWidth',1.5);
h2 = plot(j2,w2/sigma2,'-+','Color',[253,190,188]/255,'MarkerSize',5,'LineWidth',1.5);
%moderate noise
h3 = plot(j3,w3/sigma3,'-+','Color',"#4C8577",'MarkerSize',5,'LineWidth',1.5);
h4 = plot(j4,w4/sigma4,'-+','Color',"#53AB30",'MarkerSize',5,'LineWidth',1.5);
h5 = plot(j5,w5/sigma5,'-+','Color',"#96AB30",'MarkerSize',5,'LineWidth',1.5);
h6 = plot(j6,w6/sigma6,'-+','Color',"#C5DCA0",'MarkerSize',5,'LineWidth',1.5);
h7 = plot(j7,w7/sigma7,'-+','Color',"#E1F25C ",'MarkerSize',5,'LineWidth',1.5);

%high noise
h8 = plot(j8,w8/sigma8,'-+','Color',[137, 207, 240]/255,'MarkerSize',5,'LineWidth',1.5);
%add lables
h9 = plot(0.0097, w1(23)/sigma1,'o','MarkerFaceColor',[57 181 74]/255,'MarkerSize',10,'LineWidth',1.5,'Color','k');
%h10 = plot(0.0400, w1(end)/sigma1,'^','MarkerFaceColor',[255 191 0]/255,'MarkerSize',10,'LineWidth',1.5,'Color','k');
hold off;
uistack(h1, 'top');
uistack(h9,'top');
%uistack(h10,'top');
box on;
xlabel('j_+')
ylabel('\chi_{AA}^s/\sigma')
ax = gca;
ax.FontSize=18;
set(gca,'yscale','log')
axis([0 0.042 5 3000])
lgd=legend([h0,h1,h2,h3,h4,h5,h6,h7,h8], {'$\sigma=4.0\!\times\!10^{-6}$','$\sigma=2.5\!\times\!10^{-5}$','$\sigma=1.0\!\times\!10^{-4}$','$\sigma=2.5\!\times\!10^{-3}$','$\sigma=6.4\!\times\!10^{-3}$','$\sigma=1.0\!\times\!10^{-2}$','$\sigma=2.25\!\times\!10^{-2}$','$\sigma=4.0\!\times\!10^{-2}$','$\sigma=2.5\!\times\!10^{-1}$'},'Location','northeast','Interpreter','latex'); 
lgd.FontSize=14;
ImageID='/Users/Phantom/Desktop/Code/data/xover_all.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/xover_all.eps','ContentType','vector','BackgroundColor','none');