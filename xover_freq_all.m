clear all;



m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/ws_all_2.mat');
sigma1=0.005^2;
j1=m.jp_all;
freq1=m.freq_s;
std_freq1=m.std_freq_s;


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.05/ws_all_2.mat');
sigma2=0.05^2;
j2=m.jp_all;
freq2=m.freq_s;
std_freq2=m.std_freq_s;


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.1/ws_all_2.mat');
sigma3=0.1^2;
j3=m.jp_all;
freq3=m.freq_s;
std_freq3=m.std_freq_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.15/ws_all_2.mat');
sigma4=0.15^2;
j4=m.jp_all;
freq4=m.freq_s;
std_freq4=m.std_freq_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.5/ws_all_2.mat');
sigma5=0.5^2;
j5=m.jp_all;
freq5=m.freq_s;
std_freq5=m.std_freq_s;


clf;


figure(1)
hold on;
h0 = plot(0.0097, freq1(7)/sigma1,'o','MarkerFaceColor',[57 181 74]/255,'MarkerSize',12,'LineWidth',1.5,'Color','k');
h1 = plot(j1, freq1/sigma1, '-diamond', 'MarkerSize',6,'MarkerFaceColor','k','LineWidth',2.5,'Color','k');
h2 = plot(j2, freq2/sigma2, '-o','MarkerSize',5, 'MarkerFaceColor',"#4C8577",'LineWidth',1.5,'Color',"#4C8577");
h3 = plot(j3, freq3/sigma3, '-o','MarkerSize',5, 'MarkerFaceColor',"#96AB30",'LineWidth',1.5,'Color',"#96AB30");
h4 = plot(j4, freq4/sigma4, '-o','MarkerSize',5, 'MarkerFaceColor',"#C5DCA0",'LineWidth',1.5,'Color',"#C5DCA0");
h5 = plot(j5, freq5/sigma5, '-^','MarkerSize',5, 'MarkerFaceColor',[137,207,240]/255,'LineWidth',1.5,'Color',[137,207,240]/255);
hold off;
uistack(h1, 'top');
uistack(h0, 'top');
box on;
xlabel('j_+')
ylabel('{\rm Var}^s[\Omega_{\rm A}](L)/\sigma')
lgd=legend([h1,h2,h3,h4,h5], {'$\sigma=2.5\!\times\!10^{-5}$','$\sigma=2.5\!\times\!10^{-3}$','$\sigma=1.0\!\times\!10^{-2}$','$\sigma=2.25\!\times\!10^{-2}$','$\sigma=2.5\!\times\!10^{-1}$'},'Location','northeast','Interpreter','latex');
axis([0 0.030 0.7*10^(-1) 10])
lgd.FontSize=18;
ax = gca;
ax.FontSize=18;
yscale('log');
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_all_2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_all_2.eps','ContentType','vector','BackgroundColor','none');
