clear all;


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.002/ws_all.mat');
sigma0=0.002^2;
j0=m.jp_all;
w0=m.width_s;
std_w0=m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/ws_all.mat');
sigma1=0.005^2;
j1=m.jp_all;
w1=m.width_s;
std_w1=m.std_width_s;


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.01/ws_all.mat');
sigma2=0.01^2;
j2=m.jp_all;
w2=m.width_s;
std_w2=m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.05/ws_all.mat');
sigma3=0.05^2;
j3=m.jp_all;
w3=m.width_s;
std_w3=m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.08/ws_all.mat');
sigma4=0.08^2;
j4=m.jp_all;
w4=m.width_s;
std_w4=m.std_width_s;


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.1/ws_all.mat');
sigma5=0.1^2;
j5=m.jp_all;
w5=m.width_s;
std_w5=m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.15/ws_all.mat');
sigma6=0.15^2;
j6=m.jp_all;
w6=m.width_s;
std_w6=m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.2/ws_all.mat');
sigma7=0.2^2;
j7=m.jp_all;
w7=m.width_s;
std_w7=m.std_width_s;
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.5/ws_all.mat');
sigma8=0.5^2;
j8=m.jp_all;
w8=m.width_s;
std_w8=m.std_width_s;
%}

T=2000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = 800;
Dt=1:dt*Evo:(T-tw)+1;
Lx=2^12;

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
h9 = plot(0.0097, w1(17)/sigma1,'o','MarkerFaceColor',[57 181 74]/255,'MarkerSize',10,'LineWidth',1.5,'Color','k');
h10 = plot(0.0400, w1(end)/sigma1,'^','MarkerFaceColor',[255 191 0]/255,'MarkerSize',10,'LineWidth',1.5,'Color','k');
hold off;
uistack(h1, 'top');
uistack(h9,'top');
uistack(h10,'top');
box on;
xlabel('j_+')
ylabel('\chi_{AA}^s/\sigma')
ax = gca;
ax.FontSize=18;
set(gca,'yscale','log')
axis([0 0.042 5 3000])
lgd=legend([h0,h1,h2,h3,h4,h5,h6,h7,h8], {'$\sigma=4.0\!\times\!10^{-6}$','$\sigma=2.5\!\times\!10^{-5}$','$\sigma=1.0\!\times\!10^{-4}$','$\sigma=2.5\!\times\!10^{-3}$','$\sigma=6.4\!\times\!10^{-3}$','$\sigma=1.0\!\times\!10^{-2}$','$\sigma=2.25\!\times\!10^{-2}$','$\sigma=4.0\!\times\!10^{-2}$','$\sigma=2.5\!\times\!10^{-1}$'},'Location','northeast','Interpreter','latex'); %'$\sqrt{\sigma}=0.50$',
lgd.FontSize=14;
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_all_1.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_all_1.eps','ContentType','vector','BackgroundColor','none');

%{
figure(2)
% Combine data sets
y0=sigma0*ones(length(w0),1);
y1=sigma1*ones(length(w1),1);
y2=sigma2*ones(length(w2),1);
y3=sigma3*ones(length(w3),1);
y4=sigma4*ones(length(w4),1);
y5=sigma5*ones(length(w5),1);
y6=sigma6*ones(length(w6),1);
y7=sigma7*ones(length(w7),1);

x_combined = [j0'; j1(1:end-4)'; j2'; j3'; j4'; j5';j6';j7']; 
y_combined = [y0; y1(1:end-4); y2; y3; y4; y5; y6; y7];
z_combined = [w0';w1(1:end-4)'; w2'; w3'; w4'; w5'; w6'; w7'];


% Define a common grid in the x-y plane
x_min = min(x_combined);
x_max = max(x_combined);
y_min = min(y_combined);
y_max = max(y_combined);

[X, Y] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max, 100));

% Interpolate combined z-values onto the common grid
Z_interp_combined = griddata(x_combined, y_combined, z_combined, X, Y, 'cubic');
clf;
% Display the interpolated data as a single 2D color map
imagesc(X(1,:), Y(:,1), Z_interp_combined);
axis([0.002 0.03 0.002 0.25])
set(gca,'YDir','normal') 
ax = gca;
ax.FontSize=14;
xlabel('j_+','FontSize',20) %'Interpreter','latex'
ylabel('\sigma','FontSize',20)
colormap('parula');
c=colorbar;
c.Label.String='\chi^s_{AA}';
c.Label.FontSize=18;
c.Ticks=[0 0.5 1 1.5 2];


ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/psi_phase_diagram.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/psi_phase_diagram.eps','ContentType','vector','BackgroundColor','none');
%}