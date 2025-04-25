clear all;
clf;
%%Plot phase diagram

%%sigma=0.05
jp1 = [0.030 0.040 0.050 0.070 0.100 0.150 0.200 0.250]'; 
beta1 = [0.429 0.348 0.311 0.289 0.272 0.256 0.249 0.243]'; 
sigma1= 0.05*ones(length(beta1),1);

%%sigma=0.1
jp2 = [0.040 0.050 0.060 0.070 0.100 0.150 0.200 0.250]';
beta2 = [0.728 0.534 0.418 0.371 0.270 0.254 0.246 0.242]';
sigma2= 0.1*ones(length(beta2),1);

%%sigma=0.15
jp3=[0.07 0.10 0.12 0.14 0.15 0.20 0.25]';
beta3=[0.574 0.383 0.296 0.277 0.248 0.232 0.232]'; 
sigma3= 0.15*ones(length(beta3),1);

%%sigma=0.2
jp4 = [0.100 0.120 0.150 0.200 0.250]';
beta4 = [0.520 0.410 0.342 0.271 0.239]'; 
sigma4= 0.2*ones(length(beta4),1);


%sigma=0.25
jp5 = [0.120 0.150 0.200 0.250]';
beta5 = [0.513 0.419 0.327 0.282]';
sigma5= 0.25*ones(length(beta5),1);


% Combine data sets
x_combined = [jp1; jp2; jp3; jp4; jp5];  %32
y_combined = [sigma1; sigma2; sigma3; sigma4; sigma5];
z_combined = [beta1; beta2; beta3; beta4; beta5];


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
%axis([0 0.35 0 0.25])
set(gca,'YDir','normal') 
ax = gca;
ax.FontSize=14;
xlabel('j_+','FontSize',20) %'Interpreter','latex'
ylabel('\sigma','FontSize',20)
colormap('parula');
c=colorbar;
c.Label.String='\beta';
c.Label.FontSize=20;
c.Ticks=[0.25 0.35 0.45 0.55 0.65];


ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/phase_diagram.fig';
%saveas(gcf,ImageID);