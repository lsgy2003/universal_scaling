%%Data from sigma=0.1, jp=0.200; sigma=0, jp=0.002

clear all;

sigma=0.1;
jp=0.200;
%T=2000;
T=50000;
T1=100000;
dt = 0.01;
Evo=100;
Nt = T/dt;
tw = [400 1000 10000];
Dt=1:dt*Evo:(T-tw(3))+1;
Dt1=1:dt*Evo:(T1-tw(3))+1;
Lx=[2^6 2^7 2^8 2^9 2^10]; 
lgLx =log(Lx);



m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/2^6_tw=10000.mat');
corr1=m.corr(:);
Dcorr1 = -log(corr1);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/2^7_tw=10000.mat');
corr2=m.corr(:);
Dcorr2 = -log(corr2);
%}
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/2^8_tw=10000.mat');
corr3=m.corr(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/2^9_tw=10000.mat');
corr4=m.corr(:);
Dcorr4 = -log(corr4);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/2^10_tw=10000.mat');
corr5=m.corr(:);
Dcorr5 = -log(corr5);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/fit_error_tw=10000');
width1 = m.width_s;
std1=m.std_s;


clf;

figure(1)
%
loglog(Dt, Dcorr1);
hold on;
%}
loglog(Dt, Dcorr2);
hold on;
loglog(Dt, Dcorr3);
hold on;
loglog(Dt, Dcorr4);
hold on;
loglog(Dt1, Dcorr5);
hold on;

P_g=polyfit(log(Dt(100:10000)),log(Dcorr5(100:10000)),1);
fprintf('%d\n',P_g);
x = (100:10000);
y = x.^P_g(1)*exp(1)^P_g(2);
plot(x,y,'k--','LineWidth',1);
txt = {' \beta=',P_g(1)/2};%,P_g(1)/2
text(10^(2),10^(-1),txt,'FontSize',20);
hold on;

xlabel('t','FontSize',20)
ylabel('\chi(t_0,t_0+t)|','FontSize',20)
ax = gca;
ax.FontSize=16;
lgd=legend('Lx=2^6','Lx=2^7','Lx=2^8','Lx=2^9','Lx=2^{10}','Location','northwest'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=100000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/size_scaling_tw=10000.fig';
saveas(gcf,ImageID);

%%Calculate width_s and T_s

avg_corr1 = mean(Dcorr1(5000:end));
std_corr1 = std(Dcorr1(5000:end));

avg_corr2 = mean(Dcorr2(5000:end));
std_corr2 = std(Dcorr2(5000:end));

avg_corr3 = mean(Dcorr3(10000:end));
std_corr3 = std(Dcorr3(10000:end));

avg_corr4 = mean(Dcorr4(20000:end));
std_corr4 = std(Dcorr4(20000:end));

avg_corr5 = mean(Dcorr5(40000:end));
std_corr5 = std(Dcorr5(40000:end));

width_s = [avg_corr1 avg_corr2 avg_corr3 avg_corr4 avg_corr5];%avg_corr1 avg_corr7 
std_s = [std_corr1 std_corr2 std_corr3 std_corr4 std_corr5];%std_corr1 std_corr7 

%%alpha fitting with error bar

% Define the model function
model = @(params, x)  x.^(2*params(1))* exp(1)^params(2);

% Sample data with errors
xdata = Lx;
ydata = width_s;
y_err = std_s;  % Example errors in y

% Define the weighted error function
weighted_error_function = @(params, x, y, y_err) ((model(params, x) - y) ./ y_err).^2;
%weighted_error_function = @(params, x, y, y_err) (model(params, x) - y) ./ y_err;


% Initial guess for the parameters
initial_guess = [0.5, 0];

%Fit the model, calculate the standard error of the fit parameters, use the Jacobian:
[params_alpha, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(params, x) weighted_error_function(params, x, ydata, y_err), initial_guess, xdata, ydata);
%[params_alpha, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(model, initial_guess, xdata, ydata);


% Optimal values for the parameters
disp('Optimal parameters:');
disp(params_alpha);
alpha=params_alpha(1);

m = length(ydata);  % Number of data points
n = length(params_alpha);  % Number of parameters

sigma2 = sum(residual.^2) / (m - n);  % Residual variance
cov_matrix = sigma2*inv(J.' * J);

% Compute the covariance matrix of the parameters
%cov_matrix = inv(J.' * J);'

% Extract standard errors from the covariance matrix
standard_errors = sqrt(diag(cov_matrix));

% Display the standard errors
disp('Standard errors of fit parameters:');
disp(standard_errors);

% Convert sparse matrix to full matrix
standard_errors = full(standard_errors);
standard_errors_alpha = standard_errors(1);

%Plot the alpha fitting
figure(2)
x = 2^5:2^11;
y = x.^(2*params_alpha(1))*exp(1)^params_alpha(2);
loglog(x,y,'k--','LineWidth',1);
%plot(x,y,'k--','LineWidth',1);
txt = {' \alpha=', alpha, '\pm', standard_errors_alpha};%,P_g(1)/2
text(10^3,10^(-4),txt,'FontSize',14);
hold on;
errorbar(Lx(:),width_s(:),std_s(:),'o');
%set(gca, 'XScale','log', 'YScale','log')
hold off;
xlabel('L','FontSize',14)
ylabel('-log|C|_s','FontSize',14)
title('L^{2\alpha}','FontSize',16)


%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/alpha_scaling_tw=10000.fig';
saveas(gcf,ImageID);

save('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/fit_error_tw=10000','width_s','std_s','alpha','standard_errors_alpha');
%}


%%The final EW scaling plot
figure(3)
errorbar(Lx, width1, std1, 'o','MarkerSize',5,'LineWidth',1,'Color',[0 0.4470 0.7410]);
set(gca, 'XScale','log', 'YScale','log')
hold on;

P_g=polyfit(log(Lx(:)),log(width1(:)),1);
fprintf('%d\n',P_g);
x = 2^5:2^11;
y = x.^P_g(1)*exp(1)^P_g(2);
plot(x,y,'LineWidth',2,'Color',[0 0.4470 0.7410]);
txt2={'$\propto L^{2\alpha}$'};
text(0.3*10^(3),2.5,txt2,'Interpreter','latex','FontSize',20,'Color',[0 0.4470 0.7410]);%'Interpreter','latex',
txt={'$\alpha=0.51 \pm 0.01$'};
text(0.3*10^(3),1.5,txt,'Interpreter','latex','FontSize',18,'Color',[0 0.4470 0.7410]);%'Interpreter','latex',
hold on;
xlabel('L')%'Interpreter','latex'
ylabel('-log|C_{AA}(L)|_s')
ax = gca;
ax.FontSize=16;

P_g=polyfit(lgLx(:),width1(:),1);
fprintf('%d\n',P_g);
y = P_g(1)*log(x)+P_g(2);
plot(x,y,'k--','LineWidth',1);
hold on;
txt2={'$\propto 2\gamma \log L$'};
text(10^2,10^(-2),txt2,'Interpreter','latex','FontSize',20);%'Interpreter','latex',
txt = {'$\gamma=0.13$'};%,P_g(1)/2
text(10^2,0.6*10^(-2),txt,'Interpreter','latex','FontSize',18);
hold off;
%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.200/Size/power_scaling2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/power_scaling2.eps','ContentType','Vector','BackgroundColor','None');
%}

figure(4)
% Create a new figure for the legend
legendFig = figure(2);
axes('Position', [0 0 1 1], 'Visible', 'off'); % Invisible axes

% Plot dummy lines for the legend
hold on;
hLeg(1) = plot(nan, nan, 'o', 'LineWidth',2, 'Color', [0 0.4470 0.7410]);
%hLeg(2) = plot(nan, nan, '--', 'LineWidth',1, 'Color', [0, 0, 0]);




% Create the legend
    lgd = legend(hLeg, {...
    '$\sigma=0.1,\frac{D_A}{D_B}=1,j_+ = 0.200$'});
set(lgd, 'FontSize', 18, 'Location', 'NorthEast','Interpreter','latex');
%lgd.Position = [0.1, 0.1, 0.1, 0.3]; % Adjust size and position of legend within the new figure

% Save or display the legend figure
% Or use savefig for a .fig file:
savefig(legendFig, '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/power_scaling2_Legend.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/power_scaling2_Legend.eps','ContentType','vector','BackgroundColor','none');
%}
