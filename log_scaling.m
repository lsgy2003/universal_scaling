%%Data from sigma=0.1, jp=0.040; sigma=0, jp=0.002

clear all;

sigma=0;
jp=0.011;
T=10000;
%T=50000;
T1=100000;
dt = 0.01;
Evo=100;
Nt = T/dt;
tw = 1000;
Dt=1:dt*Evo:(T-tw(1))+1;
Dt1=1:dt*Evo:(T1-tw(1))+1;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);

%{
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/2^7/tw=1000.mat');
corr1=m.corr(:);
Dcorr1 = -log(corr1);
%}
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/Size/2^8_tw=1000.mat');
corr2=m.avg_corr_psi_AA(:);
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/Size/2^9_tw=1000.mat');
corr3=m.avg_corr_psi_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/Size/2^10_tw=1000.mat');
corr4=m.avg_corr_psi_AA(:);
Dcorr4 = -log(corr4);
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/Size/2^11_tw=1000.mat');
corr5=m.avg_corr_psi_AA(:);
Dcorr5 = -log(corr5);
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/Size/2^12_tw=1000.mat');
corr6=m.avg_corr_psi_AA(:);
Dcorr6 = -log(corr6);
%{
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/2^13/tw=1000.mat');
corr7=m.corr(:);
Dcorr7 = -log(corr7);
%}

clf;

figure(1)
%{
loglog(Dt, Dcorr1);
hold on;
%}
loglog(Dt, Dcorr2);
hold on;
loglog(Dt, Dcorr3);
hold on;
loglog(Dt, Dcorr4);
hold on;
loglog(Dt, Dcorr5);
hold on;
loglog(Dt, Dcorr6);
hold on;
%{
loglog(Dt, Dcorr7);
hold on;
%}
xlabel('t','FontSize',20)
ylabel('\chi(t_0,t_0+t)','FontSize',20)
ax = gca;
ax.FontSize=16;
lgd=legend('Lx=2^8','Lx=2^9','Lx=2^10','Lx=2^11','Lx=2^{12}','Location','northwest'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=20000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/Size/size_scaling_tw=1000.fig';
saveas(gcf,ImageID);

%%Calculate width_s and T_s
%{
avg_corr1 = mean(Dcorr1(200:end));
std_corr1 = std(Dcorr1(200:end));
%}
avg_corr2 = mean(Dcorr2(2000:end));
std_corr2 = std(Dcorr2(2000:end));

avg_corr3 = mean(Dcorr3(2000:end));
std_corr3 = std(Dcorr3(2000:end));

avg_corr4 = mean(Dcorr4(2000:end));
std_corr4 = std(Dcorr4(2000:end));

avg_corr5 = mean(Dcorr5(2000:end));
std_corr5 = std(Dcorr5(2000:end));

avg_corr6 = mean(Dcorr6(2000:end));
std_corr6 = std(Dcorr6(2000:end));
%{
avg_corr7 = mean(Dcorr7(200:end));
std_corr7 = std(Dcorr7(200:end));
%}
width_s = [avg_corr2 avg_corr3 avg_corr4 avg_corr5 avg_corr6];%avg_corr6 avg_corr7 
std_s = [std_corr2 std_corr3 std_corr4 std_corr5 std_corr6];%std_corr6 std_corr7 

% Define the model function
model = @(params, x)  2*params(1)*log(x)+params(2);

% Sample data with errors
xdata = Lx;
ydata = width_s;
y_err = std_s;  % Example errors in y

% Define the weighted error function
%weighted_error_function = @(params, x, y, y_err) ((model(params, x) - y) ./ y_err).^2;
%weighted_error_function = @(params, x, y, y_err) (model(params, x) - y) ./ y_err;


% Initial guess for the parameters
initial_guess = [0.23, 0];

%Fit the model, calculate the standard error of the fit parameters, use the Jacobian:
%[params_gamma, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(params, x) weighted_error_function(params, x, ydata, y_err), initial_guess, xdata, ydata);
[params_gamma, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(model, initial_guess, xdata, ydata);



% Optimal values for the parameters
disp('Optimal parameters:');
disp(params_gamma);
gamma=params_gamma(1);

m = length(ydata);  % Number of data points
n = length(params_gamma);  % Number of parameters

sigma2 = sum(residual.^2) / (m - n);  % Residual variance
%C = sigma2 * inv(J' * J);  % Covariance matrix



% Compute the covariance matrix of the parameters
cov_matrix = sigma2*inv(J.' * J);

% Extract standard errors from the covariance matrix
standard_errors = sqrt(diag(cov_matrix));

% Display the standard errors
disp('Standard errors of fit parameters:');
disp(standard_errors);

% Convert sparse matrix to full matrix
standard_errors = full(standard_errors);
standard_errors_gamma = standard_errors(1);

%Plot the alpha fitting
figure(2)
x1 = linspace(4,10,100);
y1 = 2*params_gamma(1)*x1+params_gamma(2);
%y1 = 2*0.25*x1-1.26;
plot(x1,y1,'k--','LineWidth',1);
%txt = {' \alpha=', alpha, '\pm', standard_errors_alpha};%,P_g(1)/2
%text(10^3,10^(-4),txt,'FontSize',14);
hold on;
errorbar(log(Lx(:)),width_s(:),std_s(:),'o');
%set(gca, 'XScale','log', 'YScale','log')
hold off;
xlabel('L','FontSize',14)
ylabel('\chi_{AA}^s','FontSize',14)
title('L^{2\alpha}','FontSize',16)


%%Log-fit
%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/Size/ws_jp=0.011.fig';
saveas(gcf,ImageID);

%save('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/Size/ws_jp=0.011.mat','width_s','std_s','alpha','standard_errors_alpha','gamma','standard_errors_gamma');



%%Compare three parameter sets
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0.1/jp=0.040/Size/ws_jp=0.040.mat');
width1 = m.width_s;
std1=m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/sigma=0/jp=0.011/Size/ws_jp=0.011.mat');
width0 = m.width_s;
std0=m.std_width_s;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.5/jp=0.011/Size/ws_jp=0.011.mat');
width2 = m.width_s;
std2=m.std_width_s;


clf;

%%Fitting \alpha
figure(3);
errorbar(lgLx, width1, std1, 'square','MarkerSize',5,'LineWidth',1,'Color',[0 0 0]);
hold on;
errorbar(lgLx(2:end-1), width0, std0, '^','MarkerSize',5,'LineWidth',1,'Color',[0.6350 0.0780 0.1840]);
hold on;
errorbar(lgLx(2:end), width2, std2, 'diamond','MarkerSize',5,'LineWidth',1,'Color',[85, 107, 47]/255);
hold on;

%
%%fit to log scaling 
P_g=polyfit(lgLx(:),width1(:),1);
fprintf('%d\n',P_g);
x = linspace(4,10,100);
y = P_g(1)*x+P_g(2);
%chi-square goodness test
dw = (P_g(1)*lgLx+P_g(2))-width1;
chi_log= sum((dw./std1).^2)/(length(dw)-2);
%R^2 goodness test
dw2= (P_g(1)*lgLx+P_g(2))-mean(width1);
R_log = 1-sum(dw.^2)/sum(dw2.^2);
%plot
plot(x,y,'k','LineWidth',2);
hold on;
txt={'$\propto 2\gamma \log L$'};
text(5,2.4,txt,'Interpreter','latex','FontSize',22);%'Interpreter','latex',
hold on;
txt = {'$\gamma=0.25 \pm 0.01$'};%,P_g(1)/2
text(5,2.1,txt,'Interpreter','latex','FontSize',18);
xlabel('log(L)')
ylabel('\chi_{AA}^s(L)')
ax = gca;
ax.FontSize=18;
%}

%
%fit to power law
P_g=polyfit(lgLx(:),log(width1(:)),1);
fprintf('%d\n',P_g);
y = exp(P_g(1)*x+P_g(2));
%chi-square goodness test
dw = (exp(P_g(1)*lgLx(:)+P_g(2)))'-width1;
dw2 = (exp(P_g(1)*lgLx(:)+P_g(2)))'-mean(width1);
chi_power= sum((dw./std1).^2)/(length(dw)-2);
R_power = 1-sum(dw.^2)/sum(dw2.^2);

plot(x,y,'--','LineWidth',1,'Color',[0 0.4470 0.7410]);
hold on;
txt={'$\propto L^{2\alpha}$'};
text(8.5,4,txt,'Interpreter','latex','FontSize',20,'Color',[0 0.4470 0.7410]);%'Interpreter','latex',
txt={'$\alpha=0.12$'};
text(8.5,3.8,txt,'Interpreter','latex','FontSize',18,'Color',[0 0.4470 0.7410]);%'Interpreter','latex',
%}
hold on;
%
%zero noise fit
P_g=polyfit(lgLx(2:end-1),width0(:),1);
fprintf('%d\n',P_g);
x = linspace(5,9,100);
y = P_g(1)*x+P_g(2);
plot(x,y,'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
txt = {'$\gamma=0.24 \pm 0.01$'};%,P_g(1)/2
text(4.8,1,txt,'Interpreter','latex','FontSize',18,'Color',[0.6350 0.0780 0.1840]);
hold on;
%}
%%fit to log scaling
P_g=polyfit(lgLx(2:end),width2(:),1);
fprintf('%d\n',P_g);
x = linspace(5,10,100);
y = P_g(1)*x+P_g(2);
%chi-square goodness test
dw = (P_g(1)*lgLx(2:end)+P_g(2))-width2;
chi_log= sum((dw./std2).^2)/(length(dw)-2);
%R^2 goodness test
dw2= (P_g(1)*lgLx(2:end)+P_g(2))-mean(width2);
R_log = 1-sum(dw.^2)/sum(dw2.^2);
%plot
plot(x,y,'LineWidth',2,'Color',[85, 107, 47]/256);
hold off;
txt = {'$\gamma=0.26 \pm 0.01$'};%,P_g(1)/2
text(8,1.5,txt,'Interpreter','latex','FontSize',18,'Color',[85, 107, 47]/255);
xlabel('log(L)','FontSize',18)
ylabel('\chi_{AA}^s(L)','FontSize',18)
axis([4 10 0 4.5]);
hold on;


hold off;
%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/log_scaling1.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/log_scaling1.eps','ContentType','Vector','BackgroundColor','None');
%}


figure(4)
% Create a new figure for the legend
legendFig = figure(2);
axes('Position', [0 0 1 1], 'Visible', 'off'); % Invisible axes

% Plot dummy lines for the legend
hold on;
hLeg(1) = plot(nan, nan, 'square', 'LineWidth',2, 'Color', [0, 0, 0]);
hLeg(2) = plot(nan, nan, '^', 'LineWidth',2, 'Color', [0.6350 0.0780 0.1840]);
hLeg(3) = plot(nan, nan, 'diamond', 'LineWidth',2, 'Color', [85, 107, 47]/255);



% Create the legend
    lgd = legend(hLeg, {
    '$\sigma=0.1,\frac{D_A}{D_B}=1,j_+=0.040$', ...
    '$\sigma=0, \frac{D_A}{D_B}=1,j_+=0.011$', ...
    '$\sigma=0.5,\frac{D_A}{D_B}=100,j_+=0.011$', ...
    }); %\langle{\rm Var}(\Delta \theta_{AB}(L)\rangle$' '$\sigma=0.1,\frac{D_A}{D_B}=1,j_+=0.040$', ...

set(lgd, 'FontSize', 16, 'Location', 'NorthEast','Interpreter','latex');
%lgd.Position = [0.1, 0.1, 0.1, 0.3]; % Adjust size and position of legend within the new figure

% Save or display the legend figure
% Or use savefig for a .fig file:
savefig(legendFig, '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/log_scaling1_Legend.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/log_scaling1_Legend.eps','ContentType','vector','BackgroundColor','none');
%}
