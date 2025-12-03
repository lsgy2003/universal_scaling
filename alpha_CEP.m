%%Data from sigma=0.005, jp=0.0097

clear all;


sigma=0.005;
jp=0.0400;

T=10000;
T1=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 1000 3000 7000];
Dt=1:dt*Evo:(T-tw(4))+1;
Dt1=1:dt*Evo:(T1-tw(4))+1;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^8_tw=1000.mat');
corr2=m.avg_corr_psi_AA(:);
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^9_tw=1000.mat');
corr3=m.avg_corr_psi_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^10_tw=1000.mat');
corr4=m.avg_corr_psi_AA(:);
Dcorr4 = -log(corr4);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^11_tw=1000.mat');
corr5=m.avg_corr_psi_AA(:);
Dcorr5 = -log(corr5);
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0400/New/2^12_tw=1000.mat');
corr6=m.avg_corr_psi_AA(:);
Dcorr6 = -log(corr6);

avg_corr2 = mean(Dcorr2(7000:end));
std_corr2 = std(Dcorr2(7000:end));

avg_corr3 = mean(Dcorr3(7000:end));
std_corr3 = std(Dcorr3(7000:end));

avg_corr4 = mean(Dcorr4(7000:end));
std_corr4 = std(Dcorr4(7000:end));

avg_corr5 = mean(Dcorr5(7000:end));
std_corr5 = std(Dcorr5(7000:end));
%
avg_corr6 = mean(Dcorr6(7000:end));
std_corr6 = std(Dcorr6(7000:end));

width_s = [avg_corr1 avg_corr2 avg_corr3 avg_corr4 avg_corr5];
std_s = [std_corr1 std_corr2 std_corr3 std_corr4 std_corr5];

%% Define the fitting model function
model = @(params, x)  x.^(2*params(1))* exp(1)^params(2);

%%alpha fitting
% Sample data with errors
xdata = Lx(:);
ydata = width_s;
y_err = std_s;  % Example errors in y

% Define the weighted error function
weighted_error_function = @(params, x, y, y_err) (model(params, x) - y) ./ y_err;
%weighted_error_function = @(params, x, y, y_err) ((model(params, x) - y) ./ y_err).^2;

% Initial guess for the parameters
initial_guess = [1, 1];

%Fit the model, calculate the standard error of the fit parameters, use the Jacobian:
[params_alpha, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(params, x) weighted_error_function(params, x, ydata, y_err), initial_guess, xdata, ydata);

% Optimal values for the parameters
disp('Optimal parameters:');
disp(params_alpha);
alpha=params_alpha(1);


% Compute the covariance matrix of the parameters
cov_matrix = inv(J.' * J);

% Extract standard errors from the covariance matrix
standard_errors = sqrt(diag(cov_matrix));

% Display the standard errors
disp('Standard errors of fit parameters:');
disp(standard_errors);

% Convert sparse matrix to full matrix
standard_errors = full(standard_errors);
standard_errors_alpha = standard_errors(1);

save('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/data/fit_error_tw=7000','width_s','std_s','alpha','standard_errors_alpha','beta','standard_errors_beta');


%Plot alpha fitting
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/data/fit_error_tw=7000.mat');
width1 = m.width_s;
std1=m.std_s;
alpha1=m.alpha;
std_alpha1=m.standard_errors_alpha;


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/data/fit_error_width_AA.mat');
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
fprintf('%d\n',P_g);
x = 2^7:2^13;
y = x.^P_g(1)*exp(1)^P_g(2);

%plot
loglog(x,y,'k','LineWidth',3);
hold on;
txt = {'$\alpha=1.35\pm0.05$'};
text(2*10^3,4*10^(-3),txt,'Interpreter','latex','FontSize',18);
txt = {'$-\log|C_{AA}| \propto L^{2\alpha}$'};
text(2*10^3,1.5*10^(-3),txt,'Interpreter','latex','FontSize',18);
hold on;
xlabel('L')
ylabel('-log|C_{AA}(L)|_s, w_{A}(L)_s')
ax = gca;
ax.FontSize=18;
hold on;

%fit the phase correlator
P_g=polyfit(log(Lx(:)),log(width2(:)),1);
fprintf('%d\n',P_g);
x = 2^7:2^13;
y = x.^P_g(1)*exp(1)^P_g(2);
plot(x,y,'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
txt = {'$\alpha=1.35\pm0.01$'};
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


ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/alpha_CEP_AA.fig';
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/alpha_CEP_AA.eps','ContentType','vector','BackgroundColor','none')

saveas(gcf,ImageID);
