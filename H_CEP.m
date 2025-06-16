clear all;

sigma=0.005;
jp=0.0400;

T=10000;
T1=50000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 1000 3000 7000];
Dt=1:dt*Evo:(T-tw(4))+1;
Dt1=1:dt*Evo:(T-tw(2))+1;
Lx=[2^8 2^9 2^10 2^11];

%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^8_tw=7000_2.mat');
corr2=m.avg_corr_psi_AA(:);
Dcorr2 = -log(corr2);
%}
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^9_tw=7000_2.mat');
corr3=m.avg_corr_psi_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^10_tw=7000_2.mat');
corr4=m.avg_corr_psi_AA(:);
Dcorr4 = -log(corr4);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^11_tw=7000_2.mat');
corr5=m.avg_corr_psi_AA(:);
Dcorr5 = -log(corr5);
%
%}
%
avg_corr2 = mean(Dcorr2(1000:end));
std_corr2 = std(Dcorr2(1000:end));
%}
avg_corr3 = mean(Dcorr3(1000:end));
std_corr3 = std(Dcorr3(1000:end));

avg_corr4 = mean(Dcorr4(1000:end));
std_corr4 = std(Dcorr4(1000:end));

avg_corr5 = mean(Dcorr5(1000:end));
std_corr5 = std(Dcorr5(1000:end));


%}
width_s = [avg_corr2 avg_corr3 avg_corr4 avg_corr5];%avg_corr1 avg_corr2  avg_corr7 avg_corr6
std_width_s = [std_corr2 std_corr3 std_corr4 std_corr5];%std_corr1 std_corr2 std_corr7 std_corr6
%{
%%Fit \beta
xdata = Dt1(10:100);
ydata = Dcorr6(10:100);
xdata = xdata(:);% Ensure xdata and ydata are column vectors for consistency
ydata = ydata(:);
log_xdata = log(xdata);
log_ydata = log(ydata);
p = polyfit(log_xdata, log_ydata, 1);
slope = p(1);
intercept_beta = p(2);
beta = p(1)/2;

y_fit = exp(intercept_beta + slope * log_xdata);% Calculate the covariance matrix
residuals = log_ydata - log(y_fit);
N = length(xdata);
SSE = sum(residuals.^2); 
sigma_squared = SSE / (N - 2);
cov_matrix = sigma_squared*inv([N, sum(log_xdata); sum(log_xdata), sum(log_xdata.^2)]);

standard_errors = sqrt(diag(cov_matrix));% Standard errors of the fit parameters
standard_errors_beta = standard_errors(1)/2;

disp('fit parameters with errors:');% Display the standard errors
disp(['\beta: ', num2str(beta),',', 'error of \beta: ', num2str(standard_errors_beta)]);
%}

%%Plot the New dependence of the correlation function
figure(1)
loglog(Dt, Dcorr2,'.','MarkerSize',10,'Color',[100, 149, 237]/255);
hold on;
loglog(Dt, Dcorr3,'.','MarkerSize',10,'Color',[255, 191, 0]/255);
hold on;
loglog(Dt, Dcorr4,'.','MarkerSize',10,'Color',[255, 127, 80]/255);
hold on;
loglog(Dt, Dcorr5,'.','MarkerSize',10,'Color',[204, 204, 255]/255);
hold on;
%{
%plot \beta
x = (10:100);
y = x.^(2*beta)*exp(1)^intercept_beta;
plot(x,y,'k--','LineWidth',1);
txt = {' \beta=',beta};%
text(10^(2),10^(-2),txt,'FontNew',20);
hold off;
%}

xlabel('t','FontSize',20)
ylabel('-log H_{AA}(t_0,t_0+t;L)','FontSize',20)
ax = gca;
ax.FontSize=16;
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','Location','northwest'); 
lgd.FontSize=18;
axis([10^0 3*10^3 10^(-6) 10^2]);

%{
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/g1_CEP.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/g1_CEP.eps','ContentType','vector','BackgroundColor','none');
%}

% Define the model function
model = @(params, x)  x.^(params(1))* params(2);

%%alpha fitting
% Sample data with errors
xdata = Lx;
ydata = width_s;
y_err = std_width_s;  % Example errors in y

% Define the weighted error function
%weighted_error_function = @(params, x, y, y_err) (model(params, x) - y) ./ y_err;
%weighted_error_function = @(params, x, y, y_err) ((model(params, x) - y) ./ y_err).^2;

% Initial guess for the parameters
initial_guess = [0.5, 1];

%Fit the model, calculate the standard error of the fit parameters, use the Jacobian:
%[params_alpha, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(params, x) weighted_error_function(params, x, ydata, y_err), initial_guess, xdata, ydata);
[params_alpha, resnorm, residual, exitflag, output, lambda, J] =lsqcurvefit(model, initial_guess, xdata, ydata);

% Optimal values for the parameters
disp('Optimal parameters:');
disp(params_alpha);
alpha=params_alpha(1);


m = length(ydata);  % Number of data points
n = length(params_alpha);  % Number of parameters

sigma2 = sum(residual.^2) / (m - n);  % Residual variance
% Compute the covariance matrix of the parameters
cov_matrix = sigma2*inv(J.' * J);

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
x = Lx;
y = x.^(params_alpha(1))*params_alpha(2);
loglog(x,y,'k--','LineWidth',1);
%plot(x,y,'k--','LineWidth',1);
txt = {' \alpha=', alpha, '\pm', standard_errors_alpha};%,P_g(1)/2
text(10^2,2.5,txt,'FontSize',14);
hold on;
errorbar(xdata,width_s(:),std_width_s(:),'o');
%set(gca, 'XScale','log', 'YScale','log')
hold off;
xlabel('L','FontSize',14)
%ylabel('-log C_{AA}','FontSize',14)
title('L^{2\alpha}','FontSize',16)


%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/ws_jp=0.0097_2.fig';
saveas(gcf,ImageID);

save('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/ws_jp=0.0097_2.mat','width_s','std_width_s','alpha','standard_errors_alpha');
%}


