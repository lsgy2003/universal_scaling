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
Dt1=1:dt*Evo:(T1-tw(4))+1;
Lx=[2^8 2^8 2^10 2^11 2^12];
n = [48 96 240 480 720 960 1200];
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^8_tw=7000_1_n=48.mat');
corr2=m.avg_corr_psi_AA(:);
Dcorr2 = -log(corr2);
%}
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^8_tw=7000_1_n=96.mat');
corr3=m.avg_corr_psi_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^8_tw=7000_1_n=240.mat');
corr4=m.avg_corr_psi_AA(:);
Dcorr4 = -log(corr4);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^8_tw=7000_1_n=480.mat');
corr5=m.avg_corr_psi_AA(:);
Dcorr5 = -log(corr5);
%
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^8_tw=7000_1_n=720.mat');
corr6=m.avg_corr_psi_AA(:);
Dcorr6 = -log(corr6);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^8_tw=7000_1_n=960.mat');
corr7=m.avg_corr_psi_AA(:);
Dcorr7 = -log(corr7);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/2^8_tw=7000_1_n=1200.mat');
corr8=m.avg_corr_psi_AA(:);
Dcorr8 = -log(corr8);
%}
%
avg_corr2 = mean(corr2(500:end));
std_corr2 = std(corr2(500:end));
%}
avg_corr3 = mean(corr3(500:end));
std_corr3 = std(corr3(500:end));

avg_corr4 = mean(corr4(500:end));
std_corr4 = std(corr4(500:end));

avg_corr5 = mean(corr5(500:end));
std_corr5 = std(corr5(500:end));
%
avg_corr6 = mean(corr6(500:end));
std_corr6 = std(corr6(500:end));

avg_corr7 = mean(corr7(500:end));
std_corr7 = std(corr7(500:end));

avg_corr8 = mean(corr8(500:end));
std_corr8 = std(corr8(500:end));

%}
width_s = [avg_corr2 avg_corr3 avg_corr4 avg_corr5 avg_corr6 avg_corr7 avg_corr8];%avg_corr1 avg_corr2  avg_corr7 avg_corr6
std_width_s = [std_corr2 std_corr3 std_corr4 std_corr5 std_corr6 std_corr7 std_corr8];%std_corr1 std_corr2 std_corr7 std_corr6


%%Plot the New dependence of the correlation function
figure(1)
%loglog(Dt, corr2);
%hold on;
loglog(Dt, corr2,'.','MarkerSize',10,'Color',[100, 149, 237]/255);
hold on;
loglog(Dt, corr3,'.','MarkerSize',10,'Color',[255, 191, 0]/255);
hold on;
loglog(Dt, corr4,'.','MarkerSize',10,'Color',[255, 127, 80]/255);
hold on;
loglog(Dt, corr5,'.','MarkerSize',10,'Color',[204, 204, 255]/255);
hold on;
loglog(Dt, corr6,'.','MarkerSize',10,'Color',[159, 226, 191]/255);
hold on;
loglog(Dt, corr7,'.','MarkerSize',10,'Color','k');
hold on;
loglog(Dt, corr8,'.','MarkerSize',10,'Color','r');
hold off;
%}
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
ylabel('C_{AA}(t_0,t_0+t;L)','FontSize',20)
ax = gca;
ax.FontSize=16;
lgd=legend('N=48','N=96','N=240','N=480','N=720','N=960','N=1200','Location','southwest'); 
lgd.FontSize=14;
axis([10^1 3*10^3 0.999 1.0003]);

%{
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/ensemble_CEP_1.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/ensemble_CEP_1.eps','ContentType','vector','BackgroundColor','none');
%}

% Define the model function
model = @(params, x)  x.^(params(1))* params(2);

%%alpha fitting
% Sample data with errors
xdata = n;
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
x = 200:1300;
y=0.999965*x.^0;
%loglog(x,y,'k-','LineWidth',1);
txt = {'$\sim \rm{const.}$'};
text(500,0.9999,txt,'FontSize',25,'Interpreter','latex');
hold on;
errorbar(xdata,width_s(:),std_width_s(:),'square','MarkerSize',10,'Color','k',"MarkerFaceColor",'k');
%set(gca, 'XScale','log', 'YScale','log')
hold off;
xlabel('N')
ylabel('C_{AA}^{s}(N)')
ax = gca;
ax.FontSize=16;
axis([40 1300 0.9997 1.0001]);
box on;


%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/ensemble_CEP_2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/ensemble_CEP_2.eps','ContentType','vector','BackgroundColor','none');
%}

%save('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/New/ws_jp=0.0097_1.mat','width_s','std_width_s','alpha','standard_errors_alpha');



