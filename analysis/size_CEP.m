clear all;
clf;

T=10000;
T1=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [800 1000 4000 7000];
Dt=1:dt*Evo:(T-tw(4))+1;
Dt1=1:dt*Evo:(T1-tw(4))+1;
Lx=[2^8 2^9 2^10 2^11 2^12];%2^12

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^8.mat');
corr2=m.avg_fluc.C_AA(:);
Dcorr2 = -log(corr2);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^9.mat');
corr3=m.avg_fluc.C_AA(:);
Dcorr3 = -log(corr3);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^10.mat');
corr4=m.avg_fluc.C_AA(:);
Dcorr4 = -log(corr4);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^11.mat');
corr5=m.avg_fluc.C_AA(:);
Dcorr5 = -log(corr5);
%
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^12_1.mat');
corr6=m.avg_fluc.C_AA(:);
%Dcorr6 = -log(corr6);
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^12_2.mat');
corr7=m.avg_fluc.C_AA(:);
Dcorr6 = -log(corr6*0.6+corr7*0.4);
%}

%%Plot the size dependence of the correlation function
figure(1)
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

xlabel('t','FontSize',20)
ylabel('-log|C(tw,tw+t)|','FontSize',20)
ax = gca;
ax.FontSize=16;
lgd=legend('Lx=2^8','Lx=2^9','Lx=2^{10}','Lx=2^{11}','Lx=2^{12}','Location','northwest'); %'Lx=2^7','Lx=2^8','Lx=2^{11}','Lx=2^{12}','Lx=2^{13}','tw=60000600','tw=3200','tw=6400', ,'Lx=2^{10}','Lx=2^{11}','Lx=2^{12}'
lgd.FontSize=18;
axis([10^0 10^4 10^(-7) 10^(-2)])

avg_corr2 = mean(Dcorr2(2000:end));
std_corr2 = std(Dcorr2(2000:end));
%}
avg_corr3 = mean(Dcorr3(2000:end));
std_corr3 = std(Dcorr3(2000:end));

avg_corr4 = mean(Dcorr4(2000:end));
std_corr4 = std(Dcorr4(2000:end));

avg_corr5 = mean(Dcorr5(2000:end));
std_corr5 = std(Dcorr5(2000:end));

avg_corr6 = mean(Dcorr6(2000:end));
std_corr6 = std(Dcorr6(2000:end));

width_s = [avg_corr2 avg_corr3 avg_corr4 avg_corr5 avg_corr6];%avg_corr1 avg_corr2  avg_corr7 avg_corr6
std_s = [std_corr2 std_corr3 std_corr4 std_corr5 std_corr6];%std_corr1 std_corr2 std_corr7 std_corr6


x = log(Lx(:));
y = log(width_s(:));

[p,S] = polyfit(x,y,1);
[yfit,delta] = polyval(p,x,S);

slope = p(1);
alpha = slope/2;

% standard error of slope
n = length(x);
yres = y - polyval(p,x);
s2 = sum(yres.^2)/(n-2);
Sxx = sum((x-mean(x)).^2);
slope_err = sqrt(s2/Sxx);
standard_errors_alpha = slope_err/2;

fprintf('alpha = %.3f +/- %.3f\n', alpha, standard_errors_alpha);
save('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/fit_error.mat','width_s','std_s','alpha','standard_errors_alpha');

%Plot the alpha fitting
figure(2)
L = 2^7:2^13;
chi = L.^p(1)*exp(1)^p(2);
loglog(L,chi,'k','LineWidth',2);
txt = {'$\propto L^{2\alpha}$'};%,P_g(1)/2
text(0.4*10^3,10^(-3),txt,'FontSize',20,'Interpreter','latex');
txt = {'$\alpha=1.35 \pm 0.02$'};%,P_g(1)/2
text(0.4*10^3,5*10^(-4),txt,'FontSize',20,'Interpreter','latex');
hold on;
errorbar(Lx(:),width_s(:),std_s(:),'o','MarkerSize',5,'LineWidth',2,'Color','k');
set(gca, 'XScale','log', 'YScale','log')
%set(gca, 'XScale','log', 'YScale','log')
hold off;
xlabel('L')
ylabel('\chi^s_{AA}(L)')
ax = gca;
ax.FontSize=18;
%title('L^{2\alpha}','FontSize',16)

%
ImageID='/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/alpha.fig';
saveas(gcf,ImageID);

%{
%% Old weighted-fitting method 

% Define the model function
model = @(params, x)  x.^(2*params(1))* exp(1)^params(2);

%%alpha fitting
% Sample data with errors
xdata = Lx(1:5);
ydata = width_s(1:5);
y_err = std_s(1:5);  % Example errors in y

% Define the weighted error function
weighted_error_function = @(params, x, y, y_err) ((model(params, x) - y) ./ y_err).^2;

% Initial guess for the parameters
initial_guess = [1.5, 0];

%Fit the model, calculate the standard error of the fit parameters, use the Jacobian:
[params_alpha, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(params, x) weighted_error_function(params, x, ydata, y_err), initial_guess, xdata, ydata);


% Optimal values for the parameters
disp('Optimal parameters:');
disp(params_alpha);
alpha=params_alpha(1);

%
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

save('/Users/Phantom/Documents/MATLAB/Flocking/Density/NOld_data/jm=-0.25/Grad/fit_error_DA=100.mat','width_s','std_s','alpha','standard_errors_alpha');
%}
