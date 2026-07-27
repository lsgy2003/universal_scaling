clear all;

sigma=0.5;
jp=0.002;

T=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
ts=0:dt*Evo:T;
Lx=[2^8 2^9 2^10 2^11 2^12 2^13]; %2^7 2^11 2^12 2^13
lgLx =log(Lx);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^8.mat');
freq0=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^9.mat');
freq1=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^10.mat');
freq2=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^11.mat');
freq3=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^12.mat');
freq4=m.omega1(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Omega/sigma=0.5/Size/2^13.mat');
freq5=m.omega1(:);


clf;
figure(1)
loglog(ts,freq0,'Color',[100, 149, 237]/255);
hold on;
loglog(ts,freq1,'Color',[255, 191, 0]/255);
hold on;
loglog(ts,freq2,'Color',[255, 127, 80]/255);
hold on;
loglog(ts,freq3,'Color',[204, 204, 255]/255);
hold on;
loglog(ts,freq4,'Color',[159, 226, 191]/255);
hold on;
loglog(ts,freq5,'Color',[0 0 0]);
hold on;
xlabel('$t$','Interpreter','latex')
ylabel('$\langle \overline{\Omega}_{A}(L,t)\rangle$','Interpreter','latex')
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','L=2^{13}','Location','southwest'); %,'L=2^{14}'
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;


%saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_L_sigma=0.5_1.fig');
%exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_L_sigma=0.5_1.eps','ContentType','vector','BackgroundColor','none');


figure(2)
y1=[freq0(end),freq1(end),freq2(end),freq3(end),freq4(end),freq5(end)];
plot(Lx,y1,'b--o','LineWidth',1.5);
xlabel('$L$','Interpreter','latex')
ylabel('$\langle \overline{\Omega}_{A}(L,t=10^4)\rangle$','Interpreter','latex')
ax = gca;
ax.FontSize=18;
xlim([0,8500]);

%saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_L_sigma=0.5_2.fig');
%exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_L_sigma=0.5_2.eps','ContentType','vector','BackgroundColor','none');


figure(3)

avg_freq0 = mean(freq0(2000:end));
std_freq0 = std(freq0(2000:end));

avg_freq1 = mean(freq1(2000:end));
std_freq1 = std(freq1(2000:end));

avg_freq2 = mean(freq2(2000:end));
std_freq2 = std(freq2(2000:end));

avg_freq3 = mean(freq3(2000:end));
std_freq3 = std(freq3(2000:end));

avg_freq4 = mean(freq4(2000:end));
std_freq4 = std(freq4(2000:end));

avg_freq5 = mean(freq5(2000:end));
std_freq5 = std(freq5(2000:end));

width_s = [avg_freq0 avg_freq1 avg_freq2 avg_freq3 avg_freq4 avg_freq5];%avg_freq1 avg_freq7 
std_s = [std_freq0 std_freq1 std_freq2 std_freq3 std_freq4 std_freq5];%std_freq1 std_freq7 


% Define the model function
model = @(params, x)  x.^(params(1))* exp(1)^params(2);

% Sample data with errors
xdata = Lx;
ydata = width_s;
y_err = std_s;  % Example errors in y

% Define the weighted error function
weighted_error_function = @(params, x, y, y_err) ((model(params, x) - y) ./ y_err).^2;
%weighted_error_function = @(params, x, y, y_err) (model(params, x) - y) ./ y_err;


% Initial guess for the parameters
initial_guess = [-0.5, 0];

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
figure(3)
x = 2^7:2^14;
y = x.^(params_alpha(1))*exp(1)^params_alpha(2);
loglog(x,y,'k--','LineWidth',1);
%plot(x,y,'k--','LineWidth',1);
txt = {'$\eta=-0.49 \pm 0.01$'};%,P_g(1)/2
text(0.5*10^3,2*10^(-2),txt,'FontSize',20,'Interpreter','latex');
txt = {'$\propto L^\eta$'};
text(0.5*10^3,1.5*10^(-2),txt,'FontSize',20,'Interpreter','latex');
hold on;
errorbar(Lx(:),width_s(:),std_s(:),'o');
%set(gca, 'XScale','log', 'YScale','log')
hold off;
xlabel('$L$','Interpreter','latex')
ylabel('$\langle \overline{\Omega}_{A}(L)\rangle$','Interpreter','latex')
ax = gca;
ax.FontSize=18;


saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_L_sigma=0.5_3.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/omega_L_sigma=0.5_3.eps','ContentType','vector','BackgroundColor','none');
