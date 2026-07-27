clear all;

T=10000;
T1=50000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 1000 3000 7000];
Dt=1:dt*Evo:(T-tw(4))+1;
Dt1=1:dt*Evo:(T1-tw(4))+1;
n = 2.^(5:9);
%
m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/SI/S3C.mat');
corr2 = m.avg_fluc_N.C_AA(:,1);
Dcorr2 = -log(corr2);

corr3=m.avg_fluc_N.C_AA(:,2);
Dcorr3 = -log(corr3);

corr4=m.avg_fluc_N.C_AA(:,3);
Dcorr4 = -log(corr4);

corr5=m.avg_fluc_N.C_AA(:,4);
Dcorr5 = -log(corr5);
%
corr6=m.avg_fluc_N.C_AA(:,5);
Dcorr6 = -log(corr6);


avg_corr2 = mean(corr2(2000:end));
std_corr2 = std(corr2(2000:end));
%}
avg_corr3 = mean(corr3(2000:end));
std_corr3 = std(corr3(2000:end));

avg_corr4 = mean(corr4(2000:end));
std_corr4 = std(corr4(2000:end));

avg_corr5 = mean(corr5(2000:end));
std_corr5 = std(corr5(2000:end));
%
avg_corr6 = mean(corr6(2000:end));
std_corr6 = std(corr6(2000:end));
%}
width_s = [avg_corr2 avg_corr3 avg_corr4 avg_corr5 avg_corr6];%avg_corr1 avg_corr2  avg_corr7 avg_corr6
std_width_s = [std_corr2 std_corr3 std_corr4 std_corr5 std_corr6];%std_corr1 std_corr2 std_corr7 std_corr6

clf;

%%Plot the New dependence of the correlation function
figure(1)

loglog(Dt, corr2,'.','MarkerSize',10,'Color',[100, 149, 237]/255);
hold on;
loglog(Dt, corr3,'.','MarkerSize',10,'Color',[255, 191, 0]/255);
hold on;
loglog(Dt, corr4,'.','MarkerSize',10,'Color',[255, 127, 80]/255);
hold on;
loglog(Dt, corr5,'.','MarkerSize',10,'Color',[204, 204, 255]/255);
hold on;
loglog(Dt, corr6,'.','MarkerSize',10,'Color',[159, 226, 191]/255);
hold off;
%}

xlabel('t','FontSize',20)
ylabel('C_{AA}(t_0,t_0+t;L)','FontSize',20)
ax = gca;
ax.FontSize=16;
lgd=legend('N=2^6','N=2^7','N=2^8','N=2^9','N=2^{10}','Location','southwest'); 
lgd.FontSize=14;
axis([10^1 3*10^3 0.999 1.0003]);

ImageID='/Users/Phantom/Desktop/Code/data/ensemble_CEP_C_AA_1.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/ensemble_CEP_C_AA_1.eps','ContentType','vector','BackgroundColor','none');


% Define the model function
model = @(params, x)  x.^(params(1))* params(2);

%%alpha fitting
% Sample data with errors
xdata = n;
ydata = width_s;
y_err = std_width_s;  % Example errors in y

% Initial guess for the parameters
initial_guess = [0.5, 1];

%Fit the model, calculate the standard error of the fit parameters, use the Jacobian:
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
x = 20:800;
y=0.999965*x.^0;
%loglog(x,y,'k-','LineWidth',1);
txt = {'$\sim \rm{const.}$'};
text(100,0.99998,txt,'FontSize',25,'Interpreter','latex');
hold on;
errorbar(xdata,width_s(:),std_width_s(:),'square','MarkerSize',10,'Color','k',"MarkerFaceColor",'k');
set(gca, 'XScale','log', 'YScale','log')
hold off;
xlabel('N')
ylabel('C_{AA}^{s}(N)')
ax = gca;
ax.FontSize=16;
axis([20 800 0.99994 1]);
box on;


%
ImageID='/Users/Phantom/Desktop/Code/data/ensemble_CEP_C_AA_2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/ensemble_CEP_C_AA_2.eps','ContentType','vector','BackgroundColor','none');
%}

%save('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/ws_jp=0.0097_1.mat','width_s','std_width_s','alpha','standard_errors_alpha');



