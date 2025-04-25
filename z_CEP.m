%%Data from sigma=0.005, jp=0.0097

clear all;

Lx=[2^9 2^10 2^11 2^12];
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.001;
T=10000;
Evo=1000;
tw=1001;
StopT = T - tw+1;
Nt=StopT/dt;
t=Nt/Evo;
Lx=[2^8 2^9 2^10 2^11 2^12]; %2^7 2^11 2^12 2^13


clf; 
%plot the amplitude spectra of frequency
figure(1)
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^9_tw=1001.mat');
P2=m.P(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^10_tw=1001.mat');
P3=m.P(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^11_tw=1001.mat');
P4=m.P(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^12_tw=1001.mat');
P5=m.P(:);

f= 1/(dt*Evo)*((0:(t/2)))/t;

plot(f(2:end),P2(2:end),'LineWidth',2,'Color',[255, 191, 0]/255) ;
hold on;
plot(f(2:end),P3(2:end),'LineWidth',2,'Color',[255, 127, 80]/255) ;
hold on;
plot(f(2:end),P4(2:end),'LineWidth',2,'Color',[204, 204, 255]/255) ;
hold on;
plot(f(2:end),P5(2:end),'LineWidth',2,'Color',[159, 226, 191]/255) ;
hold on;
plot(peaks,peaks_amp,'square','MarkerSize',8,'MarkerFaceColor','r')
hold on;

axis([0. 0.015 0 0.02]);
xlabel('f')
ylabel('Amplitude')
lgd=legend('L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northeast');
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;

ImageID = '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/fft_scaling_tw=1001_2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/fft_scaling_tw=1001_2.eps','ContentType','vector','BackgroundColor','none');


%% z fitting

peaks=[0.0096 0.0048 0.00244 0.00122]; %sigma=0.005, jp=0.0097
peaks_amp=[0.00003 0.0004 0.0042 0.0190];
peaks_err=[0.0001 0.0001 0.00004 0.00002];

% Define the model function
model = @(params, x)  x.^(-params(1))* exp(1)^params(2);

% Sample data with errors
xdata = Lx;
ydata = peaks;
y_err = peaks_err;  % Example errors in y

% Define the weighted error function
weighted_error_function = @(params, x, y, y_err) (model(params, x) - y) ./ y_err;

% Initial guess for the parameters
initial_guess = [1, 1];

%Fit the model, calculate the standard error of the fit parameters, use the Jacobian:
[params_z, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(params, x) weighted_error_function(params, x, ydata, y_err), initial_guess, xdata, ydata);

% Optimal values for the parameters
disp('Optimal parameters:');
disp(params_z);
z=params_z(1);

% Compute the covariance matrix of the parameters
cov_matrix = inv(J.' * J);

% Extract standard errors from the covariance matrix
standard_errors = sqrt(diag(cov_matrix));

% Display the standard errors
disp('Standard errors of fit parameters:');
disp(standard_errors);

% Convert sparse matrix to full matrix
standard_errors = full(standard_errors);
standard_errors_z = standard_errors(1);

filename='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/fft_peaks_tw=1001.mat';
save(filename,'sigma','jp','Lx','peaks','peaks_amp','peaks_err','z','standard_errors_z');



m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/data/fft_peaks_tw=1001.mat');
peak1 = m.peaks;
std1=m.peaks_err;
z1=m.z;
std_z1=m.standard_errors_z;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/data/fft_peaks_width.mat');
peak2 = m.peaks;
std2=m.peaks_err;
z2=m.z;
std_z2=m.standard_errors_z;


figure(2)

errorbar(Lx(2:end), peak1, std1, 'o','MarkerSize',5,'LineWidth',1,'Color',[0 0 0]);
set(gca, 'XScale','log', 'YScale','log')
hold on;
errorbar(Lx,peak2, std2, 'diamond','MarkerSize',5,'LineWidth',1,'Color',[0.6350 0.0780 0.1840]);
hold on;

%%fit the full correlator
P_g=polyfit(log(Lx(2:end)),log(peak1(:)),1);
fprintf('%d\n',P_g);
x = 1.5*(2^8:2^12);
y = x.^P_g(1)*exp(1)^P_g(2);

%plot
loglog(x,y,'k','LineWidth',2);
hold on;
txt = {'$z=0.99\pm0.01$'};
text(0.9*10^2,1.5*10^(-2),txt,'Interpreter','latex','FontSize',18);
txt = {'$f_{-\log|C_{AA}|} \propto L^{-z}$'};
text(0.9*10^2,1*10^(-2),txt,'Interpreter','latex','FontSize',18);
txt={'$f\propto k^z \propto L^{-z}$'};
text(1.3*10^3,3*10^(-2),txt,'Interpreter','latex','FontSize',20,'FontWeight','bold');
hold on;
xlabel('L')
ylabel('f_{-log|C_{AA}|}(L), f_{w_{A}}(L)')
ax = gca;
ax.FontSize=18;
hold on;

%fit the phase correlator
P_g=polyfit(log(Lx(:)),log(peak2(:)),1);
fprintf('%d\n',P_g);
x = 1.5*(2^7:2^12);
y = x.^P_g(1)*exp(1)^P_g(2);
plot(x,y,'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
txt = {'$z=0.98\pm0.01$'};
text(1.1*10^3,1.5*10^(-2),txt,'Interpreter','latex','FontSize',18,'Color',[0.6350 0.0780 0.1840]);
txt = {'$f_{w_{A}} \propto L^{-z}$'};
text(1.1*10^3,1.1*10^(-2),txt,'Interpreter','latex','FontSize',18,'Color',[0.6350 0.0780 0.1840]);


hold off;


ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/z_CEP1.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/z_CEP1.eps','ContentType','vector','BackgroundColor','none');

