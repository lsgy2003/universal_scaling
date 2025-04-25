clear all;

sigma=0.005;
jp=0.0097;

%temporal correlation functions
Lx=[2^8 2^9 2^10 2^11 2^12]; 
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.001;
T=10000;
Evo=1000;
Nt=T/dt;
t=Nt/Evo;


%%z fitting


%
peaks=[0.0360 0.0190 0.0096 0.0048 0.0024]; %width jp=0.0097
peaks_err=[0.001 0.001 0.0001 0.0001 0.0001]; 
peaks_amp = [1.48*10^(-6) 2.21*10^(-5) 0.0003 0.0031 0.019];
%}

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

%filename='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/fft_peaks_width.mat';
%save(filename,'sigma','jp','Lx','peaks','peaks_amp','peaks_err','z','standard_errors_z');

clf;
figure(1)

x = 2^7:2^13;
y = x.^(-params_z(1))*exp(1)^params_z(2);
loglog(x,y,'k--','LineWidth',1);
txt = {'$z=0.975$'};%,P_g(1)
text(10^3,2*10^(-2),txt,'Interpreter','latex','FontSize',30);
txt2={'$f\sim L^{-z} \sim k^z$'};
text(10^3,4*10^(-2),txt2,'Interpreter','latex','FontSize',30);%'Interpreter','latex',
hold on;
errorbar(Lx,peaks,peaks_err,'square','MarkerSize',10,'MarkerFaceColor','r');%'MarkerFaceColor','r','MarkerSize',10
hold off;
xlabel('L')
ylabel('$f_{w_{A}}$','Interpreter','latex')
ax = gca;
ax.FontSize=30;
%axis([10^2 10^4 2*10^(-4) 5*10^(-2)]);

%{
ImageID = '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/fit_fft_width_1.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/fit_fft_width_1.eps','ContentType','vector','BackgroundColor','none');
%}
figure(2)


m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^8_width.mat');
P1=m.P(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^9_width.mat');
P2=m.P(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^10_width.mat');
P3=m.P(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^11_width.mat');
P4=m.P(:);

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^12_width.mat');
P5=m.P(:);

f= 1/(dt*Evo)*((0:(t/2)))/t;
plot(f(2:end),P1(2:end),'LineWidth',2,'Color',[100, 149, 237]/255) ;
hold on;
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

axis([0. 0.04 0 0.025]);
xlabel('$f_{w_{A}}$','Interpreter','latex')
ylabel('Amplitude')
lgd=legend('L=2^8','L=2^9','L=2^{10}','L=2^{11}','L=2^{12}','Location','northeast'); 
lgd.FontSize=16;
ax = gca;
ax.FontSize=18;
hold off;
%
ImageID = '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/fft_scaling_width_2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/fft_scaling_width_2.eps','ContentType','vector','BackgroundColor','none');
%}


