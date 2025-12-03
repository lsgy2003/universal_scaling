%%Fit the peak

clear all;
%temporal correlation functions
Lx=2^8;
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

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^8_tw=1001.mat');
P1=m.P(:);

clf;
figure(1)
f= 1/(dt*Evo)*((0:(t/2)))/t;
plot(f(2:end),P1(2:end),'LineWidth',1,'Color',[61,101,165]/255);
axis([0. 0.04 0 Inf]);
xlabel('f')
ylabel('|P(f)|')
ax = gca;
ax.FontSize=20;
hold on;

%choose the peak region to fit
f_peak=f(23:65)';
P_peak=P1(23:65);


%fit to a Gassuian peak + an exp decay background
F=@(c,xdata)c(3)/(c(1)*sqrt(2*pi))*exp(-1/2*((xdata-c(2))/c(1)).^2)+c(4)*exp(-c(5)*xdata)+c(6);
c0 = [0.0005 0.0049 0.000 0.000 150.0 0.000];
[c,resnorm,~,exitflag,output] = lsqcurvefit(F,c0,f_peak,P_peak);
x1 = linspace(0.002,0.008,200);
disp(c);

plot(x1,F(c,x1),'r','LineWidth',2);
eqn = sprintf('$f_{peak} = %.4f$',c(2));
text(0.01,3*10^(-4),eqn,'Interpreter','latex','FontSize',20);
hold off;

%ImageID = '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/jp=0.0097/FFT/2^8_fft.fig';
%saveas(gcf,ImageID);
%exportgraphics(gcf,'2^8_fft.eps','ContentType','vector','BackgroundColor','none');



