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
tw=1000;
StopT = T - tw+1;
Nt=StopT/dt;
t=Nt/Evo;

f= 1/(dt*Evo)*((0:(t/2)))/t;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^8_fft.mat');
P=m.avg_fluc.P_AA(:);

clf;
figure(1)
plot(f,P(:),'LineWidth',1,'Color',[61,101,165]/255);
axis([0. 0.04 0 Inf]);
xlabel('f')
ylabel('|P(f)|')
ax = gca;
ax.FontSize=20;
hold on;

%choose the peak region to fit
f_peak=f(90:300);
P_peak=P(90:300)';


%fit to a Gassuian peak + an exp decay background
F=@(c,xdata)c(3)/(c(1)*sqrt(2*pi))*exp(-1/2*((xdata-c(2))/c(1)).^2)+c(4)*exp(-c(5)*xdata)+c(6);
c0 = [0.004 0.0190 6*10^(-6) 0.000 0.0 -2*10^(-6)];
[c,resnorm,~,exitflag,output] = lsqcurvefit(F,c0,f_peak,P_peak);
x1 = linspace(0.01,0.03,200);
disp(c);

plot(x1,F(c,x1),'r','LineWidth',2);
eqn = sprintf('$f_{peak} = %.4f$',c(2));
text(0.01,3*10^(-4),eqn,'Interpreter','latex','FontSize',20);
hold off;

ImageID = '/Users/Phantom/Desktop/Code/data/jm=-0.25/CEP/2^8_fft.fig';
%saveas(gcf,ImageID);
%exportgraphics(gcf,'2^8_fft.eps','ContentType','vector','BackgroundColor','none');



