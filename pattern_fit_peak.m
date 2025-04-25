%%Fit the peak

clear all;
%temporal correlation functions
Lx=2^8;
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

dt = 0.01;
T=800;
Evo=100;
Nt = T/dt;
t=Nt/Evo;
f= 1/(dt*Evo)*((0:(t/2)))/t;

m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/lin_stability/2^8_freq.mat');

clf;
figure(1)

Y = fft(m.Output.density{1:1}(2:end));
P2 = abs(Y/t);
P1 = P2(1:t/2+1); %single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
plot(f(2:end),P1(2:end),'LineWidth',1,'Color',[61,101,165]/255);
axis([0. 0.3 0 0.012]);
xlabel('f')
ylabel('Amplitude')
hold on;
errorbar(0.07125,0.011,0.0001,'square','MarkerSize',10,'MarkerFaceColor','r');
hold on;
txt={'$f_{\overline{|\vec{P}_A|}}=0.0713$'};
text(0.1,0.008,txt,'Interpreter','latex','FontSize',30);%'Interpreter','latex',
ax = gca;
ax.FontSize=18;
hold off;
%
ImageID = '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/density_freq.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/density_freq.eps','ContentType','Vector','BackgroundColor','None');
%}

figure(2)
Y = fft(m.Output.phase_dif{1,1}(2:end));
P2 = abs(Y/t);
P1 = P2(1:t/2+1); %single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
plot(f(2:end),P1(2:end),'LineWidth',1,'Color',[61,101,165]/255);
axis([0. 0.3 0 0.07]);
xlabel('f')
ylabel('Amplitude')
hold on;
errorbar(0.07125,0.058,0.0001,'square','MarkerSize',10,'MarkerFaceColor','r');
hold on;
txt={'$f_{\overline{|\Delta \theta|}}=0.0713$'};
text(0.1,0.04,txt,'Interpreter','latex','FontSize',30);%'Interpreter','latex',
ax = gca;
ax.FontSize=18;
hold off;

%
ImageID = '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/phase_dif_freq.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/phase_dif_freq.eps','ContentType','Vector','BackgroundColor','None');
%}


