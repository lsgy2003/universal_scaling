clear all;
clf;


figure(1)
%%sigma=0.25
jp1 = [0.120 0.150 0.200 0.250 0.300 0.350];
beta1 = [0.513 0.419 0.327 0.282 0.241 0.230];
dbeta1 = [0.021 0.012 0.023 0.023 0.004 0.009];

x1 = linspace(0.1,0.4,500);
F=@(c,xdata)c(1)*exp(-c(2)*xdata)+c(3);
c0 = [1 10 0.25];
[c,resnorm,~,exitflag,output] = lsqcurvefit(F,c0,jp1,beta1);
hold on;
semilogx(x1,F(c,x1),'r','LineWidth',1,'Color','[1.0000 0.8000 0.2000]');
eqn = sprintf('$\\beta = %.3fe^{-%.3fLx}+%.3f$',c(1),c(2),c(3));
%text(0.2,0.6,eqn,'Interpreter','latex','FontSize',16);
hold on;
%}


%%sigma=0.2
jp2 = [0.100 0.120 0.150 0.200 0.250 0.300 0.350];
beta2 = [0.520 0.410 0.342 0.271 0.239 0.237 0.235]; 
dbeta2 = [0.011 0.011 0.015 0.010 0.005 0.006 0.006]; 

x2 = linspace(0.06,0.4,500);
F=@(c,xdata)c(1)*exp(-c(2)*xdata)+c(3);
c0 = [5 20 0.25];
[c,resnorm,~,exitflag,output] = lsqcurvefit(F,c0,jp2,beta2);
hold on;
semilogx(x2,F(c,x2),'g','LineWidth',1,'Color','[0.9500 0.6500 0.0500]');
eqn = sprintf('$\\beta = %.3fe^{-%.3fLx}+%.3f$',c(1),c(2),c(3));
%text(0.2,0.6,eqn,'Interpreter','latex','FontSize',16);
hold on;

%%sigma=0.15
jp3 = [0.07 0.10 0.12 0.14 0.15 0.20 0.25 0.30 0.35];
beta3 = [0.574 0.383 0.296 0.277 0.248 0.232 0.232 0.229 0.228];
dbeta3 = [0.024 0.008 0.011 0.004 0.008 0.006 0.004 0.002 0.005];

x3 = linspace(0.04,0.4,500);
F=@(c,xdata)c(1)*exp(-c(2)*xdata)+c(3);
c0 = [5 50 0.25];
[c,resnorm,~,exitflag,output] = lsqcurvefit(F,c0,jp3,beta3);
hold on;
semilogx(x3,F(c,x3),'b','LineWidth',1,'Color','[0.8000 0.5000 0.0000]');
eqn = sprintf('$\\beta = %.3fe^{-%.3fLx}+%.3f$',c(1),c(2),c(3));
%text(0.2,0.6,eqn,'Interpreter','latex','FontSize',16);
hold on;

%%sigma=0.1
jp4 = [0.040 0.050 0.060 0.070 0.100 0.150 0.200 0.250 0.300 0.350];
beta4 = [0.728 0.534 0.418 0.371 0.270 0.254 0.246 0.242 0.242 0.243];
dbeta4 = [0.036 0.028 0.027 0.022 0.008 0.005 0.003 0.003 0.003 0.005];

%emperical fit
x4 = linspace(0.030,0.4,500);
F=@(c,xdata)c(1)*exp(-c(2)*xdata)+c(3);
c0 = [5 50 0.25];
[c,resnorm,~,exitflag,output] = lsqcurvefit(F,c0,jp4,beta4);
hold on;
semilogx(x4,F(c,x4),'c','LineWidth',1,'Color','[0.6000 0.3500 0.0000]');
eqn = sprintf('$\\beta = %.3fe^{-%.3fLx}+%.3f$',c(1),c(2),c(3));
%text(0.2,0.6,eqn,'Interpreter','latex','FontSize',16);
hold on;

%%sigma=0.05
jp5 = [0.030 0.040 0.050 0.070 0.100 0.150 0.200 0.250 0.300 0.350];
beta5 = [0.429 0.348 0.311 0.289 0.272 0.256 0.249 0.243 0.240 0.221];
dbeta5 = [0.001 0.001 0.001 0.001 0.006 0.003 0.002 0.002 0.001 0.009];

x5 = linspace(0.010,0.4,500);
F=@(c,xdata)c(1)*exp(-c(2)*xdata)+c(3);
c0 = [1 50 0.25];
[c,resnorm,~,exitflag,output] = lsqcurvefit(F,c0,jp5,beta5);
hold on;
semilogx(x5,F(c,x5),'m','LineWidth',1,'Color','[0.4000 0.2500 0.0000]');
eqn = sprintf('$\\beta = %.3fe^{-%.3fLx}+%.3f$',c(1),c(2),c(3));
%text(0.2,0.6,eqn,'Interpreter','latex','FontSize',16);
hold on;


filename='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Xover/corr_beta.mat';
save(filename,'jp1','beta1','jp2','beta2','jp3','beta3','jp4','beta4','jp5','beta5'); %,'jp2','beta2','jp3','beta3'\


errorbar(jp1, beta1, dbeta1, '*', 'MarkerSize',8, 'LineWidth',1, 'Color', [1.0000 0.8000 0.2000])
hold on;
errorbar(jp2, beta2, dbeta2, 'o', 'MarkerSize',8, 'LineWidth',1, 'Color', [0.9500 0.6500 0.0500])
hold on;
errorbar(jp3, beta3, dbeta3, '+', 'MarkerSize',8, 'LineWidth',1, 'Color', [0.8000 0.5000 0.0000])
hold on;
errorbar(jp4, beta4, dbeta4, 'x', 'MarkerSize',8, 'LineWidth',1, 'Color', [0.6000 0.3500 0.0000])
hold on;
errorbar(jp5, beta5, dbeta5, 's', 'MarkerSize',8, 'LineWidth',1, 'Color', [0.4000 0.2500 0.0000])
hold on;

%Comparison to EW scaling 
x = linspace(0,0.41,500);
y = 0.25*(ones(1,length(x)));
loglog(x,y,'--','LineWidth',1,'Color','[0.5290 0.6940 0.1250]')
txt = {' \beta=1/4'};
text(0.03,0.2,txt,'FontSize',18,'Color','[0.5290 0.6940 0.1250]');
hold off;

xlabel('j_+')
ylabel('\beta')
ax = gca;
ax.FontSize=16;
axis([0 0.41 0.1 1])
%lgd=legend('\sigma=0.25','\sigma=0.20','\sigma=0.15','\sigma=0.10','\sigma=0.05','Location','northeast'); %'tw=1600','tw=3200','tw=6400',
%lgd.FontSize=18;
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_beta.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_beta.eps','ContentType','vector','BackgroundColor','none');


figure(2)
% Create a new figure for the legend
legendFig = figure(2);
axes('Position', [0 0 1 1], 'Visible', 'off'); % Invisible axes

% Plot dummy lines for the legend
hold on;
hLeg(1) = plot(nan, nan, '*', 'LineWidth',2, 'Color', [1.0000 0.8000 0.2000]);
hLeg(2) = plot(nan, nan, 'o', 'LineWidth',2, 'Color', [0.9500 0.6500 0.0500]);
hLeg(3) = plot(nan, nan, '+', 'LineWidth',2, 'Color', [0.8000 0.5000 0.0000]);
hLeg(4) = plot(nan, nan, 'x', 'LineWidth',2, 'Color', [0.6000 0.3500 0.0000]);
hLeg(5) = plot(nan, nan, 's', 'LineWidth',2, 'Color', [0.4000 0.2500 0.0000]);



% Create the legend
    lgd = legend(hLeg, {
    '$\sigma=0.25$','$\sigma=0.20$','$\sigma=0.15$','$\sigma=0.10$','$\sigma=0.05$'
    }); %\langle{\rm Var}(\Delta \theta_{AB}(L)\rangle$' '$\sigma=0.1,\frac{D_A}{D_B}=1,j_+=0.040$', ...

set(lgd, 'FontSize', 18, 'Location', 'NorthEast','Interpreter','latex');
%lgd.Position = [0.1, 0.1, 0.1, 0.3]; % Adjust size and position of legend within the new figure

% Save or display the legend figure
% Or use savefig for a .fig file:
savefig(legendFig, '/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_beta_Legend.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_beta_Legend.eps','ContentType','vector','BackgroundColor','none');
