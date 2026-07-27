clear all;

Lx=[2^8 2^9 2^10 2^11 2^12]; 
lgLx =log(Lx);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/EW/fit_error.mat');
width1 = m.width_s;
std1=m.std_s;


clf;
%
%%Fitting \alpha
figure(1);
errorbar(lgLx, width1, std1, 'o','MarkerSize',5,'LineWidth',1,'Color',[0 0.4470 0.7410]);
hold on;

%%fit to log scaling
P_g=polyfit(lgLx(2:end),width1(2:end),1);
fprintf('%d\n',P_g);
x = linspace(5,9,100);
y = P_g(1)*x+P_g(2);
%chi-square goodness test
dw = (P_g(1)*lgLx+P_g(2))-width1;
chi_log= sum((dw./std1).^2)/(length(dw)-2);
%R^2 goodness test
dw2= (P_g(1)*lgLx+P_g(2))-mean(width1);
R_log = 1-sum(dw.^2)/sum(dw2.^2);
%plot
plot(x,y,'k--','LineWidth',1);
hold on;
txt={'$\propto 2\gamma \log L$'};
text(8,1*10^(-4),txt,'Interpreter','latex','FontSize',20);%'Interpreter','latex',
txt = {'$\gamma=0.0001$'};%,P_g(1)/2
text(8,3*10^(-5),txt,'Interpreter','latex','FontSize',18);
hold on;
xlabel('log(L)')%'Interpreter','latex'
ylabel('\chi_{AA}^s(L)')
ax = gca;
ax.FontSize=16;
%axis([3 7.5 -1 2.5])


%fit to power law
%P_g=polyfit(lgLx(2:end),log(width1(2:end)),1);
%fprintf('%d\n',P_g);
x = linspace(5,9,100);
P_g(1)=m.alpha*2;
P_g(2)=m.offset;
y = exp(P_g(1)*x+P_g(2));
plot(x,y,'LineWidth',2,'Color',[0 0.4470 0.7410]);
hold on;
txt={'$\propto L^{2\alpha}$'};
text(6,1.5*10^(-4),txt,'Interpreter','latex','FontSize',20,'Color',[0 0.4470 0.7410]);%'Interpreter','latex',
txt={'$\alpha=\frac{1}{2}$'};
text(6,1*10^(-4),txt,'Interpreter','latex','FontSize',20,'Color',[0 0.4470 0.7410]);%'Interpreter','latex',

%
ImageID='/Users/Phantom/Desktop/Code/data/power_scaling3.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/power_scaling3.eps','ContentType','Vector','BackgroundColor','None');

figure(3)
% Create a new figure for the legend
legendFig = figure(3);
axes('Position', [0 0 1 1], 'Visible', 'off'); % Invisible axes

% Plot dummy lines for the legend
hold on;
hLeg(1) = plot(nan, nan, 'o', 'LineWidth',2, 'Color', [0 0.4470 0.7410]);
%hLeg(2) = plot(nan, nan, '--', 'LineWidth',1, 'Color', [0, 0, 0]);


% Create the legend
    lgd = legend(hLeg, {...
    '$\sigma=0.005,j_+ = 0.0400$'});
set(lgd, 'FontSize', 18, 'Location', 'NorthEast','Interpreter','latex');
%lgd.Position = [0.1, 0.1, 0.1, 0.3]; % Adjust size and position of legend within the new figure

% Save or display the legend figure
% Or use savefig for a .fig file:

savefig(legendFig, '/Users/Phantom/Desktop/Code/data/power_scaling3_Legend.fig');
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/power_scaling3_Legend.eps','ContentType','vector','BackgroundColor','none');
%}
