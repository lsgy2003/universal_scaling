clear all;

Lx=[2^8 2^9 2^10 2^11 2^12 2^13]; 
lgLx =log(Lx);

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral1/fit_error.mat');
width1 = m.width_s;
std1=m.std_width_s;


m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral2/fit_error.mat');
width2 = m.width_s;
std2=m.std_width_s;

m=load('/Users/Phantom/Desktop/Code/data/jm=-0.25/Chiral3/fit_error.mat');
width3 = m.width_s;
std3=m.std_width_s;


clf;

%%Fitting \alpha
figure(1);
errorbar(lgLx(2:end), width1(2:end), std1(2:end), 'diamond','MarkerSize',5,'LineWidth',1,'Color','k');
hold on;
P_g=polyfit(lgLx(2:end),width1(2:end),1);
fprintf('%d\n',P_g);
x = linspace(5,9.5,100);
y = 0.5*x-2.30;
plot(x,y,'LineWidth',2,'Color','k');
hold on;


errorbar(lgLx(3:end), width2(3:end), std2(3:end), 'diamond','MarkerSize',5,'LineWidth',1,'Color', [255, 127, 80]/255);
hold on;
P_g=polyfit(lgLx(3:end),width2(3:end),1);
fprintf('%d\n',P_g);
x = linspace(6,9.5,100);
y = 0.5*x-2.36;
%plot
plot(x,y,'LineWidth',2,'Color',[255, 127, 80]/255);
hold on;
xlabel('log(L)')
ylabel('\chi_{AA}^s(L)')
ax = gca;
ax.FontSize=18;
%}

%
errorbar(lgLx(:), width3(:), std3(:), 'diamond','MarkerSize',5,'LineWidth',1,'Color',[85, 107, 47]/255);
hold on;
P_g=polyfit(lgLx(:),width3(:),1);
fprintf('%d\n',P_g);
x = linspace(5,9.5,100);
y = 0.5*x-2.64;
%plot
plot(x,y,'LineWidth',2,'Color',[85, 107, 47]/255); 
hold on;

txt={'$\propto 2\gamma \log L$'};
text(5.5,2,txt,'Interpreter','latex','FontSize',22,'Color','k');%'Interpreter','latex',
hold on;
txt = {'$\gamma=\frac{1}{4}$'};%,P_g(1)/2
text(5.5,1.5,txt,'Interpreter','latex','FontSize',18,'Color','k');


%fit to power law
P_g=polyfit(lgLx(2:end),log(width1(2:end)),1);
fprintf('%d\n',P_g);
x = linspace(5,9.5,100);
y = exp(P_g(1)*x+P_g(2));
plot(x,y,'--','LineWidth',1,'Color',[0 0.4470 0.7410]);
hold on;
txt={'$\propto L^{2\alpha}$'};
text(8,3,txt,'Interpreter','latex','FontSize',20,'Color',[0 0.4470 0.7410]);%'Interpreter','latex',
txt={'$\alpha=0.18$'};
text(8,2.5,txt,'Interpreter','latex','FontSize',18,'Color',[0 0.4470 0.7410]);%'Interpreter','latex',
axis([5,9.5,-0.5,3.5])
%}
hold off;
%
ImageID='/Users/Phantom/Desktop/Code/data/log_scaling.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/log_scaling.eps','ContentType','Vector','BackgroundColor','None');
%}

%
figure(2)
% Create a new figure for the legend
legendFig = figure(2);
axes('Position', [0 0 1 1], 'Visible', 'off'); % Invisible axes

% Plot dummy lines for the legend
hold on;
hLeg(1) = plot(nan, nan, 'diamond', 'LineWidth',2, 'Color', 'k');
hLeg(2) = plot(nan, nan, 'diamond', 'LineWidth',2, 'Color', [255, 127, 80]/255);
hLeg(3) = plot(nan, nan, 'diamond', 'LineWidth',2, 'Color', [85, 107, 47]/255);

% Create the legend
    lgd = legend(hLeg, {
    '$\sigma=2.5\!\times\!10^{-1},j_+=0.0097$', ...
    '$\sigma=2.25\!\times\!10^{-2},j_+=0.0140$', ...
    '$\sigma=4.0\!\times\!10^{-2},j_+=0.0080$', ...
    }); %\langle{\rm Var}(\Delta \theta_{AB}(L)\rangle$' '$\sigma=0.1,\frac{D_A}{D_B}=1,j_+=0.040$', ...

set(lgd, 'FontSize', 16, 'Location', 'NorthEast','Interpreter','latex');
%lgd.Position = [0.1, 0.1, 0.1, 0.3]; % Adjust size and position of legend within the new figure

% Save or display the legend figure
% Or use savefig for a .fig file:

savefig(legendFig, '/Users/Phantom/Desktop/Code/data/log_scaling1_Legend.fig');
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/log_scaling1_Legend.eps','ContentType','vector','BackgroundColor','none');
%}
