%% Robustness of CEP exponents
clear; clc; close all;

%% --- Data ---
ref1 = 1.36;   %alpha_C_{AA}, sigma=0.001
ref1_err = 0.03;

ref2 = 1.35;   %alpha_w_A, sigma=0.001
ref2_err = 0.03;

labels = { ...
    '$4 \times 10^{-6}$', ...
    '$2.5 \times 10^{-5}$', ...
    '$75$', ...
    '$50$', ...
    '$-0.35$', ...
    '$-0.15$', ...   
    '+ Cross diffusion', ...
    '+ Gradient nonlinearity'
    };

% alpha_C and error
alpha_Cp  = [1.35 1.35  1.36 1.32  1.38 1.39  1.36 1.35]; %Extracted from data files under data/jm=-0.35, data/jm=-0.15, and data/jm=-0.25/Extension
err_Cp    = [0.02 0.02  0.04 0.06  0.04 0.01  0.03 0.03];

% alpha_w and error
alpha_w   = [1.37 1.35  1.39 1.35  1.42 1.40  1.35 1.35];%Extracted from data files under data/jm=-0.35, data/jm=-0.15, and data/jm=-0.25/Extension
err_w     = [0.02 0.03  0.03 0.04  0.02 0.01  0.03 0.03];

group_id = [1 1  2 2  3 3  4 4];

colors = [
    0.000 0.447 0.741;  % blue
    0.850 0.325 0.098;  % orange
    0.466 0.674 0.188;  % green
    0.494 0.184 0.556;  % purple
];

markers = {'o','s','d','^','v'};

%% --- Y positions with gaps between groups ---
y = zeros(size(labels));

pos = 1;
for i = 1:numel(labels)
    y(i) = pos;
    pos = pos + 1;

    if i < numel(labels) && group_id(i+1) ~= group_id(i)
        pos = pos + 0.8;
    end
end

%% --- Figure ---
figure('Color','w','Position',[100 100 600 1200]);

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

plot_panel(alpha_Cp, err_Cp, y, labels, group_id, colors, markers, ...
    ref1, ref1_err, '$\alpha_{C_{AA}}$');

ax = gca;

ax.Position = [0.3 0.12 0.65 0.82];

text(1.19, 0.2, '\textbf{Noise Strength}  ', 'FontWeight','bold', 'FontSize',18, 'HorizontalAlignment','right','Color',[0.000 0.447 0.741],'Interpreter','latex');
text(1.19, 0.5, '(ref. $\sigma=10^{-6}$)', 'FontWeight','bold', 'FontSize',16, 'HorizontalAlignment','right','Color',[0.000 0.447 0.741],'Interpreter','latex');


text(1.19, 2.9, '\textbf{Diffusion}  ', 'FontWeight','bold', 'FontSize',18, 'HorizontalAlignment','right','Color',[0.850 0.325 0.098],'Interpreter','latex');
text(1.19, 3.2, '(ref. $D_A=100$)', 'FontWeight','bold', 'FontSize',16, 'HorizontalAlignment','right','Color',[0.850 0.325 0.098],'Interpreter','latex');

text(1.19, 5.7, '\textbf{Nonreciprocity}  ', 'FontWeight','bold', 'FontSize',18, 'HorizontalAlignment','right','Color',[0.466 0.674 0.188],'Interpreter','latex');
text(1.19, 6.0, '(ref. $j_- = -0.25$)', 'FontWeight','bold', 'FontSize',16, 'HorizontalAlignment','right','Color',[0.466 0.674 0.188],'Interpreter','latex');


text(1.19, 8.5, '\textbf{Model Extension}  ', 'FontWeight','bold', 'FontSize',18, 'HorizontalAlignment','right','Color',[0.494 0.184 0.556],'Interpreter','latex');


saveas(gcf,'/Users/Phantom/Desktop/Code/data/CEP_robust1.fig');
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/CEP_robust1.eps','ContentType','vector','BackgroundColor','none');


figure('Color','w','Position',[100 100 600 1200]);


tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
plot_panel(alpha_w, err_w, y, labels, group_id, colors, markers, ...
    ref2, ref2_err, '$\alpha_{w_A}$');
yticklabels([])

saveas(gcf,'/Users/Phantom/Desktop/Code/data/CEP_robust2.fig');
exportgraphics(gcf,'/Users/Phantom/Desktop/Code/data/CEP_robust2.eps','ContentType','vector','BackgroundColor','none');



%% --- Helper function ---
%% --- Helper function ---
function plot_panel(alpha, err, y, labels, group_id, colors, markers, ref, ref_err, xlab)

    hold on;

    ymin = min(y)-0.8;
    ymax = max(y)+0.8;

    % reference band
    patch([ref-ref_err ref+ref_err ref+ref_err ref-ref_err], ...
          [ymin ymin ymax ymax], ...
          [0.85 0.85 0.85], ...
          'EdgeColor','none','FaceAlpha',0.35);

    % reference line
    xline(ref,'k-','LineWidth',1.19);

    % data points with horizontal error bars
    for i = 1:numel(alpha)
        gid = group_id(i);

        errorbar(alpha(i), y(i), err(i), ...
            'horizontal', ...
            'LineStyle','none', ...
            'Color',colors(gid,:), ...
            'Marker',markers{gid}, ...
            'MarkerSize',8, ...
            'MarkerFaceColor',colors(gid,:), ...
            'MarkerEdgeColor',colors(gid,:), ...
            'LineWidth',1.3, ...
            'CapSize',8);
    end

    % group separators
    for i = 1:numel(group_id)-1
        if group_id(i+1) ~= group_id(i)
            ysep = 0.5*(y(i)+y(i+1));
            yline(ysep,'--','Color',[0.75 0.75 0.75],'LineWidth',1);
        end
    end

    xlim([1.190 1.50]);
    ylim([ymin ymax]);

    yticks(y);
    yticklabels(labels);

    xlabel(xlab,'Interpreter','latex','FontSize',20);
    set(gca,'TickLabelInterpreter','latex','FontSize',16);
    set(gca,'YDir','reverse')
    set(gca,'XAxisLocation','top')
    box on;

end