%n = 17; %length(jp)


    m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/random_initial/xi=0/D=0/avg.mat');
    jp0 = m.jp;
    y0 =  m.avg_freq;


    m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/random_initial/xi=0/D=1/avg.mat');
    jp1 = m.jp_all;
    y1 =  m.avg_freq;



%Compare the frequency
figure(2)
plot(jp0(1:end-1),y0(1:end-1),'k-^','LineWidth',2);
hold on;
plot(jp1(1:end-1),y1(1:end-1),'r--*','LineWidth',2);
hold on;
xlabel('$j_+$','Interpreter','latex')
ylabel('$\overline{\Omega}_A$','Interpreter','latex')
ax = gca;
ax.FontSize=16;
lgd=legend('D_A=D_B=0','D_A=D_B=1','Location','northeast');
lgd.FontSize = 14;
hold off;

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/diffusion_freq.fig';
saveas(gcf,ImageID);
%{
%Compare the kink numbers
figure(3)
plot(jp(:),z0(:),'b-^','LineWidth',2);
hold on;
plot(jp(:),z1(:),'--*','LineWidth',3,'Color',colors_p(1,:));
hold on;
plot(jp(:), z2(1:18),'--*','LineWidth',2,'Color',colors_p(2,:));
hold on;
plot(jp(:), z3(1:18),'--*','LineWidth',2,'Color',colors_p(3,:));
hold on;
xlabel('j+','FontSize',14)
ylabel('N','FontSize',14)
lgd=legend('D=0','D=1','D=5','D=10','Location','southeast');
lgd.FontSize = 14;
hold off;

ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/random_initial/xi=0/diffusion_kink.fig';
saveas(gcf,ImageID);
%}