%%Pattern formation phase diagram

clear all;

Lx=2^10;
dx=1/2;
x = (-Lx/2:dx:Lx/2)';
Nx=Lx/dx+1;

clf;

%%Plot the formation of pattern vs. pattern
f=figure(1);
f.Position =[500 500 500 500];
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/lin_stability/no_pattern.mat');

subplot(5,1,1);
plot(x,m.Output.Freq{1,2}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -0.3, 0.3])
title('$t=2$','Interpreter','latex','FontSize',15)
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)

subplot(5,1,2);
plot(x,m.Output.Freq{1,5}{1,1},'k-','LineWidth',2);
title('$t=5$','Interpreter','latex','FontSize',15)
axis([-Lx/2, Lx/2, -0.3, 0.3])
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)

subplot(5,1,3)
plot(x,m.Output.Freq{1,10}{1,1},'k-','LineWidth',2);
title('$t=10$','Interpreter','latex','FontSize',15)
axis([-Lx/2, Lx/2, -0.3, 0.3])
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)

subplot(5,1,4)
plot(x,m.Output.Freq{1,15}{1,1},'k-','LineWidth',2);
title('$t=15$','Interpreter','latex','FontSize',15)
axis([-Lx/2, Lx/2, -0.3, 0.3])
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)

subplot(5,1,5)
plot(x,m.Output.Freq{1,200}{1,1},'k-','LineWidth',2);
title('$t=200$','Interpreter','latex','FontSize',15)
axis([-Lx/2, Lx/2, -0.3, 0.3])
ylabel('\Omega_A','FontSize',15)
xlabel('x','FontSize',20)

%saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/no_pattern.fig');
%exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/no_pattern.eps','ContentType','Vector','BackgroundColor','None');


h=figure(2);
h.Position =[500 500 500 500];
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/lin_stability/pattern.mat');
subplot(5,1,1)
plot(x,m.Output.Freq{1,2}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -0.3, 0.3])
set(gca,'xtick',[])
title('$t=2$','Interpreter','latex','FontSize',15)
ylabel('\Omega_A','FontSize',15)

subplot(5,1,2)
plot(x,m.Output.Freq{1,400}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -0.3, 0.3])
set(gca,'xtick',[])
title('$t=400$','Interpreter','latex','FontSize',15)
ylabel('\Omega_A','FontSize',15)

subplot(5,1,3)
plot(x,m.Output.Freq{1,600}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -1, 1])
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)
title('$t=600$','Interpreter','latex','FontSize',15)

subplot(5,1,4)
plot(x,m.Output.Freq{1,1000}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -1, 1])
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)
title('$t=1000$','Interpreter','latex','FontSize',15)

subplot(5,1,5)
plot(x,m.Output.Freq{1,2000}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -1, 1])
ylabel('\Omega_A','FontSize',15)
xlabel('x','FontSize',20)
title('$t=2000$','Interpreter','latex','FontSize',15)

%saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern.fig');
%exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/pattern.eps','ContentType','Vector','BackgroundColor','None');


w=figure(3);
w.Position =[500 500 500 500];
m=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/lin_stability/unstable.mat');
subplot(5,1,1)
plot(x,m.Output.Freq{1,2}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -0.3, 0.3])
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)
title('$t=2$','Interpreter','latex','FontSize',15)

subplot(5,1,2)
plot(x,m.Output.Freq{1,100}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -0.3, 0.3])
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)
title('$t=100$','Interpreter','latex','FontSize',15)

subplot(5,1,3)
plot(x,m.Output.Freq{1,150}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -1, 1])
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)
title('$t=150$','Interpreter','latex','FontSize',15)

subplot(5,1,4)
plot(x,m.Output.Freq{1,250}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -1, 1])
set(gca,'xtick',[])
ylabel('\Omega_A','FontSize',15)
title('$t=250$','Interpreter','latex','FontSize',15)

subplot(5,1,5)
plot(x,m.Output.Freq{1,1000}{1,1},'k-','LineWidth',2);
axis([-Lx/2, Lx/2, -1, 1])
ylabel('\Omega_A','FontSize',15)
xlabel('x','FontSize',20)
title('$t=1000$','Interpreter','latex','FontSize',15)

saveas(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/unstable.fig');
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/unstable.eps','ContentType','Vector','BackgroundColor','None');
