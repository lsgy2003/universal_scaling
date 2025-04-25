clear all;

sigma=0.005;

ini=load('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/uni_ini.mat');
jp_all=ini.jp_all;
jp1=jp_all(jp_all<=0.0097);
jp2=jp_all(jp_all>=0.0097);
jp3=0.0097;

T=5000;
T1=10000;
%T1=10000;
dt = 0.001;
Evo=1000;
Nt = T/dt;
tw = [100 400 1000];
%tw = [1000 4000 8000];
Dt=1:dt*Evo:(T-tw(3))+1;
Dt1=1:dt*Evo:(T1-tw(3))+1;
Lx=2^12;

clf;
%
width_s=zeros(1,length(jp_all));


for i=1:39
    m=load(sprintf('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/Size/jp=%.4f_tw=1000.mat',jp_all(i)));
    corr=m.corr(:);
    Dcorr = -log(corr);
    avg_corr = mean(Dcorr(2000:end));
    width_s(i) = avg_corr;
    loglog(Dt,Dcorr);
    hold on;
end 

for i=40:47
    m=load(sprintf('/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/Size/jp=%.4f_tw=1000.mat',jp_all(i)));
    corr=m.corr(:);
    Dcorr = -log(corr);
    avg_corr = mean(Dcorr(2000:end)); %calculate the averaged correlation after the saturation
    width_s(i) = avg_corr;
    loglog(Dt1,Dcorr);
end
hold off;
xlabel('t','FontSize',20)
ylabel('\chi_{AA}(t)','FontSize',20)
%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_sigma=0.005_1.fig';
saveas(gcf,ImageID);

filename='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Turing/sigma=0.005/ws_all.mat';
save(filename,'Lx','width_s','jp_all','sigma'); 
%}
%
figure(2)
plot(jp1,width_s(1:17)/sigma^2,'-square','MarkerFaceColor',[255 29 37]/255,'MarkerSize',8,'LineWidth',1,'Color','k');
hold on;
plot(jp2,width_s(17:end)/sigma^2,'-square','MarkerFaceColor',[0 113 188]/255,'MarkerSize',8,'LineWidth',1,'Color','k');
hold on;
plot(jp3,width_s(17)/sigma^2,'o','MarkerFaceColor',[57 181 74]/255,'MarkerSize',10,'LineWidth',1,'Color','k');
hold on;
plot(jp2(end),width_s(end)/sigma^2,'^','MarkerFaceColor',[255 191 0]/255,'MarkerSize',12,'LineWidth',1,'Color','k');
hold off;

set(gca,'yscale','log')
xlabel('j_+')
ylabel('\chi_{AA}^s/\sigma^2')
ax = gca;
ax.FontSize=18;
axis([0 0.041 10^(1) 5*10^(3)]);
txt={'$\sigma=0.005$'};
text(0.03,2.5*10^(3),txt,'Interpreter','latex','FontSize',20,'FontWeight','bold');
%
ImageID='/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_sigma=0.005_2.fig';
saveas(gcf,ImageID);
exportgraphics(gcf,'/Users/Phantom/Documents/MATLAB/Flocking/Density/New_data/Figure/xover_sigma=0.005_2.eps','ContentType','vector','BackgroundColor','none');
