close all
clear all
load prem.mat
radius(end) = 10;
earthR = 6.371*1e3;
p=linspace(0.001,0.1128,201)';
depthf = -earthR*log(radius./earthR);
alphaf = (earthR./radius).*alpha;

[Xp,Tp] = getdxdt(depthf,alphaf,p);
Xp = 2*Xp; Tp=2*Tp;
Xdeg=(180/(pi*earthR)).*Xp;
tord = Tp-Xdeg.*p;


figure(1)
    plot(alphaf,depthf)
    set(gca, 'YDir', 'reverse')
    xlabel('velocity','FontSize',16)
    ylabel('Depth','FontSize',16)
    title('Depth versus Velocity Curve for Flat Earth transformation','FontSize',16)
    ylim([0,12000])
    xlim([0,70])

figure(2)
    plot(Xdeg,Tp/60,'*');
    title('T(X) curve showing shadow zone','FontSize',16)
    xlabel('X(p) [degrees]','FontSize',16)
    ylabel('T(p) [minutes]','FontSize',16)
    %ylim([0,25])
    xlim([0,180])

figure(3)
subplot(2,2,1)
    plot(Xdeg,Tp-Xdeg./0.1,'-');
    title('T(X) curve zoomed in with Redution Velocity','FontSize',16)
    xlabel('X(p) [degrees]','FontSize',16)
    ylabel('T(p)-X(p)/Vr [seconds]','FontSize',16)
    ylim([60,90])
    xlim([12,35])
subplot(2,2,3)    
    plot(Xdeg,p,'-');
    title('X(p) versus p','FontSize',16)
    xlabel('X(p) [degrees]','FontSize',16)
    ylabel('p','FontSize',16)
    %ylim([40,90])
    xlim([10,35])
subplot(2,2,[2,4]) 
    plot(p,tord,'*');
    title('\tau versus p','FontSize',16)
    xlabel('p','FontSize',16)
    ylabel('\tau','FontSize',16)
    %ylim([40,90])
    xlim([p(50),p(end)])
    
    
    figure(10)
close 10