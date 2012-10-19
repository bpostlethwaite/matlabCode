close all

% MARMOD, parameters.
depth=[0.0,1.5,6.0,6.5,10.0];
alpha=[4.5,6.8,7.0,8.0,8.1];
beta=[2.4,3.75,3.85,4.6,4.7];
rho=[2.0,2.8,2.9,3.1,3.1];
p=linspace(0.1236,0.2217,100)';

[Xp,Tp] = getdxdt(depth,alpha,p);
Xp = Xp*2; Tp=Tp*2;

tor = Tp-Xp.*p;


figure(1)
    subplot(3,1,1)
    plot(Xp,Tp-Xp/8);
    title('T(X) curve ''opened up'' with reduction velocity x/v','FontSize',16)
    xlabel('X(p)','FontSize',16)
    ylabel('T(p)','FontSize',16)

    subplot(3,1,2)
    plot(p,Xp)
    title('X(p) curve','FontSize',16)
    xlabel('p','FontSize',16)
    ylabel('X(p)','FontSize',16)
    
    subplot(3,1,3)
    plot(p,tor)
    title('\tau(p) Delay Time Function','FontSize',16)
    xlabel('p','FontSize',16)
    ylabel('\tau(p)','FontSize',16)