% Vectored reflection and transmission
% Plots vectorized reflection and transmission coeffs


clear all
close all
% Define Parameters
plotswitch = 1;
rad1  = pi/180*(90:-1:0)';
rho1  = 2.9;
beta1 = 3.9;
rho2  = 3.38;
beta2 = 4.49;
u     = 1/beta1;
p     = u*sin(rad1);
cos2  = sqrt(1-p.^2.*beta2^2);

% Define Coeffiecients and Vectors
Ivec  = [sin(rad1),-cos(rad1)];
R     = (rho1*beta1*cos(rad1) - rho2*beta2*cos2)./(rho1*beta1*cos(rad1) + rho2*beta2*cos2);
Rvec  = [sin(rad1),cos(rad1)];
T     = (2*rho1*beta1*cos(rad1))./(rho1*beta1*cos(rad1) + rho2*beta2*cos2);
Tvec  = [(sqrt(T.^2 - (cos2.*T).^2)),(-cos2.*T)];

% Plots
if plotswitch
for ii = 1:length(rad1)
    figure(1)
    H=line([-.1,0.32,0.32],[0.25,0.25,0.6]);
        set(H,'Color',[0 0 0],'LineWidth',2)
    hold on
    quiver(0,.6,Ivec(ii,1),Ivec(ii,2),0.3)
        xlim([-.1,2])
        ylim([-.8,.6])
    hold on
    quiver(0,0,Rvec(ii,1),Rvec(ii,2),abs(R(ii)))
    hold on
    
    if isreal(Tvec(ii,1))
        quiver(0,0,Tvec(ii,1),Tvec(ii,2),abs(T(ii)))
    else
        quiver(0,0,real(Tvec(ii,1)),0,abs(T(ii)))
        
    end
    
    pause(0.5)
end
end