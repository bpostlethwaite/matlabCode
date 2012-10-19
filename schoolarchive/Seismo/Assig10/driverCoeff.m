clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Parameters
t = 0:0.025:10-0.025;
omega = 2*pi;
beta = [3900,4490];
theta = [0:5:85]';
p = sind(theta)/beta(1);
rho  = [2900,3380];
N = length(t);
ff = [0:N/2,-N/2+1:-1]/(N*0.025);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ricker Wave and Angle
ft = (1-omega^2.*(t-5).^2/2).*exp(-omega^2.*(t-5).^2/4);
cos2 = sqrt(1-p.^2*beta(2)^2);

% Compute Coefficients
for ii = 1:length(theta)
    top    = (rho(1)*beta*cosd(theta(ii)) - rho(2)*beta(2)*cos2(ii));
    bottom = (rho(1)*beta*cosd(theta(ii)) + rho(2)*beta(2)*cos2(ii));
    SS(ii) = top/bottom;
end

% Multiply Real Part
for ii = 1:13
    ftC(ii,:) = SS(ii) .* ft;
end

% Perform fft, ifft and multiply Imaginary section
FT = fft(ft);

for ii = 14 : length(SS)
    FTC(ii,1:N/2) = FT(1:N/2) * (SS(ii));
    FTC(ii,N/2+1:length(FT)) = FT(N/2+1:length(FT)) * conj((SS(ii)));
    ftC(ii,:) = real(ifft(FTC(ii,:)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params for part b
radius  = 3480.0;
alpha   = [13.72,8.06];
beta    = [7.26,0.001];
rho     = [5.57,9.90];
p = 1/alpha(1)*sind(theta);

% Get Coefficients
for ii = 1 : length(p)
    [RdPP(ii),RdSP(ii),RdPS(ii),RdSS(ii),TuPP(ii),TuSP(ii),TuPS(ii),TuSS(ii)] =...
        rtcoef(alpha(1),beta(1),rho(1),alpha(2),beta(2),rho(2),p(ii));
end
CoeffR = [RdPP;RdPS;RdSS;TuPP;TuSS];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS
figure(1)
subplot(2,1,1)
    plot(t,ftC(1:13,:))
    xlim([4,6])
    legend([repmat('theta = ',13,1) num2str(theta(1:13))],'Location','BestOutside')
    ylabel('Amplitude')
    xlabel('\theta')
    title('Ricker Wave at Real Phase')
subplot(2,1,2)
    plot(t,ftC(14:end,:))
    xlim([4,6])
    legend([repmat('theta = ',5,1) num2str(theta(14:end))],'Location','BestOutside')
    ylabel('Amplitude')
    xlabel('\theta')
    title('Ricker Wave at Imaginary Phase')
    
figure(2)
    plot(theta,CoeffR)
    legend('RdPP','RdPS','RdSS','TuPP','TuSS','Location','Best')
    title('Amplitude & Phase of downgoing P Waves')
    ylabel('Relative Amplitude')
    xlabel('Theta [degrees]')

figure(3)
    section(ftC,2,1,-1)


