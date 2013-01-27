clear all
close all
% Define some parameters, time n such
t = 0:0.01:20;
f = 1; b = 3/40;
A(1,:) = sin(pi.*t/20).^2;
A(2,:) = [sin(pi*t(t<=2)/4).^2, ones(1,length(t(t>2 & t<18))), sin(pi*(t(t>=18)-20)/4).^2];
for ii = 1:2
vib = A(ii,:).*sin(2*pi.*(f + b*t).*t);
vib = [vib,zeros(1,length(vib)+1)];
newt= linspace(0,40,length(vib));
% Perform ffts, multiply in freq domain = convolution in time domain
vibf     = fft(vib);
autovibf = vibf.*conj(vibf);
autovib  = fftshift(ifft((autovibf)));
% Double autocorrelation
doubleautof = autovibf.*conj(autovibf);
doubleauto = fftshift(ifft(doubleautof));
% PLOTS
figure(ii)
subplot(3,1,1)    
    plot(newt,vib)
        xlabel('time [s]')
        ylabel('Amplitude')
        title('Vibrosis signal (sweep function)')    
subplot(3,1,2)
    plot(newt,(autovib))
        xlim([18,22])
        xlabel('time [s]')
        ylabel('Amplitude')
        title('AutoCorrelated Vibrosis signal (sweep function)')
subplot(3,1,3)
    plot(newt,doubleauto)
        xlim([18,22])
        xlabel('time [s]')
        ylabel('Amplitude')
        title('Double AutoCorrelated Vibrosis signal (sweep function)')
end    