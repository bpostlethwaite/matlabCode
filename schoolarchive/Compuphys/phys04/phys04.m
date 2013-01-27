close all

load('phys04.mat');

m = 0.1;
k = 6.4;
mk = sqrt(m/k);


for ii = 3:3
    figure(1)
    plot(RK2{ii}(:,1),log10(abs(RK2{ii}(:,2) - mk*cos(1/mk*RK2{ii}(:,1) + 3*pi/2))),'*r')
    figure(2)
    plot((RK4{ii}(:,1)),log10(abs((RK4{ii}(:,2) - mk*cos(1/mk*RK4{ii}(:,1) + 3*pi/2))./(mk*cos(1/mk*RK4{ii}(:,1) + 3*pi/2)))),'+b')
end
    title('Log10 Relative error for RK4 at h = 0.001')
    xlabel('Time [s]/\omegao')
    ylabel('Log 10 Relative Error (RK4 - 1/\omegacos(\omegat - \phi))/\omegacos(\omegat - \phi)')


figure(3)   %Log delta error versus time
for ii = 1:4
    hold on
    plot(RK2{ii}(:,1),log10(abs(RK2{ii}(:,5))),'*r',...
    RK4{ii}(:,1),log10(abs(RK4{ii}(:,5))),'+b')    
end
    title('Log \DeltaE versus time for RK2 & RK4 and 4 stepsizes')
    xlabel('Time [s]/\omegao')
    ylabel('Log \DeltaE')

    
figure(4)
    plot(RK4{3}(:,1),RK4{3}(:,2))
    title('Position as a function of time [t]/\omegao for RK4 at h = 0.001')
    xlabel('Time [s]/\omegao')
    ylabel('position')


figure(5)
    plot(RK4{3}(:,2),RK4{3}(:,3))
    axis square
    title('Phase Space Diagram')
    xlabel('velocity')
    ylabel('Position')
    xlim([-.2 .2])
    ylim([-.2 .2])
    grid on
