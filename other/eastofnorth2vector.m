%% eastofnorth2vector
clear all
close all
ang = 0:359;


x = sind(ang);
y = cosd(ang);

rads = ang*pi/180;


newang = 180 / pi * acos( ( x .*0 + y .* 1 ) ./ sqrt(x.^2 + y.^2) );
newang(x < 0) = 360 - newang(x < 0) ;


disp([ang(:),newang(:)])

%{

for ii = 1:length(x)
    figure(1)
    plot( x(ii), y(ii) )
    hold on
    axis square
    xlim([-1.2,1.2])
    ylim([-1.2,1.2])
    pause(0.01)
end

%}