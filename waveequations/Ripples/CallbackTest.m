clear all; close all
%set(gca,'ButtonDownFcn',@myCallback);
%then with a press of some button I callback with...
%pt = get(gca,'CurrentPoint');
%x = round(pt(1,1))
%y = round(pt(1,2))
%z = rfh_cb = @newfig; % Create function handle for newfig function

%figure('ButtonDownFcn',fh_cb);

x = linspace(0,4*pi,1000);
k = 1;
w = 1;
dt = 0.01;
t = 0;
f = sin(k*x - w*t);
fig = figure(1);
%h = plot(x,f,'YDataSource','f');
h = plot(f,'YDataSource','f');
while 1
    t = t + dt;
    f = sin(k*x - w*t);
	refreshdata(h,'caller') % Evaluate y in the function workspace
	drawnow; 
    point = get(gca,'CurrentPoint');
    title(sprintf('Current Point is %i %i',point(1),point(2)))
end

