% Figure 3.4

x=-1:.001:1;
x0=0.5;
x1=x0-(sin(2*x0)+0.1)/(2*cos(2*x0));
x2=x1-(sin(2*x1)+0.1)/(2*cos(2*x1));

figure
plot(x,sin(2*x)+0.1)
axis ([-1 1 -1 1.5])
hold on
xlabel('x')

plot(linspace(-1,1),linspace(0,0),'k');
plot(linspace(x0,x1),linspace(sin(2*x0)+0.1,0),'g');
plot(linspace(x1,x2),linspace(sin(2*x1)+0.1,0),'g');

plot(linspace(x0,x0),linspace(0,sin(2*x0)+0.1),'r:');
plot(linspace(x1,x1),linspace(0,sin(2*x1)+0.1),'r:');

plot ([x0,x1,x2],[0,0,0],'+')

gtext('x_0')
gtext('x_1')
gtext('x_2')


  
 
  
 
 
