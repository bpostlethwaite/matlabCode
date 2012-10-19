% Figure 9_5 : draw the figure of the book cover
    
x = -1:.01:3.1;
y = 0:.01:1.1;
[X,Y] = meshgrid(x,y);
phi = .5 * ((1.5 - X.*(1-Y)).^2 + (2.25 - X.*(1-Y.^2)).^2 + ...
      (2.625 - X.*(1-Y.^3)).^2 );
mesh(x,y,phi)
xlabel('x_1')
ylabel('x_2')
zlabel('\phi')