% Example 10.4 -- Figure 10.5 
clear all
clf

%evaluation mesh
x = .2:.01:5.2;

%quadrtic interpolant
xi = [1,2,4]; yi = [1,3,3];
[coef,table] = divdif(xi,yi);

%evaluate quadratic
y2 = evalnewt(x,xi,coef);

%add data point
xi = [xi,5]; yi = [yi,4];

%construct and evaluate cubic interpolant
[coef,table] = divdifadd(xi,yi,table);
y3 = evalnewt(x,xi,coef);

plot (x,y2,'b',x,y3,'g')
hold on
xlabel('x')
ylabel('p')
gtext('p_2(x)')
gtext('p_3(x)')