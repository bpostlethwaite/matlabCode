% Example 16.5 : primitive forward Euler for stiff problem

h = [.0005,.001]*pi;
N = 1./[.0005,.001]; N = N/2;

y0 = 1; 

y = y0; t = 0;
for i=1:N(1)
    y = y + h(1)*(-1000*(y - cos (t)) - sin (t));
    t = t + h(1);
end
y

y = y0; t = 0;
for i=1:N(2)
    y = y + h(2)*(-1000*(y - cos (t)) - sin (t));
    t = t + h(2);
end
y