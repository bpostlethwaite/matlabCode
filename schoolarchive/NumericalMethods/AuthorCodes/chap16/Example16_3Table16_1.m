% Example 16.3 -- Table 16.1
% y' = y, y(0) = 1; primitive program

h = .1; t = 0:h:.6;
ye = exp(t);
y(1) = ye(1);
for i=1:6
    y(i+1) = y(i) + h*y(i);
end
ye
y
err1 = abs(y-ye)

clear all

h = .2; t = 0:h:.6;
ye = exp(t);
y(1) = ye(1);
for i=1:3
    y(i+1) = y(i) + h*y(i);
end
y
err2 = abs(y-ye)