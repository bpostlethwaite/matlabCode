% Example 16.6 : primitive backward Euler for stiff problem
% y' = -1000(y - cos t) - sin t, y(0) = 1


y0 = 1; t0 = 0;
nn = [1000,500];

% backward Euler
for j = 1:2
    h = pi/2/nn(j); y = y0; t = t0;
    for i = 1:nn(j)
        t = t + h;
        y = (y + h*(1000*cos(t) - sin(t)))/(1+h*1000);
    end
    y
end


