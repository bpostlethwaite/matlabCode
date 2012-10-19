% Example 2.8 -- Figure 2.3

x = [];

% Generate all positive numbers of the system (2,3,-2,3)
for i = 1:0.25:1.75,
    for j = -2:3,
        x = [x i*2^j];
    end
end

x=[x -x 0];     % Add all negative numbers and 0
x = sort(x);    % Sort
y = zeros(1,length(x));
plot(x,y,'+')
