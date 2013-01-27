% Example 4.16 -- Figure 4.6:  simple cubic polynomial interpolation

t = [0,0.1,0.8,1]'; b = [1,-0.9,10,9]';
A = zeros(4,4); % This tells Matlab we have a 4 by 4 matrix in mind
powers = 0:3; 
for j=1:4
  A(:,j) = t.^powers(j);
end

x = A \ b; % This solves the system Ax = b

tt = -0.1:.01:1.1; 
pt = x(1) + x(2).*tt + x(3).*tt.^2 + x(4).*tt.^3;
plot(tt,pt,'LineWidth',2) 
hold on 
plot(t',b','ro','LineWidth',2) 
xlabel('t')
ylabel('v')