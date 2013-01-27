% Example 16.23 -- Figure 6.14 :  Shooting for simple BVP example
% u'' + e^u = 0, u(0) = u(1) = 0.
% Thus, y_1' = y_2, y_2' = -exp(y_1), y_1(0) = 0.
% For y_2(0) = c solve ivode and match y_1(1) = 0.

c0 = input(' guess missing initial value for y_2 : ');

% solve IVP for plotting purposes
tspan = [0 1];
y0 = [0 ; c0];
[tot0,yot0] = ode45(@func, tspan, y0);

% solve BVP using fzero to solve g(c) = y_1(1;c) = 0. Note no use of g-derivative
c = fzero(@func14,c0);

% plot obtained BVP solution
y0 = [0 ; c];
[tot,yot] = ode45(@func, tspan, y0);
plot(tot0,yot0(:,1)','m--',tot,yot(:,1)','b')
xlabel('t')
ylabel('y_1')
legend('initial trajectory','BVP trajectory')