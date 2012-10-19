function g = func14(c)
% algebraic function for simple shooting example
[tot,yot] = ode45(@func, [0 1], [0 c]);
g = yot(end,1);