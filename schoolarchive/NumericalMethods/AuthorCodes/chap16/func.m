function f = func(t,y)
% functions to integrate using rk4

isw = 5;
if isw == 1
    % predator - prey
    a = .25; b = -.01; c = -1; d = .01;
    f(1) = a*y(1) + b*y(1)*y(2);
    f(2) = c*y(2) + d*y(1)*y(2);
    %f = shiftdim(f);
    
elseif isw == 2
    % astronoical example
    mu = .012277471; muh = 1 - mu;
    D1 = ((y(1) + mu)^2 + y(3)^2)^(3/2);
    D2 = ((y(1) - muh)^2 + y(3)^2)^(3/2);
    f(1) = y(2);
    f(2) = y(1) + 2*y(4) - muh*(y(1) + mu)/D1 - mu*(y(1)-muh)/D2;
    f(3) = y(4);
    f(4) = y(3) - 2*y(2) - muh*y(3)/D1 - mu*y(3)/D2;
    
elseif isw == 3
    % BVP example
    f(1) = y(2);
    f(2) = -exp(y(1));
    
elseif isw == 4
    
    %f = 2*t* sqrt(max(1 - y^2,0));
    f(1) = y(2);
    f(2) = (1/t)*y(2) - 4*(t^2)*y(1);
    
elseif isw == 5
    
    % function for tables 16.2 - 16.4
    f = -y(1)^2;
    
end

f = shiftdim(f);