% Examp16.4 -- Figure 16.3 : HIRES

isw = 1;
u = [1,0,0,0,0,0,0,.0057]';  % initial data

if isw == 1 % use forward Euler
  t = 0:.001:322; h = .001;
  y = u * ones(1,length(t));
  for i = 2:length(t)
    y(:,i) = y(:,i-1) + h*hires(t(i-1),y(:,i-1));
  end
  figure(isw)
  plot(t,y(6,:))
  xlabel('t')
  ylabel('y_6')
    
elseif isw == 2 % use Matlab's ode45
  tspan = 0:1:322;
  [t,y] = ode45(@hires,tspan,u);
  figure(isw)
  plot(t,y)
end   

 