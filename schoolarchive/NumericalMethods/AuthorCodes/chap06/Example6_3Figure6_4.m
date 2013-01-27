% Example 6.3 : low degree polynomial data fitting

clear all
% data
m = 21;
tt = 0:1/(m-1):1;
bb = cos(2*pi*tt);

% find polynomial coefficients
for n=1:5
  coefs{n} = lsfit(tt,bb,n);
end

% Evaluate and plot
t = 0:.01:1;
z = ones(5,101);
for n=1:5
  z(n,:) = z(n,:) * coefs{n}(n);
  for j=n-1:-1:1
    z(n,:) = z(n,:).*t + coefs{n}(j);
  end
end
plot(t,z,tt,bb,'ro')
xlabel('t')
ylabel('p_{n-1}')