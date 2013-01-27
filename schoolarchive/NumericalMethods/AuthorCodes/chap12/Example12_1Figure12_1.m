% Example 12.1 - Figure 12.1 : best polynomial approximation to cos(2 pi x)

clear all
close all
% construct Hilbert matrix for n up to 4 (i.e. size 5)
n = 4; np1 = n+1;
B  = zeros(np1,np1);
for j=1:np1
    for k=1:np1
        B(j,k) = 1/(j+k-1);
    end
end

b = [0, 0, 1/(2*pi^2), 3/(4*pi^2), 1/pi^2- 3/(2*pi^4)]';

x = 0:.01:1; % evaluation mesh
u = cos(2*pi*x); % exact
figure(1)
clf
plot(x,u,'g')
hold on
xlabel('x')
ylabel('p_n')

xx = ones(length(x),1);
% find and plot approximations
for l = 1:np1 
    c = B(1:l,1:l) \ b(1:l);
    v = xx*c;
    if l == 2||l == 4||l == 5, plot(x,v), end
    xx = [xx,xx(:,l).*x'];
end

legend('exact','n=0 or 1','n=2 or 3','n=4')
