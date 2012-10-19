% Figure 15.1 : plot areas under curve for different basic rules

% define f as cubic
clear all
clf
xi = [0,.2,.6,1]';
yi = [3,2,2,2]';
A = [xi.^3;xi.^2;xi.^1;xi.^0]; A = reshape(A,4,4);
c = A \ yi;

x = 0:.01:1;
f = polyval(c,x);

% exact integral
subplot(2,2,1)
fill([1,0,x],[0,0,f],'c')
hold on
xlabel('x')
axis([-.1 1.2 0 3])
plot(x,f)


subplot(2,2,3)
% trap
fill([0,1,1,0],[0,0,2,3],'c')
hold on
xlabel('x')
axis([-.1 1.2 0 3])
plot(x,f)

subplot(2,2,4)
% mid
ymid = polyval(c,.5);
fill([0,1,1,0],[0,0,ymid,ymid],'c')
hold on
xlabel('x')
axis([-.1 1.2 0 3])
plot(x,f)

subplot(2,2,2)
% simpson
xsi = [0,.5,1]';
ysi = [3,ymid,2]';
As = [xsi.^2;xsi.^1;xsi.^0]; As = reshape(As,3,3);
cs = As \ ysi;
fs = polyval(cs,x);
fill([1,0,x],[0,0,fs],'c')
hold on
xlabel('x')
axis([-.1 1.2 0 3])
plot(x,f)
