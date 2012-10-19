% Figure 10.1 : plot reasonable and unreasonable interpolations
% Do not look at this program before you have read through Section 10.4

clear all
close all
% data
xx = [1,2,4,5,6];
yy = [1,1.8,2,1.8,.5];

% poly interpolation based on data plus: "reasonable"
xxp = [xx,3]; yyp = [yy,2];
coef = divdif(xxp,yyp);

% evaluate on [.5,6.1] ]and plot
x = .5:.01:6.1;
y = evalnewt(x,xxp,coef);

figure(1)
plot(xx,yy,'bd',x,y,'g')
hold on
axis([0 7 0 2.5])
xlabel('x')
ylabel('v')

% poly interpolation based on data plus: "unreasonable"
xxp = [xx,1.5,2.9,4.4,5.3]; yyp = [yy,2.3,.4,2.4,1.7];
coef = divdif(xxp,yyp);

% evaluate on [.5,6.1] ]and plot
x = .5:.01:6.1;
y = evalnewt(x,xxp,coef);

figure(2)
plot(xx,yy,'bd',x,y,'g')
hold on
axis([0 7 0 2.5])
xlabel('x')
ylabel('v')