clear all
close all

load data.mat

x1 = data(1:5,:)'; x1 = x1(:);
x2 = data(6:10,:)'; x2 = x2(:);
x3 = data(11:15,:)'; x3 = x3(:);
x4 = data(16:20,:)'; x4 = x4(:);
x5 = data(21:25,:)'; x5 = x5(:);
x6 = data(26:30,:)'; x6 = x6(:);
y  = data(31:35,:)'; y = y(:);

X = [ones(length(x1),1),x1,x2,x3,x4,x5,x6];

Beta = (X'*X)\X'*y;

yy = X*Beta;

plot(x1,y,'b*',x2,y,'k*',x3,y,'r*',x4,y,'y*',x5,y,'g*',x6,y,'m*')
hold on

stepwise(X(:,2:end),y)

%Rank X6 X4 X3 X1 X2 X5


