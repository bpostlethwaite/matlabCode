%clear all
close all
C1 = load('/home/ben/workspace/programming/school/final/traffic.dat');

t3 = C1(:,1);
p3 = C1(:,2:end);
x = linspace(0,1000,length(p3(1,:)));

figure(1)
mesh(x,t3,p3)

figure(2)
plot(t3,p3(:,501))