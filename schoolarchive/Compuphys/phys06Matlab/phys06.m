clear all
close all
%

% Test1 = load('/home/ben/workspace/programming/data/ConstantHeatTest1.dat');
% tTest1 = Test1(:,1);
% Test1 = Test1(:,2:end);
% 
% xTest1 = linspace(0,1,length(Test1(1,:)));
% 
% for ii = 1:length(Test1(:,1))
%     
%  figure(1)
%     hold on
%     plot(xTest1,Test1(ii,:))
% end

% Melting of silver = 1235

A1 = load('/home/ben/workspace/programming/data/propane_1_1.dat');
tA = A1(:,1);
LA = A1(:,2:5);
A1 = A1(:,6:end);
xA = linspace(1,30,length(A1(1,:)));

B1 = load('/home/ben/workspace/programming/data/propane_1_2.dat');
tB = B1(:,1);
LB = B1(:,2:5);
B1 = B1(:,6:end);
xB = linspace(0,30,length(B1(1,:)));

%for ii = 1:length(A1(:,1))
    
B1 = A1;
LB = LA;

% figure(1)
% plot(xA,A1')
%     legend([repmat('t = ',length(tA),1) num2str(tA)],'Location','BestOutside')
%     xlim([1,20])
%     title('Propane Torch on Silver Sheet')
%     xlabel('Radius [mm]')
%     ylabel('Temp [K]')


figure(2)
    ttB = tB(1:10:end);
    plot(xB,B1(1:10:end,:)')
    legend([repmat('t = ',length(ttB),1) num2str(ttB)],'Location','BestOutside')
    xlim([1,20])
    title('Propane Torch on Silver Sheet')
    xlabel('Radius [mm]')
    ylabel('Temp [K]')

figure(3)
hold on
plot(tA,LA,'v')
ylim([0,2.5]);
title('Latent Heat of Melting versus Time')
xlabel('Time [s]')
ylabel('Joules')
legend('r = 1mm','r = 2mm','r = 3mm' , 'r = 4mm','Location','Best')

% x1 = 1.74*99.4718;
% x2 = 0.8*x1;
% x3 = 0.6*x1;
% x4 = 0.4*x1;
% Heat = pi*(x1 + 3*x2 + 5*x3 + 7*x4);
% W = [x1,x2,x3,x4];
% S = W/0.429;
% W = 99.4718;
% S = W/0.429
% 
Circ = 2*pi*[13,14,15,16,17];
TempForce = [300,310,320,310,300] - 293;
TempForce = sum(TempForce .* Circ);

DiskTemp = TempForce / (30^2 * pi) + 293;
SteadyStateA = mean(A1,2);
SteadyStateB = mean(B1,2);

% figure(3)
% plot(tA,log10(abs(SteadyStateA - DiskTemp)/DiskTemp),'r+',...
%     tB,log10(abs(SteadyStateB - DiskTemp)/DiskTemp),'bv')
% legend('Explicit (Forward) Relative Error','Implicit (Backward) Relative Error')
% title('Mean Temp Relative Error for times for given initial temp Distribution')
% xlabel('Time [s]')
% ylabel('log error')

