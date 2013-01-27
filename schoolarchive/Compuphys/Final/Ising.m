clear all
close all
%A = 1Temp, 2number of steps,  3accepted,  4mean energy, 5var E, 6mean
%magnetization, 7var M , 8abs M, 9susceptibility, 10Heat Capacity
load('Ising1')
p = 9;
a = 1;

G{1} = load('/home/ben/workspace/programming/data/IsingLat1.dat');
G{2} = load('/home/ben/workspace/programming/data/IsingLat8.dat');
G{3} = load('/home/ben/workspace/programming/data/IsingLat75.dat');
G{4} = load('/home/ben/workspace/programming/data/IsingLat7.dat');
G{5} = load('/home/ben/workspace/programming/data/IsingLat4.dat');
G{6} = load('/home/ben/workspace/programming/data/IsingLat0.dat');


T1 = A(:,1);
M1 = A(:,p);
T2 = B(:,1);
M2 = B(:,p);
T3 = C(:,1);
M3 = C(:,p);


% plot(T1,M1,'+',T2,M2,'v',T3,M3,'*')
% legend('10x10 Grid','30x30 Grid','50x50 Grid','Location','Best')
% xlabel('Heat Capacity')
% title('Heat Capacity')
% ylabel('Temperature')


Jh = 1;
Jv = [1,0.8,0.75,0.7,0.4,0];

sk = Jv/Jh;

for ii = 1:6
T(ii,:) = G{ii}(:,1);
M(ii,:) = G{ii}(:,p);
end
figure(2)
for ii = 1:6
subplot(3,3,ii)
plot(T(ii,:),M(ii,:),'+')
title(sprintf('Jv/Jh = %2.2f',sk(ii)))
end

%sprintf('G%i',skull(ii))(:,1)