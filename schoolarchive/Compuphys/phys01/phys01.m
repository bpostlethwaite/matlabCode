%A=getdatas('phys01');
%load('physvar');
E = 1.1929093e-7;
close all
X = log10(A(:,1));
Y = log10((abs(A(:,2)-A(:,4))./A(:,4)));
Z=log10( (abs(A(:,3)-A(:,4))./A(:,4)));

AA=[X , X.^0];
my=(AA'*AA)\AA'*Y;
mz=(AA'*AA)\AA'*Z;

y1=10.^(my(1).*X+my(2));
y2=10.^(mz(1).*X+mz(2));

randomwalk =  sqrt(A(:,1)).*E;
correlated =A(:,1).*E;



figure(1)
loglog(A(:,1), (abs(A(:,2)-A(:,4))./A(:,4)),A(:,1), (abs(A(:,3)-A(:,4))./A(:,4)),A(:,1),y1,'r--',A(:,1),y2,'r--')
xlabel('N');  ylabel('relative error');
title('LogLog plot of relative error')
legend('Sum 1 relative error','Sum 2 relative error','corrisponding linear fit','Location','Best')

figure(2)
loglog(A(:,1),y1,'b',A(:,1),y2,'k',A(:,1),randomwalk,'g--',A(:,1),correlated,'r--');
xlabel('N');  ylabel('relative error');
title('LogLog plot of relative error vs Random walk and Correlated machine errors')
legend('Sum 1 relative error Linear Fit','Sum 2 relative error Linear Fit','Random Walk Error','Correlated Error','Location','Best')

% figure(2)
% plot(log(A(:,1)), (abs(A(:,2)-A(:,4))./A(:,4))) 
% xlabel('N'); ylabel('relative error');
% title('LogLog plot of relative error')
% legend('Sum 1 relative error','Location','Best')
% 
% figure(3)
% plot(log(A(:,1)), (abs(A(:,3)-A(:,4))./A(:,4)),'g--')
% legend('Sum 2 relative error','Location','Best')
% xlabel('N'); ylabel('relative error');
% title('LogLog plot of relative error')

%watch this

