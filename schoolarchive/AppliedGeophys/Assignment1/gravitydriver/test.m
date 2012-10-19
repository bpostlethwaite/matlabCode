% n=10;
% e = ones(n,1);
% dx = spdiags([-e,e],[0,1],n-1,n);
% I = speye(n);
% beta = 1e-2;
% W = [kron(I,kron(I,dx)); kron(I,kron(dx,I)); kron(dx,kron(I,I));...
%     beta*kron(I,kron(I,I))];

% clear all; 
% kappa =2;
% gamma = 10;
% n = 9;
% g= linspace(0,100,n+1)';
% for ii = 1: n +1
%     z((ii-1)*100 + 1 :  ii*100,1) = g(ii); 
% end
% fz =10./(z.^kappa + gamma);
% 
% fz2 = ((z-50).^2)./3000;
% 
% Wz = diag(fz);
% 
% figure(1)
% subplot(1,2,1)
% plot(z,fz);
% subplot(1,2,2)
% plot(z,fz2)
% 


figure(4)
plot(alphaDat{1},alphaDat{2},alphaDat{1},alphaDat{1}.*alphaDat{3})
