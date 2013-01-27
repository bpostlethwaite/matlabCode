function [W,Wz] = getW(n,depthTarget)

e = ones(n,1);
dx = spdiags([-e,e],[0,1],n-1,n);
I = speye(n);
beta = 1e-2;
W = [kron(I,kron(I,dx)); kron(I,kron(dx,I)); kron(dx,kron(I,I));...
    beta*kron(I,kron(I,I))];

kappa =2;
gamma = 10;
g= linspace(0,100,n)';
for ii = 1: n
<<<<<<< HEAD
    z((ii-1)*100 + 1 :  ii*100,1) = g(ii); 
end
%fz =10./(z.^kappa + gamma);
%Wz = diag(fz);
fz2 = ((z-depthTarget).^2)./3000;
=======
    z( (ii-1) * 100 + 1 :  ii*100,1) = g(ii); 
end
%fz =10./(z.^kappa + gamma);
%Wz = diag(fz);
fz2 = ((z - depthTarget).^2)./3000;
>>>>>>> origin/master
Wz = diag(fz2);
end