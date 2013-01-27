function out = simulation(DIV,A1,A2,A3,q,u,m)

a1 = A1*(1./m);
a2 = A2*(1./m);
a3 = A3*(1./m);

Ainv1 = spdiags(1./a1,0,size(A1,1),size(A1,1));
Ainv2 = spdiags(1./a2,0,size(A2,1),size(A2,1));
Ainv3 = spdiags(1./a3,0,size(A3,1),size(A3,1));

Sinv = blkdiag(Ainv1,Ainv2,Ainv3);

out = DIV*Sinv*DIV'*u - q;
