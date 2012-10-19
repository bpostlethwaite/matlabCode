% Question 1c
clear all
close all


N = 2^7 - 1;
h = 1/(N+1);
for w = [1,10]
    x = linspace(0,2*pi,N);
    [X,Y] = meshgrid(x,x);
    g = h^2*ones(N^2,1);
    
    G = numgrid('S',N+2);
    D = delsq(G) - spdiags(w^2*ones(N^2,1),0,N^2,N^2);
    [L,U] = luinc(D,'0');
    [u,flag,relres,iter] = minres(D,g,1e-6,3000);
    %[u,flag,relres,iter] = gmres(D,g,20,1e-6,1000,L,U);
    fprintf('w of %i finds residual norm of %2.3g in %i iterations\n',w,relres,iter)
end
