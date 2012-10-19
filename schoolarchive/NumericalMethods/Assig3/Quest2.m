clear all
close all
clc

tol = 1e-10;
nmax = 100;
ii = 1;
h = 2^-7;
N = 1/h - 1;
G = numgrid('S',N);
D = delsq(G);
x = 0:1/(N-1):1;
solver = char(['backslash';'min   res']);
for jj = 1:2
    %u = ones(N-2,N-2);
    %u = buildPyramid(N-2,1);
    u = 3*rand(N-2,N-2);
    u(1,:) = 0; u(:,1) = 0; u(end,:) = 0; u(:,end) = 0;
    u = u(:);
    
    for k = 1:nmax
        [F,J] = buildJacoII(u,D,N,h);
        if jj == 1
            p = -J \ F;
        else
            p = minres(-J,F,1e-7,1000);
        end
        u = u + p;
        if norm(p) < tol*(1+norm(u))
            F = buildJacoII(u,D,N,h);
            break
        end
    end

fprintf('\n\nfor %s solver find scaled 2 norm of %1.4d\n',solver(jj,:),norm(u)/N)
fprintf('and exponential infinite norm is %1.4d\n',norm(exp(u),inf))
fprintf('taking %i iterations for a solution\n\n',k)

Q = G;
Q(G>0) = u;
%U = reshape(u,N,N);
figure(jj)
surf(x,x,Q)
%surf(x,x,U)
shading interp
end

