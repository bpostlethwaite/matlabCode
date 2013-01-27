function [x,res] = poismg(A,b,x,level)
%
% function [x,res] = poismg(A,b,x,level)
%
% multigrid V-cycle to solve simplest Poisson on a square
% The uniform grid is N by N, N = 2^l-1 some l > 2,
% b is the right hand side; homogeneous Dirichlet;
% A has been created by   A = delsq(numgrid('S',N+2));

coarsest = 3;              % coarsest grid
nu1 = 2;                   % relaxations before coarsening grid
nu2 = 2;                   % relaxations after return to finer level
omeg = .8;                 % relaxation damping parameter

if level == coarsest
    x = A\b;               % solve exactly on coarsest level

else % begin multigrid cycle
    
    % relax using damped Jacobi
    Dv = diag(A);         % diagonal part of A as vector
    for i=1:nu1
      r = b - A*x;
      x = x + omeg*r./Dv;
    end
    
    % restrict residual from r to rc on coarser grid
    r = b - A*x; 
    N = sqrt(length(b));
    r = reshape(r,N,N);
    Nc = (N+1)/2 - 1; nc = Nc^2;    % coarser grid dimensions
    Ac = delsq(numgrid('S',Nc+2));  % coarser grid operator
    rc = r(2:2:N-1,2:2:N-1) + .5*(r(3:2:N,2:2:N-1)+r(1:2:N-2,2:2:N-1) +...
        r(2:2:N-1,3:2:N)+r(2:2:N-1,1:2:N-2)) + .25*(r(3:2:N,3:2:N)+...
        r(3:2:N,1:2:N-2)+r(1:2:N-2,3:2:N)+r(1:2:N-2,1:2:N-2));
    rc = reshape(rc,nc,1);
    
    % descend level. Use V-cycle
    vc = zeros(size(rc));            % initialize correction to 0
    [vc,r] = poismg(Ac,rc,vc,level-1); % samesame on coarser grid
    
    % prolongate correction from vc to v on finer grid
    v = zeros(N,N);
    vc = reshape(vc,Nc,Nc);
    v(2:2:N-1,2:2:N-1) = vc;
    vz = [zeros(1,N);v;zeros(1,N)];   % embed v with a ring of 0s
    vz = [zeros(N+2,1),vz,zeros(N+2,1)];
    v(1:2:N,2:2:N-1) = .5*(vz(1:2:N,3:2:N)+vz(3:2:N+2,3:2:N));
    v(2:2:N-1,1:2:N) = .5*(vz(3:2:N,1:2:N)+vz(3:2:N,3:2:N+2));
    v(3:2:N-2,3:2:N-2) = .25*(v(2:2:N-3,2:2:N-3)+v(2:2:N-3,4:2:N-1)+...
        v(4:2:N-1,4:2:N-1)+v(4:2:N-1,2:2:N-3));
    
    % add to current solution
    n = N^2;
    x = x + reshape(v,n,1);
    
    % relax using damped Jacobi
    for i=1:nu2
      r = b - A*x;
      x = x + omeg*r./Dv;
    end
     
end
res = norm(b - A*x);
    
