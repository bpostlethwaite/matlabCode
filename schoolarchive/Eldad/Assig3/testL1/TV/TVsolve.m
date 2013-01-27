function[x] = TVsolve(A,b,regpar,n,maxiter,x0)
% [x] = TVsolve(A,b,regpar,n)
%

% generate gradients
G       = getGradientMatrix(n);
Av     = getAverageMatrix(n);
eps   = 1e-5;
TV    = @(z) (sum(sqrt( Av*((G*z).^2) + eps)));
dTV  = @(z) G'* sdiag(Av'*(1./sqrt( Av*((G*z).^2) + eps) ))*G*z;
HTV  = @(z) G'* sdiag(Av'*(1./sqrt( Av*((G*z).^2) + eps) ))*G;


% initialize
if ~exist('x0'), x = A'*b/norm(A*ones(size(A,2),1)); else x=x0; end
gamma = 100;
iter = 1;
while 1,
    
    J   = 0.5*norm(A*x-b)^2 + regpar*TV(x);
    dJ = A'*(A*x-b) + regpar*dTV(x);
    H  = HTV(x);
    fprintf('%3d.0   %3.2e      %3.2e       %3.2e\n',iter,J,norm(A*x-b), TV(x));
    Hess = @(z) (A'*(A*z)) + regpar*H*z;
    M       = @(z) (H + gamma*speye(prod(n)))\z;
    s =  pcg(Hess,-dJ,1e-2,5,M);
    
    muls = 1; lsiter = 1;
    if max(abs(s)) > 1, s = s/max(abs(s)); end
    while 1,
       xt  = x + muls*s;
       Jt   = 0.5*norm(A*xt-b)^2 + regpar*TV(xt);
       fprintf('%3d.%1d   %3.2e      %3.2e       %3.2e\n',iter,lsiter,Jt,norm(A*xt-b), TV(xt));
       
       if Jt < J, break; end
       muls = muls/2; lsiter = lsiter+1;
       if lsiter>10, disp('LSB'); return; end
    end
    x = xt;
    iter = iter+1;
    if iter > maxiter, disp('max iter achieved'); return; end
    if norm(s,'inf') < 1e-6, disp('norm(s)<1e-3'); return; end
end
        
    


return

function[G] = getGradientMatrix(n)

e1 = ones(n(1),1); e2 = ones(n(2),1);
h1 = 1/n(1); h2 = 1/n(2);

G1 = 1/h1*spdiags([-e1  e1],[0,1],n(1)-1,n(1));
G2 = 1/h2*spdiags([-e2  e2],[0,1],n(2)-1,n(2));

G = [kron(speye(n(2)),G1); kron(G2,speye(n(1)))];

return

function[Av] = getAverageMatrix(n)

e1 = ones(n(1),1); e2 = ones(n(2),1);

A1 = 0.5*spdiags([e1  e1],[-1,0],n(1),n(1)-1); A1([1,end]) = 1;
A2 = 0.5*spdiags([e2  e2],[-1,0],n(2),n(2)-1); A2([1,end]) = 1;

Av = [kron(speye(n(2)),A1), kron(A2,speye(n(1)))];

return

function[T] = sdiag(t)

T = spdiags(t(:),0,numel(t),numel(t));

