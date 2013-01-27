kk = 1;
for n = [8,16,32,64,128,256,512]
    h = 1/n;
    % nodal grid
    tN = linspace(0,1,n+1); tN = tN(:);
    % cell centered grid
    tC = tN(1:end-1) + diff(tN)/2; tC = tC(:);
    % Boundary points
    tB = [tN(1); tN(end)];
    % define the fields and fluxes
    a = 1e0; b = pi/2+0.01;
    u = @(t)(atan(a*t - 1/2*a).*t - 1/2*atan(a*t - 1/2*a) ...
        -1/(2*a)*log((a*t - 1/2*a).^2 + 1) + b*t);
    up = @(t)(atan(a*(t-0.5)) + b);
    m = @(t)(1./up(t).*exp(t));
    J = @(t)(up(t).*m(t));
    rhs = @(t)(exp(t));
    % compute operators
    e = ones(n+1,1);
    Dn2c = 1/h*spdiags([-e,e],[0,1],n,n+1);
    Dc2n = 1/h*spdiags([-e,e],[-1,0],n+1,n);
    Ac2n = 1/2*spdiags([e,e],[-1,0],n+1,n);
    s = Ac2n*(1./m(tC));
    % correct for conductivity outside of the grid
    % (for second order at the boundary)
    s(1) = (1/m(h/2) + 1/m(-h/2))/4;
    s(end) = (1/m(1-h/2) + 1/m(1+h/2))/4;
    S = sdiag(1./s);
    % Boundary condition matrix
    B = sparse(n+1,2); B(1,1) = 1/h; B(end,end) = -1/h;
    % check truncation error
    r1 = Dc2n*u(tC) - S*J(tN) - B*u(tB);
    r2 = Dn2c * J(tN) - rhs(tC);
    % now sol0ve and compare solution error
    A = Dn2c*S*Dc2n;
    b = rhs(tC) + (Dn2c*S)*(B*u(tB));
    uN = A\b; uA = u(tC);
    r3 = uA-uN;
    % print some info
    rho(kk,1:4) = [h,norm(r1,'inf'),norm(r2,'inf'),norm(r3,'inf')]; kk = kk+1;
    fprintf(['%3.2e %3.2e %3.2e %3.2e\n',...
        h,norm(r1,'inf'),norm(r2,'inf'),norm(r3,'inf')]);
end
% check convergence rate
fprintf('\n\n\n Convergence summary\n\n');
fprintf('h |r1|_inf |r2|_inf |error|_inf \n\n')
fprintf('%3.2e %3.2e %3.2e %3.2e\n',(rho(1:end-1,:)./rho(2:end,:))')
