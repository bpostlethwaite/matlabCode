% generate the DC resistivity matrix
% on a uniform mesh using nodal finite difference
%  Solve on different mesh and see how the error decays
%

clear all


%for i=3:8
    i = 4;
    n = 2^i;
    h = 1/n;

    % staggered mesh
    [x,y] = ndgrid(0:h:1,0:h:1);                 % nodal (for u)
    [xc,yc] = ndgrid(h/2:h:1-h/2,h/2:h:1-h/2);   % cell center (for sigma)
    [xe1,ye1] = ndgrid(0:h:1,h/2:h:1-h/2);       % y face for ux
    [xe2,ye2] = ndgrid(h/2:h:1-h/2,0:h:1);       % x face for uy

    
    % Conductivity function
    sigma = cos(2*pi*xc).^2 .* cos(2*pi*yc).^2 + 1;
    %sigma = ones(n,n);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % derivative matrix in 1D
    % deivatives on nodes
    en   = ones(n+1,1);
    ddxn = @(m,k) 1/k*spdiags([-en en],[0,1],m,m+1); 
    % deivatives on cells
    ec   = ones(n,1);
    ddxc = @(m,k) 1/k*spdiags([-ec ec],[-1,0],m+1,m); 
    % Average matrix in 1D
    e   = ones(n,1);
    av = @(m) spdiags([e e]/2,[-1,0],m+1,m); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2D
    I = speye(n+1);
    Ic = speye(n);
    Dn  = ddxn(n,h);
    Dc  = ddxc(n,h);
    V  = av(n);

    % correct for BC
    V([1,end])  = 1;
    Dc([1,end]) = [2; -2]/h;


    GRAD = [kron(I,Dn); kron(Dn,I)];
    DIV  = [kron(I,Dc), kron(Dc,I)];
    Avrg = [kron(Ic,V); kron(V,Ic)];

    SigmaMid = Avrg * sigma(:);
    ne       = size(GRAD,1);
    A = DIV* spdiags(SigmaMid,0,ne,ne) * GRAD;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % testing
    u    = cos(2*pi*x) .* cos(2*pi*y);
    ux   = -2*pi*sin(2*pi*xe2(:)) .* cos(2*pi*ye2(:));
    uy   = -2*pi*sin(2*pi*ye1(:)) .* cos(2*pi*xe1(:));
    uxx = -4*pi^2*cos(2*pi*x) .* cos(2*pi*y);
    uyy = -4*pi^2*cos(2*pi*y) .* cos(2*pi*x);
    q   = uxx+uyy;

    q = q(:);

    % error in the gradient equation
    r1 = GRAD*u(:) - [ux(:); uy(:)];
    % Error in the div equation
    r2 = DIV*[ux(:); uy(:)] - q(:);
    % total error in discretization
    r  = A*u(:)-q(:);

    fprintf('%3d     %3.2e     %3.2e     %3.2e\n',n, ...
              norm(r1,'inf'), norm(r2,'inf'),norm(r,'inf'));


    r = reshape(abs(r),n+1,n+1);
    imagesc(r); colorbar
    
    pause(1);
%end