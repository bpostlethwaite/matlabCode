% generate the DC resistivity matrix
% on a uniform mesh using nodal finite difference
%  Solve on different mesh and see how the error decays
%

clear all
close all

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 	%Testing
    i = 4;
    n(1) = 2^i;
    n(2) = 2*n(1);
    h(1) = 1/n(1);
    h(2) = 1/n(2);
    scale = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     scale = 1200;
%     n(1) = 300;
%     n(2) = 2*n(1);
%     h(1) = scale/n(1);
%     h(2) = scale/n(2);
    
    % staggered mesh
    [x,y] = ndgrid(0:h(1):scale,0:h(2):scale);                 % nodal (for u)
    [xc,yc] = ndgrid(h(1)/2:h(1):scale-h(1)/2,h(2)/2:h(2):scale-h(2)/2);   % cell center (for sigma)
    [xe1,ye1] = ndgrid(0:h(1):scale,h(2)/2:h(2):scale-h(2)/2);       % y face for ux
    [xe2,ye2] = ndgrid(h(1)/2:h(1):scale-h(1)/2,0:h(2):scale);       % x face for uy

    % Conductivity function
    %sigma = cos(2*pi*xc).^2 .* cos(2*pi*yc).^2 + 1;
    
    sigma = ones(n(1),n(2))*0.1;
    sigma0 = sigma;
    %sigma(25:50,288:313) = 1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % derivative matrix in 1D
    % derivatives on nodes
    ddxn = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[0,1],m,m+1); 
    % derivatives on cells
    ddxc = @(m,k) 1/k*spdiags([-ones(m,1) ones(m,1)],[-1,0],m+1,m); 
    % Average matrix in 1D
    av = @(m) spdiags([ones(m,1) ones(m,1)]/2,[-1,0],m+1,m); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2D
    Ix = speye(n(1)+1);
    Icx = speye(n(1));
    Dnx  = ddxn(n(1),h(1));
    Dcx  = ddxc(n(1),h(1));
    Vx  = av(n(1));
    Iy = speye(n(2)+1);
    Icy = speye(n(2));
    Dny  = ddxn(n(2),h(2));
    Dcy  = ddxc(n(2),h(2));
    Vy  = av(n(2));

    % correct for BC
    Vx([1,end])  = 1;
    Dcx([1,end]) = [2; -2]/h(1);
    Vy([1,end])  = 1;
    Dcy([1,end]) = [2; -2]/h(2);

    % Use opposite I so scales match up
    GRAD = [kron(Iy,Dnx); kron(Dny,Ix)];
    DIV  = [kron(Iy,Dcx), kron(Dcy,Ix)];
    Avrg = [kron(Icy,Vx); kron(Vy,Icx)];

    SigmaMid  = Avrg * sigma(:);
    SigmaMid0 = Avrg * sigma0(:);
    ne        = size(GRAD,1);
    A  = DIV * spdiags(SigmaMid,0,ne,ne) * GRAD;
    A0 = DIV * spdiags(SigmaMid0,0,ne,ne) * GRAD;
    % Set a point to pin down solutions (make non-singular)
    A(end,end) = A(end,end) - h(1)*h(2);
    A0(end,end) = A0(end,end) - h(1)*h(2);
    
    % Set Dipole Dipole charges
    % for various spacings
    spacing = [50,100,500,1000];
    centre = round(0.5*n(2));
    
    

    
%     for ii = 1:length(spacing)
%         cdist = round(0.5*spacing(ii)/h(2));
%         q = sparse(n(1)+1,n(2)+1);
%         q(1,centre-cdist) = 1/(h(1)*(h(2)));
%         q(1,centre+cdist) = -1/(h(1)*(h(2)));
%         %Now solve Ax=b for x, x=A\b and reshape u
%         U  = A\q(:);
%         U0 = A0\q(:);
%         u  = reshape(U,n(1)+1,n(2)+1);
%         u0 = reshape(U0,n(1)+1,n(2)+1);
%         % Get d and ratio of d
%         d  = 1/h(2) * diff(u(1,:));
%         d0 = 1/h(2) * diff(u0(1,:));
%         figure(5)
%         subplot(4,1,ii)
%             plot(y(1,2:end),d./d0)
%             title(sprintf('Ratio Data for spacing %i metres',spacing(ii)))
%             xlabel('Metres')
%             ylabel('d/do (ratio)')
%             
%     end
    
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
    r = reshape(abs(r),n(1)+1,n(2)+1);
    imagesc(r); colorbar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    
    
%end




