function [u,J] = dc3dfwd_multq(n1,n2,n3,m,q,Q,numq)

% function [u,J] = dc3dfwd(n1,n2,n3,m)
%   DC Resistivity 3D Forward Modeling Code
%   For an input conductivity m and sources q, with dimensions n1,n2,n3
%   the code calculates the potentials u and the sensitivity matrix J.


    %% For testing, use these initial values (expand comment section)

    %{
% Define dimensions
N = 8;
n1=N;n2=N;n3=N;
n = n1*n2*n3;
h1 = 1/n1;
h2 = 1/n2;
h3 = 1/n3;

% Define q
numq = 2;
qtot = zeros(n,numq);
for i = 1:numq;
    q = zeros(n1,n2,n3);
    q(1,5,2) = -1/(h1*h2*h3);
    q(1,5,i+3) = 1/(h1*h2*h3);
    q = q(:);
    qtot(:,i) = q;    
end
q = qtot(:);

% Define Q
num_surf = n2*n3;
Q = speye(num_surf);
k = zeros(1,n1);k(1)=1;
Q = kron(Q,k);
Q = kron(speye(numq),Q);

% Define true model
m = ones(n,1);
%m = kron(ones(numq,1),m);
%}

    %% Assign dimensions to variables

    n = n1*n2*n3;
    m = m(:);
    q = q(:);
    h1 = 1/n1;
    h2 = 1/n2;
    h3 = 1/n3;
        
    %% Form derivative matrices (center to node, and node to center)
    
    % Cell centres to nodes:
    Dc2n_1 = 1/h1*spdiags([-ones(n1,1),ones(n1,1)],-1:0,n1+1,n1);
    Dc2n_2 = 1/h2*spdiags([-ones(n2,1),ones(n2,1)],-1:0,n2+1,n2);
    Dc2n_3 = 1/h3*spdiags([-ones(n3,1),ones(n3,1)],-1:0,n3+1,n3);
    
    % Nodes to cell centres
    Dn2c_1 = 1/h1*spdiags([-ones(n1,1),ones(n1,1)],0:1,n1,n1+1);
    Dn2c_2 = 1/h2*spdiags([-ones(n2,1),ones(n2,1)],0:1,n2,n2+1);
    Dn2c_3 = 1/h3*spdiags([-ones(n3,1),ones(n3,1)],0:1,n3,n3+1);
    
    %% Boundary Conditions for GRAD operator - Dirichlet
    
    % For testing, place 2's in active direction, 0 in inactive
    Dc2n_1(1,1) = 2*Dc2n_1(1,1);
    Dc2n_1(end,end) = 2*Dc2n_1(end,end);
    Dc2n_2(1,1) = 2*Dc2n_2(1,1);
    Dc2n_2(end,end) = 2*Dc2n_2(end,end);
    Dc2n_3(1,1) = 2*Dc2n_3(1,1);
    Dc2n_3(end,end) = 2*Dc2n_3(end,end);
    
    %% Create 3D DIV operator using Kroenicker products.
    
    D1 = kron(speye(n3),kron(speye(n2),Dn2c_1));
    D2 = kron(speye(n3),kron(Dn2c_2,speye(n1)));
    D3 = kron(Dn2c_3,kron(speye(n2),speye(n1)));
    DIV = [D1 D2 D3];
    
    %% Create 3D GRAD operator using Kroenicker products
    
    G1 = kron(speye(n3),kron(speye(n2),Dc2n_1));
    G2 = kron(speye(n3),kron(Dc2n_2,speye(n1)));
    G3 = kron(Dc2n_3,kron(speye(n2),speye(n1)));
    GRAD = [G1;G2;G3];
    
    %% Create Averaging Matrix
    
    Ac2n1 = spdiags(0.5*ones(n1+1,2),-1:0,n1+1,n1);
    Ac2n1([1,end]) = 1*[1,1];
    Ac2n2 = spdiags(0.5*ones(n2+1,2),-1:0,n2+1,n2);
    Ac2n2([1,end]) = 1*[1,1];
    Ac2n3 = spdiags(0.5*ones(n3+1,2),-1:0,n3+1,n3);
    Ac2n3([1,end]) = 1*[1,1];
    
    %% Form S, m_2, m_in2 and Averaging operators
    
    Acf1 = kron(speye(n3),kron(speye(n2),Ac2n1));
    Acf2 = kron(speye(n3),kron(Ac2n2,speye(n1)));
    Acf3 = kron(Ac2n3,kron(speye(n2),speye(n1)));
    S = @(m)((blkdiag(diag(Acf1*(1./m)),diag(Acf2*(1./m)),diag(Acf3*(1./m)))));
    
    %% Form S_in and S_in2 operators
    
    [sx,sy] = size(S(m));
    S_in1 = @(m)((1./diag(S(m))));
    S_in2_prep = @(m)(1./ (diag(S(m)).*diag(S(m))) );
    S_in = @(m)(spdiags((S_in1(m)),0,sx,sy));
    S_in2 = @(m)(spdiags((S_in2_prep(m)),0,sx,sy));
    m_2 = @(m)((1./(m.^2)));
    m_in2 = @(m)(spdiags([m_2(m)],0,n,n));
    Av = [Acf1;Acf2;Acf3];
    
    %% Form Boundary matrix - Derichlet B.C.
    
    B_1 = zeros(n1+1,2); B_1(1,1) = 2/h1; B_1(end,end) = -2/h1;
    B_2 = zeros(n2+1,2); B_2(1,1) = 2/h2; B_2(end,end) = -2/h2;
    B_3 = zeros(n3+1,2); B_3(1,1) = 2/h3; B_3(end,end) = -2/h3;
    Bx = kron(speye(n3),kron(speye(n2),B_1));
    By = kron(speye(n3),kron(B_2,speye(n1)));
    Bz = kron(B_3,kron(speye(n2),speye(n1)));
    ub = kron(ones(n3,1),kron(ones(n2,1),[0;0]));
    Bc = blkdiag(Bx,By,Bz);
    Bcmq = kron(ones(1,numq),Bc);
    ubc = [ub;ub;ub];
    ubcmq = kron(ones(numq,1),ubc);
    
    %% Calculate u_numerical
    
    % Define A Matrix
    A = DIV*S_in(m)*GRAD;
    Amq = kron(speye(numq),A);
    
    % Calculate u
    u = Amq \ (q(:));%+DIVmq*S_in(m)*Bcmq*ubcmq);
    u = u(:);
    
    %%  View potentials u
  
    %{
    u1 = u(1:512);
    u2 = u(513:end);
    u = reshape(u1,n1,n2,n3);
    uplot = squeeze(u(1,:,:));
    imagesc(uplot)
    %}
    
    %% Calculate Sensitivities:
    
    % Form operators for Sensitivities
    DIVmq = kron(speye(numq),DIV);
    GRADmq = kron(speye(numq),GRAD);
    DIV_u = DIVmq'*u;
    ndiv_u = length(DIV_u);
    GRAD_u = @(u)(spdiags(GRADmq*u,0,ndiv_u,ndiv_u));
    S_inmq = @(m)(kron(speye(numq),S_in(m)));
    S_in2mq = @(m)(kron(speye(numq),S_in2(m)));
    Avmq = kron(ones(numq,1),Av);
    
    % Form functions for derivatives
    C = @(u,m,ubc)(DIVmq*S_inmq(m)*GRADmq*u - q);% + DIVmq*S_in(m)*Bcmq*ubcmq);
    dCdu = @(u,m)(DIVmq*S_inmq(m)*GRADmq);
    dCdm = @(u,m)(DIVmq*GRAD_u(u)*S_in2mq(m)*Avmq*m_in2(m));
    %dCdu_bc = @(u,m)(DIV*S_inmq(m)*Bc);
    f = C(u,m,ubc);
    
    % Calculate Sensitivity Matrix
    % J = Q*dudm or J = -Q*(dCdu)^-1*(dCdm)  
    J = -Q*(dCdu(u,m) \ dCdm(u,m));
    
    %% View Sensitivities
    
    %{
    J = -(dCdu(u,m)\dCdm(u,m));
    J_1 = full(J(1:n,1)); % To see first set of sensitivities
    J_2 = full(J(n+1:end,1)); % To see second set
    J1 = reshape(J_1,n1,n2,n3);
    for i = 1:n1
        uplot = squeeze(J1(i,:,:));
        figure(1)
        subplot(3,3,i)
        imagesc(uplot)
        xlabel('x');ylabel('y')
        title(['depth-slice :',num2str(i)])
        colorbar
        pause(1)
    end
    %}
    
    %% Derivative tests for Sensitivities
        
    %{
    % Derivative test for dCdu: passed!
    v = randn(n*numq,1);
    fprintf('Derivative test for dCdu \n');
    fprintf('h          diff1      diff2\n')
    for i = 1:10
        h = 10^(-i);
        fp = C(u+h*v,m,ubc);
        diff1 = norm(fp - f);
        diff2 = norm(fp - f - h*dCdu(u,m)*v);
        fprintf('%3.2e  %3.2e %3.2e \n',h,diff1,diff2)
    end
    
    % Derivative test for dCdm: passed!
    v1 = randn(n,1);
    fprintf('Derivative test for dCdm \n');
    fprintf('h          diff1      diff2\n')
    for i = 1:10
        h = 10^(-i);
        fp = C(u,m+h*v1,ubc);
        diff1 = norm(fp - f);
        diff2 = norm(fp - f - h*dCdm(u,m)*v1);
        fprintf('%3.2e  %3.2e %3.2e \n',h,diff1,diff2)
    end
    
    % Derivative test for dCdu_bc:  don't need since ubc = 0
%     v = randn(length(ubc),1);
%     fprintf('Derivative test for dCdu_bc \n');
%     fprintf('h          diff1      diff2\n')
%     for i = 1:10
%         h = 10^(-i);
%         fp = C(u,m,ubc+h*v);
%         diff1 = norm(fp - f);
%         diff2 = norm(fp - f - h*dCdu_bc(u,m)*v);
%         fprintf('%3.2e  %3.2e %3.2e \n',h,diff1,diff2)
%     end
    %}  
    
end

